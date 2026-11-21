"""
Tests for SessionStore: shared secret key and optional Redis-backed session
state, so that several MAPtoGO replicas behind a load balancer can serve the
same user, while a standalone instance still runs without Redis.
"""

import logging
import time

import fakeredis
import pytest
import redis

from SessionStore import (
    InMemoryOnlyStore,
    RedisSessionStore,
    apply_snapshot,
    create_session_store,
    decrypt_uuid,
    encrypt_uuid,
    get_session,
    load_secret_key,
    snapshot_session,
)


class FakeSessionEntry:
    """Stand-in for analysis.SessionEntry with the fields the store persists."""

    def __init__(self):
        self.time = time.time()
        self.ui = {}
        self.use_bounded_fatty_acyls = False
        self.data_loaded = False


# --- secret key ---------------------------------------------------------------

def test_secret_key_from_file_wins_over_env(tmp_path):
    key_file = tmp_path / "secret"
    key_file.write_text("from-file\n")
    environ = {
        "MAPTOGO_SECRET_KEY_FILE": str(key_file),
        "MAPTOGO_SECRET_KEY": "from-env",
    }
    assert load_secret_key(environ) == "from-file"


def test_secret_key_from_env():
    assert load_secret_key({"MAPTOGO_SECRET_KEY": "from-env"}) == "from-env"


def test_secret_key_fallback_is_random_and_warns(caplog):
    with caplog.at_level(logging.WARNING):
        first = load_secret_key({})
        second = load_secret_key({})
    assert first and second and first != second
    assert "MAPTOGO_SECRET_KEY" in caplog.text


def test_empty_secret_key_file_is_rejected(tmp_path):
    key_file = tmp_path / "secret"
    key_file.write_text("  \n")
    with pytest.raises(ValueError):
        load_secret_key({"MAPTOGO_SECRET_KEY_FILE": str(key_file)})


def test_uid_cookie_from_one_replica_decrypts_on_another():
    environ = {"MAPTOGO_SECRET_KEY": "shared-secret"}
    replica_a_key, replica_b_key = load_secret_key(environ), load_secret_key(environ)
    cookie = encrypt_uuid("3544d671-95e7-4aa3-841e-81a56c616652", replica_a_key)
    assert decrypt_uuid(cookie, replica_b_key) == "3544d671-95e7-4aa3-841e-81a56c616652"


def test_uid_cookie_does_not_decrypt_with_other_key():
    cookie = encrypt_uuid("some-uid", "key-a")
    with pytest.raises(Exception):
        decrypt_uuid(cookie, "key-b")


# --- store selection ------------------------------------------------------------

def test_without_redis_url_the_store_is_in_memory_only():
    store = create_session_store(ttl_seconds=60, environ={})
    assert isinstance(store, InMemoryOnlyStore)
    store.save("uid", {"ui": {"a": 1}})  # must be a harmless no-op
    assert store.load("uid") is None


def test_with_redis_url_the_store_uses_redis():
    store = create_session_store(ttl_seconds=60, environ={"MAPTOGO_REDIS_URL": "redis://localhost:6379/0"})
    assert isinstance(store, RedisSessionStore)


# --- redis store ------------------------------------------------------------------

def test_redis_store_round_trip_with_ttl():
    client = fakeredis.FakeRedis()
    store = RedisSessionStore(client, ttl_seconds=7200)
    state = {"ui": {"all_lipids_list": "PC 34:1\nPE 36:2", "with_lipids": True}, "use_bounded_fatty_acyls": False}
    store.save("uid-1", state)
    assert store.load("uid-1") == state
    assert 0 < client.ttl("maptogo:session:uid-1") <= 7200
    assert store.load("unknown") is None


def test_redis_store_touch_extends_ttl():
    client = fakeredis.FakeRedis()
    store = RedisSessionStore(client, ttl_seconds=7200)
    store.save("uid-1", {"ui": {}})
    client.expire("maptogo:session:uid-1", 5)
    store.touch("uid-1")
    assert client.ttl("maptogo:session:uid-1") > 5


def test_redis_outage_does_not_break_requests(caplog):
    client = redis.Redis(host="127.0.0.1", port=1, socket_connect_timeout=0.2, socket_timeout=0.2)
    store = RedisSessionStore(client, ttl_seconds=60)
    with caplog.at_level(logging.WARNING):
        store.save("uid", {"ui": {}})
        assert store.load("uid") is None
        store.touch("uid")
    assert "Redis" in caplog.text


def test_unserializable_state_is_logged_not_raised(caplog):
    store = RedisSessionStore(fakeredis.FakeRedis(), ttl_seconds=60)
    with caplog.at_level(logging.WARNING):
        store.save("uid", {"ui": {"bad": object()}})
    assert store.load("uid") is None


# --- snapshot / restore -------------------------------------------------------------

def test_snapshot_and_apply_round_trip():
    session = FakeSessionEntry()
    session.ui = {"select_organism": "NCBITaxon:9606", "select_domains": ["Molecular function"]}
    session.use_bounded_fatty_acyls = True
    restored = FakeSessionEntry()
    apply_snapshot(restored, snapshot_session(session))
    assert restored.ui == session.ui
    assert restored.ui is not session.ui
    assert restored.use_bounded_fatty_acyls is True


# --- two replicas, one shared store -----------------------------------------------------

def test_submit_on_replica_a_is_found_by_replica_b():
    shared = RedisSessionStore(fakeredis.FakeRedis(), ttl_seconds=7200)
    sessions_a, sessions_b = {}, {}

    # replica A handles POST /submit
    session_a = get_session(sessions_a, shared, "uid-1", FakeSessionEntry, create = True)
    session_a.ui["upregulated_lipids_list"] = "PS 34:2"
    session_a.ui["select_organism"] = "NCBITaxon:9606"
    shared.save("uid-1", snapshot_session(session_a))

    # the browser then lands on replica B
    session_b = get_session(sessions_b, shared, "uid-1", FakeSessionEntry)
    assert session_b is not None
    assert session_b.ui == session_a.ui
    assert sessions_b["uid-1"] is session_b
    assert session_b.data_loaded is False


def test_local_session_wins_over_store():
    store = RedisSessionStore(fakeredis.FakeRedis(), ttl_seconds=60)
    store.save("uid-1", {"ui": {"select_organism": "from-store"}})
    local = FakeSessionEntry()
    local.ui["select_organism"] = "local"
    sessions = {"uid-1": local}
    assert get_session(sessions, store, "uid-1", FakeSessionEntry) is local


def test_unknown_session_returns_none_or_new_entry():
    store = InMemoryOnlyStore()
    sessions = {}
    assert get_session(sessions, store, "uid-1", FakeSessionEntry) is None
    assert sessions == {}
    created = get_session(sessions, store, "uid-1", FakeSessionEntry, create = True)
    assert isinstance(created, FakeSessionEntry)
    assert sessions["uid-1"] is created


def test_get_session_ignores_missing_uid():
    assert get_session({}, InMemoryOnlyStore(), None, FakeSessionEntry, create = True) is None
