"""
MIT License

Copyright (c) 2026 Dominik Kopczynski  -  dominik.kopczynski {at} univie.ac.at
                   Cristina Coman  -  cristina.coman {at} univie.ac.at
                   Nils Hoffmann  -  n.hoffmann {at} fz-juelich.de
                   Robert Ahrends  -  robert.ahrends {at} univie.ac.at

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

"""
Session state that has to survive running several MAPtoGO replicas behind a
load balancer.

* The secret key that signs the Flask session and encrypts the ``uid`` cookie
  is read from ``MAPTOGO_SECRET_KEY_FILE`` (e.g. a Docker secret) or
  ``MAPTOGO_SECRET_KEY``, so every replica can read every other replica's
  cookies. Without either, a random per-process key is used (standalone mode).
* The user's inputs (``SessionEntry.ui`` and ``use_bounded_fatty_acyls``) are
  stored in Redis when ``MAPTOGO_REDIS_URL`` is set, so a request that lands on
  another replica can restore them. Without it, nothing is stored outside the
  process (standalone mode). Analysis results stay in memory: they reference
  the in-memory ontology and are recomputed when needed.
"""

import base64
import json
import logging
import os
import secrets
import time

from cryptography.fernet import Fernet
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC

logger = logging.getLogger(__name__)

SECRET_KEY_ENV = "MAPTOGO_SECRET_KEY"
SECRET_KEY_FILE_ENV = "MAPTOGO_SECRET_KEY_FILE"
REDIS_URL_ENV = "MAPTOGO_REDIS_URL"
REDIS_KEY_PREFIX = "maptogo:session:"
REDIS_TIMEOUT_SECONDS = 2


# --- secret key ----------------------------------------------------------------

def load_secret_key(environ = None):
    environ = os.environ if environ is None else environ

    key_file = environ.get(SECRET_KEY_FILE_ENV)
    if key_file:
        with open(key_file, encoding = "utf-8") as f:
            key = f.read().strip()
        if not key:
            raise ValueError(f"Secret key file '{key_file}' from {SECRET_KEY_FILE_ENV} is empty")
        return key

    key = environ.get(SECRET_KEY_ENV, "").strip()
    if key:
        return key

    logger.warning(
        f"Neither {SECRET_KEY_FILE_ENV} nor {SECRET_KEY_ENV} is set, using a random secret key. "
        "Sessions are lost on restart and cannot be shared between replicas."
    )
    return secrets.token_urlsafe(64)


# --- uid cookie encryption ----------------------------------------------------------

def derive_key(password: str, salt: bytes) -> bytes:
    kdf = PBKDF2HMAC(
        algorithm = hashes.SHA256(),
        length = 32,
        salt = salt,
        iterations = 100_000,
        backend = default_backend()
    )
    return base64.urlsafe_b64encode(kdf.derive(password.encode()))


def encrypt_uuid(uuid_str: str, password: str) -> str:
    salt = os.urandom(16)  # random salt
    key = derive_key(password, salt)
    f = Fernet(key)
    token = f.encrypt(uuid_str.encode())
    return base64.urlsafe_b64encode(salt + token).decode()


def decrypt_uuid(encrypted_data: str, password: str) -> str:
    data = base64.urlsafe_b64decode(encrypted_data.encode())
    salt = data[:16]
    token = data[16:]
    key = derive_key(password, salt)
    f = Fernet(key)
    return f.decrypt(token).decode()


# --- session stores --------------------------------------------------------------

class InMemoryOnlyStore:
    """Standalone mode: sessions live only in the process' ``sessions`` dict."""

    def save(self, uid, state):
        pass

    def load(self, uid):
        return None

    def touch(self, uid):
        pass


class RedisSessionStore:
    """Stores session inputs as JSON in Redis, expiring after ``ttl_seconds``.

    Redis errors are logged and swallowed: if Redis is unavailable, each
    replica falls back to its in-memory sessions instead of failing requests.
    """

    def __init__(self, client, ttl_seconds, prefix = REDIS_KEY_PREFIX):
        self.client = client
        self.ttl_seconds = int(ttl_seconds)
        self.prefix = prefix

    def _key(self, uid):
        return f"{self.prefix}{uid}"

    def save(self, uid, state):
        try:
            payload = json.dumps(state)
        except (TypeError, ValueError) as e:
            logger.warning(f"Session {uid} could not be serialized for Redis: {e}")
            return
        try:
            self.client.set(self._key(uid), payload, ex = self.ttl_seconds)
        except Exception as e:
            logger.warning(f"Redis unavailable, session {uid} not saved: {e}")

    def load(self, uid):
        try:
            payload = self.client.get(self._key(uid))
        except Exception as e:
            logger.warning(f"Redis unavailable, session {uid} not loaded: {e}")
            return None
        if payload is None:
            return None
        try:
            return json.loads(payload)
        except (TypeError, ValueError) as e:
            logger.warning(f"Stored session {uid} is not valid JSON, ignoring it: {e}")
            return None

    def touch(self, uid):
        try:
            self.client.expire(self._key(uid), self.ttl_seconds)
        except Exception as e:
            logger.warning(f"Redis unavailable, session {uid} expiry not extended: {e}")


def create_session_store(ttl_seconds, environ = None):
    environ = os.environ if environ is None else environ
    redis_url = environ.get(REDIS_URL_ENV, "").strip()
    if not redis_url:
        logger.info(f"{REDIS_URL_ENV} is not set, keeping sessions in memory only")
        return InMemoryOnlyStore()

    import redis  # only needed when a Redis URL is configured

    client = redis.Redis.from_url(
        redis_url,
        socket_connect_timeout = REDIS_TIMEOUT_SECONDS,
        socket_timeout = REDIS_TIMEOUT_SECONDS,
    )
    logger.info(f"Storing session inputs in Redis at {client.connection_pool.connection_kwargs.get('host', redis_url)}")
    return RedisSessionStore(client, ttl_seconds)


# --- snapshot / restore -------------------------------------------------------------

def snapshot_session(session):
    return {
        "ui": dict(session.ui),
        "use_bounded_fatty_acyls": bool(session.use_bounded_fatty_acyls),
    }


def apply_snapshot(session, state):
    session.ui = dict(state.get("ui", {}))
    session.use_bounded_fatty_acyls = bool(state.get("use_bounded_fatty_acyls", False))


def get_session(sessions, store, uid, factory, create = False):
    """Return the session for ``uid`` from memory, else restore it from the
    store. If neither has it, create a new one when ``create`` is set,
    otherwise return None."""
    if not uid:
        return None

    session = sessions.get(uid)
    if session is not None:
        return session

    state = store.load(uid)
    if state is not None:
        session = factory()
        apply_snapshot(session, state)
        session.time = time.time()
        sessions[uid] = session
        logger.info(f"Restored session from shared store: {uid}")
        return session

    if create:
        session = factory()
        sessions[uid] = session
        return session

    return None
