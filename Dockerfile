FROM python:3.12-slim-trixie AS builder
ENV UV_COMPILE_BYTECODE=1 UV_LINK_MODE=copy
# Omit development dependencies
ENV UV_NO_DEV=1
# Do not install a new Python version, use the system / base image one
ENV UV_PYTHON_DOWNLOADS=0
# Install other dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    libxml2-dev \
    libxslt1-dev \
    && rm -rf /var/lib/apt/lists/* /usr/share/doc /usr/share/man \
    && apt-get clean

# Install uv
# Download the latest installer
ADD https://astral.sh/uv/0.9.20/install.sh /uv-installer.sh
# Run the installer then remove it
RUN sh /uv-installer.sh && rm /uv-installer.sh
# Ensure the installed binary is on the `PATH`
ENV PATH="/root/.local/bin/:$PATH"

# Create our working directory to host the App
RUN mkdir /wd
WORKDIR /wd
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project

# Copy project files
COPY . /wd
COPY uv.lock /wd
WORKDIR /wd/src
RUN make
WORKDIR /wd
# Set up dependencies with uv
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked

FROM python:3.12-slim-trixie
# development or production
ARG DASH_ENVIRONMENT="development"
RUN export ${DASH_ENVIRONMENT}

# Add non-root user
RUN useradd -ms /bin/sh -u 1001 app

# Copy the application from the builder
COPY --from=builder --chown=app:app /wd /wd

ENV PATH="/wd/.venv/bin:$PATH"

USER app
WORKDIR /wd
EXPOSE 8040
CMD [ "gunicorn", "--user", "app", "--group", "app", "--preload", "--timeout", "300", "--worker-class", "gevent", "--workers=1", "--threads=4", "-b 0.0.0.0:8040", "app:server"]
