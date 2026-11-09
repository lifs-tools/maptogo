#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

usage() {
  echo "Usage: $0 [build|up|down|logs|clean]"
  echo ""
  echo "  build  Build the development Docker image"
  echo "  up     Build (if needed) and start the dev server  [default]"
  echo "  down   Stop the dev server"
  echo "  logs   Follow container logs"
  echo "  clean  Stop and remove containers AND named volumes (forces full rebuild)"
  echo ""
  echo "The app will be available at http://localhost:8040"
}

case "${1:-up}" in
  build)
    docker compose -f docker-compose.dev.yml build
    ;;
  up)
    docker compose -f docker-compose.dev.yml up --build
    ;;
  down)
    docker compose -f docker-compose.dev.yml down
    ;;
  logs)
    docker compose -f docker-compose.dev.yml logs -f
    ;;
  clean)
    echo "Removing containers and named volumes — this will trigger a full rebuild on next 'up'."
    docker compose -f docker-compose.dev.yml down -v
    ;;
  help|--help|-h)
    usage
    ;;
  *)
    echo "Unknown command: $1"
    usage
    exit 1
    ;;
esac
