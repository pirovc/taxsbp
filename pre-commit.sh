#!/bin/bash
set -euo pipefail

ruff format
ruff check --fix


coverage erase
coverage run --source=taxsbp --omit="/usr/*,tests/*" -m pytest -s -vv
coverage report
coverage html