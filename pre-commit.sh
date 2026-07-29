#!/bin/bash
set -euo pipefail

ruff format
ruff check --fix

#echo "Unit tests"
#python -m unittest discover -s tests/taxsbp/unit/
#echo "Integration tests"
#python -m unittest discover -s tests/taxsbp/integration/
#pdoc -o docs taxsbp taxsbp.taxsbp taxsbp.utils