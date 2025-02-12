#! /bin/bash
set -ex

if [[ -z $VIRTUAL_ENV ]]; then
    python3 -m venv venv
    source venv/bin/activate
fi

pip install -e ".[develop]"

pytest

ruff check
ruff format
