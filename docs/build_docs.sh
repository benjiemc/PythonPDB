#! /bin/bash

if [[ -z $VIRTUAL_ENV ]]; then
    python3 -m venv venv
    source venv/bin/activate
fi

pip install -e ".[develop]"

sphinx-apidoc -f -e -o docs/source src/python_pdb
cp examples/*.ipynb docs/source/

sphinx-build -b html ./docs public
