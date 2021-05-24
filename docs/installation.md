# Installation

Project is 'simple' repository, without setup script.
Installation is done by cloning repository and setting `PYTHONPATH`, and optionally `PATH`, environment variables.

If `ZCI_PROJECT_DIR` is direcotry of cloned project, than environment should be set to:
```
export PYTHONPATH=$ZCI_PROJECT_DIR:$PYTHONPATH
export PATH=$ZCI_PROJECT_DIR/bin:$PATH
```

## Caching

To use caching (common db), used for data fetching from Internet and some calculations, envirnoment variable
`ZCI_COMMON_DB` whould be set. E.g.
```
export ZCI_COMMON_DB=$HOME/.zcitools/common_db
```

## Other

Script `set_zci_env` in project's directory is an example of setting environment.
