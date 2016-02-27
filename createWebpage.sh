#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

rm -Rf DNA_mod_database.db

./populate_database_sql.py

DNA_mod_site/create_mod_staticsite_sql.py

