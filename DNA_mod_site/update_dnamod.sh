#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

SOURCE_DIR="$(dirname $(readlink -f $0))" # From: https://gist.github.com/tvlooy/cbfbdb111a4ebad8b93e

# Updates DNAmod database and generates static site
echo "Begining to update DNAmod and create static site"

# remove any existing database to fully re-create it
rm -f $SOURCE_DIR/../DNA_mod_database.db

$SOURCE_DIR/../populate_database_sql.py
echo "Database updated"

$SOURCE_DIR/create_mod_staticsite_sql.py

echo "Static site created"

