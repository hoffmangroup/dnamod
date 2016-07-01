#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

SOURCE_DIR="$(dirname $(readlink -f $0))" # From: https://gist.github.com/tvlooy/cbfbdb111a4ebad8b93e
CONSTANTS_SCRIPT="$SOURCE_DIR/constants.sh"

# Updates DNAmod database and generates static site
echo "Begining to update DNAmod and create static site..."

# remove any existing database to fully re-create it
rm -f $($CONSTANTS_SCRIPT 'database')

$($CONSTANTS_SCRIPT 'script_pop_db')
echo "Database updated."

$($CONSTANTS_SCRIPT 'script_create_site')

echo -e "\nSite creation complete.\n\nRun $($CONSTANTS_SCRIPT 'script_sync_site') to sync changes."

