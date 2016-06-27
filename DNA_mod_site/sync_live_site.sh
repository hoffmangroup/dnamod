#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

ERR_EXIT=64

# Copy contents of ./static/* into directory for live site

if [[ ! ($(hostname) == 'mordor' && -n $(groups | fgrep 'hoffmangroup')) ]]; then
    >&2 echo "Must be run on the Hoffman Lab cluster, by a lab member."
    exit $ERR_EXIT
fi

SOURCE_DIR="$(dirname $(readlink -f $0))" # From: https://gist.github.com/tvlooy/cbfbdb111a4ebad8b93e

# Testing mode: If this is set to true the site will be copied into a test folder in ~/DNA_Base_Database/DNA_mod_site
TESTING_MODE=false

MAIN_SITE_DIR="$SOURCE_DIR/static"

# Path variables for testing mode and real mode
TEST_PATH="$SOURCE_DIR/sync_test"
REAL_PATH='/mnt/work1/users/hoffmangroup/www/proj/dnamod'

if [[ "$TESTING_MODE" == true ]]; then
    rm -Rf "$TEST_PATH"
    COPY_PATH="$TEST_PATH"    
    mkdir "$COPY_PATH"
else
    COPY_PATH="$REAL_PATH"
fi

chmod -Rv a+rX "$MAIN_SITE_DIR/.."
rsync --progress -av "$MAIN_SITE_DIR"/* "$COPY_PATH"

