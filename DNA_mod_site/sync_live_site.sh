#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

# Copy contents of ./static/* into directory for live site

# Testing mode: If this is set to true the site will be copied into a test folder in ~/DNA_Base_Database/DNA_mod_site
TESTING_MODE=false

MAIN_SITE_DIR='/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/static'
# Path variables for testing mode and real mode
TEST_PATH='/mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site/sync_test'
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

