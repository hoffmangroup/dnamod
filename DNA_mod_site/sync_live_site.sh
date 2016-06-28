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

# push to the internal www directory, to copy over to the external

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

>&2 echo "Fixing permissions and copying to www directory."

chmod -Rv a+rX "$MAIN_SITE_DIR/.."
rsync --progress -av "$MAIN_SITE_DIR"/* "$COPY_PATH"

# push the actual site to the external directory
if [[ "$TESTING_MODE" == false ]]; then
    # commit the changes to the lab Bitbucket
    >&2 echo "Committing changes."
    hg commit --config extensions.hgspellcheck=! -m "Updated DNAmod. Consult its repository ($($(hg paths default) | sed -r 's|ssh://.*?@|https://|')) for details." proj/dnamod

    >&2 echo "Pushing to www-external."
    EXTERNAL_DIR='/mnt/work1/users/hoffmangroup/www-external/proj/dnamod'
    rsync --progress -av "$REAL_PATH"/* "$EXTERNAL_DIR"

    >&2 echo -e "\n\nEmail Qun Jin <qjin@uhnresearch.ca> to push the public webpage.\n"
fi

