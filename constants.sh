#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

# File to store and retrieve global constants

# Usage: ./constants.sh <constrant name>

# Constant name		Description
# "root_dir"	    	The directory that is the root of the repository and location of this script
# "site_dir"		The subdirectory containing the website
# "site_html_dir"       The subdirectory containing the website's static HTML files
# "site_image_dir"     The subdirectory containing the website's images
# "site_css_dir"        The subdirectory containing the website's CSS files
# "site_template_dir"   The subdirectory containing the Jinja2 template files
# "data_dir"		The subdirectory containing associated data files or databases
# "database"		The SQLite database file
# "whitelist"           The list of nucleobases that are placed on the "Verified" list
# "blacklist"           The list of entries that are excluded, including their children (by default)
# "annot_exp_alph"      The annotations for the expanded epigenetic alphabet (Viner et al. 2016, bioRxiv, doi:10.1101/043794)
# "annot_seq"		The annotations with references to sequencing methods developed to elucidate various modified bases
# "annot_nature"	The annotations with references of the nature of modified nucleobases, including whether they are natural or only synthetic.
# "python_utils"	The Python module file
# "script_pop_db"	The script to create the SQLite database
# "script_create_site"	The script to create the webpages
# "script_update_all"	The script to re-create everything
# "script_sync_site"	The script to push the webpage changes to the external directory for synchronization

ERR_EXIT=64

function not_found {
    >&2 echo "Unable to find $1"
    exit $ERR_EXIT
}

SOURCE_DIR="$(dirname $(readlink -f $0))" # From: https://gist.github.com/tvlooy/cbfbdb111a4ebad8b93e

ROOT_DIR="$SOURCE_DIR"

if $(hash hg 2>/dev/null); then
    if [[ $(hg root) != "$ROOT_DIR" ]]; then
        >&2 echo "Mercurial root does not match the location of this script. Unable to correctly determine the root directory."
        exit $ERR_EXIT
    fi
else
    >&2 echo "Mercurial is not installed; skipping repo. root check."
fi

SITE_DIR="$SOURCE_DIR/DNA_mod_site"
SITE_TEMPLATE_DIR="$SITE_DIR/templates"
SITE_HTML_DIR="$SITE_DIR/static"
SITE_IMAGE_DIR="$SITE_HTML_DIR/images"
SITE_CSS_DIR="$SITE_HTML_DIR/css"
DATA_DIR="$SOURCE_DIR/data"

if [[ ! ( -d "$ROOT_DIR" && -d "$SITE_DIR" && -d "$SITE_TEMPLATE_DIR" && -d "$SITE_HTML_DIR" && -d "$SITE_IMAGE_DIR" && -d "$SITE_CSS_DIR" && -d "$DATA_DIR" ) ]]; then
    not_found 'directories'
fi

DATABASE="$DATA_DIR/DNA_mod_database.db"

WHITELIST="$DATA_DIR/whitelist.txt"
BLACKLIST="$DATA_DIR/blacklist.txt"

ANNOT_EXP_ALPH="$DATA_DIR/expanded_alphabet.txt"
ANNOT_SEQ="$DATA_DIR/ref_annots_sequencing.txt"
ANNOT_NATURE="$DATA_DIR/ref_annots_nature.txt"

# NB: do not check for the existence of the database,
#     since we often remove and re-create it

PYTHON_UTILS_MODULE="$ROOT_DIR/dnamod_utils.py"

if [[ ! -f "$PYTHON_UTILS_MODULE" ]]; then
    not_found "${PYTHON_UTILS_MODULE##*/}"
fi

SCRIPT_POP_DB="$ROOT_DIR/populate_database_sql.py"
SCRIPT_CREATE_SITE="$SITE_DIR/create_mod_staticsite_sql.py"
SCRIPT_UPDATE_ALL="$ROOT_DIR/update_dnamod.sh"
SCRIPT_SYNC_SITE="$SITE_DIR/sync_live_site.sh"

if [[ ! ( -x "$SCRIPT_POP_DB" && -x "$SCRIPT_CREATE_SITE" && -x "$SCRIPT_UPDATE_ALL" ) ]]; then
   not_found 'scripts' 
fi

case ${1:-} in 
    root_dir)
        echo -n "$ROOT_DIR"
        ;;
    site_dir)
        echo -n "$SITE_DIR"
        ;;
    site_html_dir)
        echo -n "$SITE_HTML_DIR"
        ;;
    site_image_dir)
        echo -n "$SITE_IMAGE_DIR"
        ;;
    site_csv_dir)
        echo -n "$SITE_CSS_DIR"
        ;;
    site_template_dir)
        echo -n "$SITE_TEMPLATE_DIR"
        ;;
    data_dir)
        echo -n "$DATA_DIR"
        ;;
    database)
        echo -n "$DATABASE"
        ;;
    whitelist)
        echo -n "$WHITELIST"
        ;;
    blacklist)
        echo -n "$BLACKLIST"
        ;;
    annot_exp_alph)
        echo -n "$ANNOT_EXP_ALPH"
        ;;
    annot_seq)
        echo -n "$ANNOT_SEQ"
        ;;
    annot_nature)
        echo -n "$ANNOT_NATURE"
        ;;
    python_utils)
        echo -n "$PYTHON_UTILS_MODULE"
        ;;
    script_pop_db)
        echo -n "$SCRIPT_POP_DB"
        ;;
    script_create_site)
        echo -n "$SCRIPT_CREATE_SITE"
        ;;
    script_update_all)
        echo -n "$SCRIPT_UPDATE_ALL"
        ;;
    script_sync_site)
        echo -n "$SCRIPT_SYNC_SITE"
        ;;
    *)
        # get all possible cases by searching this script
        >&2 echo "Usage: ./constants.sh <$(sed -n '/case/,/esac/p' $(readlink -f $0) | fgrep ')' | fgrep -v '*)' | fgrep -v 'sed' | cut -d ')' -f 1 | tr -d ' ' | paste -s -d '|')>"
        exit $ERR_EXIT
esac

