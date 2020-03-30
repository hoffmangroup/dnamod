#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

# File to store and retrieve global constants.
#
# -------------------------------------------------------------------------------
# Copyright (C) 2016  Coby Viner
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# -------------------------------------------------------------------------------

# Usage: ./constants.sh <constrant name>

# Constant name		Description
# "root_dir"	    	The directory that is the root of the repository and location of this script
# "site_dir"		The subdirectory containing the website
# "site_html_dir"       The subdirectory containing the website's static HTML files
# "site_image_dir"     The subdirectory containing the website's images
# "site_css_dir"        The subdirectory containing the website's CSS files
# "site_js_dir"        The subdirectory containing the website's Javascript files
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
# "json"		The subdirectory containing the JSON index file for Elasticlunr.js searching
# "seq_annot_table" A constant holding the name of the sequencing annotation table
# "nature_annot_table" A constant holding the name of the nature annotation table
# "references_table" A constant holding the name of the references table
# "exp_alpha_table" A constant holding the name of the expanded alphabet table
# "dnamod_version" A constant holding the current version number of DNAmod. This version number captures both DNAmod database and software versions which are kept in sync for release.
# "chebi_version" A constant holding the current version number of ChEBI

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
SITE_JS_DIR="$SITE_HTML_DIR/js"
JSON="$SITE_JS_DIR/lunr.json"
DATA_DIR="$SOURCE_DIR/data"

if [[ ! ( -d "$ROOT_DIR" && -d "$SITE_DIR" && -d "$SITE_TEMPLATE_DIR" && -d "$SITE_HTML_DIR" && -d "$SITE_IMAGE_DIR" && -d "$SITE_CSS_DIR" && -d "$SITE_JS_DIR" && -d "$DATA_DIR" ) ]]; then
    not_found 'directories'
fi

DATABASE="$DATA_DIR/DNAmod.sqlite"
DATABASE_COPY="$DATA_DIR/DNAmod_copy.sqlite"
WHITELIST="$DATA_DIR/whitelist.txt"
BLACKLIST="$DATA_DIR/blacklist.txt"

ANNOT_EXP_ALPH="$DATA_DIR/nomenclature.txt"
ANNOT_SEQ="$DATA_DIR/ref_annots_sequencing.txt"
ANNOT_NATURE="$DATA_DIR/ref_annots_nature.txt"
MANUAL_ADDITIONS="$DATA_DIR/manual_additions.txt"

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

SEQ_ANNOT_TABLE="sequencing_citations"
NATURE_ANNOT_TABLE="nucleobase_nature_info"
REFERENCES_TABLE="citations"
EXP_ALPH_TABLE="expanded_alphabet"

DNAMOD_VERSION="v1.6"
CHEBI_VERSION="01-01-2020" #From: https://www.ebi.ac.uk/ols/ontologies/chebi

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
    site_css_dir)
        echo -n "$SITE_CSS_DIR"
        ;;
    site_js_dir)
        echo -n "$SITE_JS_DIR"
        ;;
    json)
        echo -n "$JSON"
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
    database_copy)
        echo -n "$DATABASE_COPY"
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
    seq_annot_table)
        echo -n "$SEQ_ANNOT_TABLE"
        ;;
    nature_annot_table)
        echo -n "$NATURE_ANNOT_TABLE"
        ;;
    references_table)
        echo -n "$REFERENCES_TABLE"
        ;;
    exp_alph_table)
        echo -n "$EXP_ALPH_TABLE"
        ;;
    manual_additions)
        echo -n "$MANUAL_ADDITIONS"
        ;;
    dnamod_version)
        echo -n "$DNAMOD_VERSION"
        ;;
    chebi_version)
        echo -n "$CHEBI_VERSION"
        ;;
    *)
        # get all possible cases by searching this script
        >&2 echo "Usage: ./constants.sh <$(sed -n '/case/,/esac/p' $(readlink -f $0) | fgrep ')' | fgrep -v '*)' | fgrep -v 'sed' | cut -d ')' -f 1 | tr -d ' ' | paste -s -d '|')>"
        exit $ERR_EXIT
esac

