#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

# Creates or re-creates DNAmod, constructing the database and then the website.
#
# -------------------------------------------------------------------------------
# Copyright (C) 2016  Ankur Jai Sood and Coby Viner
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

SOURCE_DIR="$(dirname $(readlink -f $0))" # From: https://gist.github.com/tvlooy/cbfbdb111a4ebad8b93e
CONSTANTS_SCRIPT="$SOURCE_DIR/constants.sh"

# Updates DNAmod database and generates static site
echo "Begining to update DNAmod and create static site..."

# re-build the database from scratch
$($CONSTANTS_SCRIPT 'script_pop_db' -r)
echo "Database updated."

$($CONSTANTS_SCRIPT 'script_create_site')

echo -e "\nSite creation complete.\n\nRun $($CONSTANTS_SCRIPT 'script_sync_site') to sync changes."
