#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

# Generates a template for manually-curated nucleobase annotations.
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

SOURCE_DIR="$(dirname $(readlink -f $0))" # From: https://gist.github.com/tvlooy/cbfbdb111a4ebad8b93e
CONSTANTS_SCRIPT="$SOURCE_DIR/constants.sh"

q_result=$(sqlite3 $("$CONSTANTS_SCRIPT" 'database') <<EOF

SELECT names.nameid, names.chebiname, base.commonname
FROM modbase
JOIN names ON names.nameid=modbase.nameid
JOIN base ON modbase.baseid=base.baseid
WHERE verifiedstatus = 1
ORDER BY modbase.baseid, names.chebiname COLLATE NOCASE;

EOF
)

echo "$q_result" | awk -F '|' '{if ($3 != prev_unmod) {printf("#\n# %s\n#\n", $3)}; printf("# %s\n%s\t\n", $2, $1); prev_unmod = $3}' | sed 's/CHEBI://'

