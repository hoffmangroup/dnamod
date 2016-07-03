#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

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

