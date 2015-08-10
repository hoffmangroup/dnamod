#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

# Updates DNAmod database and generates static site
echo "Begining to update DNAmod and create static site"
cd /mnt/work1/users/home2/asood/DNA_Base_Database
python populate_database_sql.py
echo "Database updated"
cd /mnt/work1/users/home2/asood/DNA_Base_Database/DNA_mod_site
python Create_mod_staticsite.py
echo "Static site created"

