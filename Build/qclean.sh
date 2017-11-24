#! /bin/sh

jobs=$(qstat | awk 'NR>2 {print $1}' | sed 's/\n/ /g')
echo Deleting all jobs : $jobs
read confirmation

if [[ $confirmation == 'yes' ]]
   qdel $jobs
fi
