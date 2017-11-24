#! /bin/sh

jobs=$(qstat | awk 'NR>2 {print $1}' | sed 's/\n/ /g')
echo -e Deleting all jobs : $jobs  (yes/no) ?   
read confirmation

if [[ $confirmation == 'yes' ]]
then
    qdel $jobs
fi
