#! /bin/sh

jobs=$(qstat | awk 'NR>2 {print $1}' | sed 's/\n/ /g')
Njobs=$(echo $jobs | wc -w)
if [[ Njobs -eq 0 ]]
then
    echo No job to delete
    exit
fi
echo -e Delete all jobs : $jobs    yes/no ?   
read confirmation

if [[ $confirmation == 'yes' ]]
then
    qdel $jobs
fi
    
