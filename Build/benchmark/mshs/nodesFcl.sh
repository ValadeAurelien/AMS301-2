#!/bin/sh

logclmin=-2
logclmax=0
nbsteps=10
infile=nodeFcl_init
outfile=nodesFcl_res
echo "#cl nodes" > $outfile

for i in $(seq 1 $nbsteps)
do
    echo $i
    echo "10**($logclmin+($logclmax $logclmin)*$i/$nbsteps)"
    cl=$(calc "10^($logclmin+($logclmax-$logclmin)*$i/$nbsteps)")
    echo $cl
    ./geo2msh $infile.geo $cl $cl 1 1 1
    nodes=$(grep -A 1 "\$Nodes" $infile_1.msh | tail -n 1)
    echo $cl $nodes >> $outfile
done

column -t $outfile > $outfile_tmp
mv $outfile_tmp $outfile
