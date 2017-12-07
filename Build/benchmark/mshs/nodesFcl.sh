#!/bin/sh

logclmin=-2.5
logclmax=0
nbsteps=50
prec=7
infile=nodesFcl_init
outfile=nodesFcl_res.txt

echo "#cl nodes" > $outfile

for i in $(seq 1 $nbsteps)
do
    echo $i
    cl=$(calc "10^($logclmin+($logclmax $logclmin*(-1))*$i/$nbsteps)" | sed 's/~//g')
    cl=${cl:0:$prec}
    ./geo2msh.sh $infile.geo $cl $cl 1 1 1
    nodes=$(grep -A 1 "\$Nodes" ${infile}_1.msh | tail -n 1)
    echo $cl $nodes >> $outfile
done

column -t $outfile > ${outfile}_tmp
mv ${outfile}_tmp $outfile
