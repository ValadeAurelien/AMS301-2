#! /bin/sh

for i in $(seq 1 $logalphanb)
do
    echo $i
    alpha=$(calc "10^($logclmin+($logalphamax $logalphain*(-1))*$i/$logalphanb)" | sed 's/~//g')
    alpha=${alpha:0:5}
    make VERBOSE=0 && mpirun -np $nbTasks ./solver $alpha $tol $maxit $Ftype $Farg $solverType $saveMshs $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName 
done
