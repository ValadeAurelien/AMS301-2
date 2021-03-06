#! /bin/sh

for i in $(seq 1 $logalphanb)
do
    alpha=$(python -c "print( 10**($logalphamin+($logalphamax $logalphamin*(-1))*$i/$logalphanb) )")
    alpha=${alpha:0:5}
    echo $i $alpha
    make VERBOSE=0 && mpirun -np $nbTasks ./solver $alpha $tol $maxit $Ftype $Farg $solverType $saveMshs $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName 
done
