#! /bin/sh

for i in $(seq 1 $logtolnb)
do
    tol=$(python -c "print( 10**($logtolmin+($logtolmax $logtolmin*(-1))*$i/$logtolnb) )")
    tol=${tol:0:5}
    echo $i $tol
    make VERBOSE=0 && mpirun -np $nbTasks ./solver $alpha $tol $maxit $Ftype $Farg $solverType $saveMshs $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName 
done
