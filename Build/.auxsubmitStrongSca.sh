#! /bin/sh

for nbTasks in $(seq $nbTasksMin $nbTasksStep $nbTasksMax)
do
    echo "#ARGS : $alpha $tol $maxit $Ftype $Farg $solverType $saveMshs $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName"
    make VERBOSE=0 && mpirun -np $nbTasks ./solver $alpha $tol $maxit $Ftype $Farg $solverType $saveMshs $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName 
done
