#! /bin/sh

for nbTasks in $(seq $nbTasksMin $nbTasksStep $nbTasksMax)
do
    make VERBOSE=0 && mpirun -np $nbTasks ./solver $alpha $tol $maxit $Ftype $Farg $solverType $saveMshs $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName 
done
