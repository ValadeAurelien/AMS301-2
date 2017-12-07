#! /bin/sh

for nbTasks in 1 2 3 4 8 12 16
do
    meshFile=$Wdir/mshs/weaksca_${nbTasks}.msh
    echo "#ARGS : $alpha $tol $maxit $Ftype $Farg $solverType $saveMshs $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName"
    make VERBOSE=0 && mpirun -np $nbTasks ./solver $alpha $tol $maxit $Ftype $Farg $solverType $saveMshs $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName 
done
