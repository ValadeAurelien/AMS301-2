#! /bin/sh

make VERBOSE=0 && mpirun -np $nbTasks ./solver $alpha $tol $maxit $Ftype $Farg $solverType $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName 
