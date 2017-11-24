#! /bin/sh

make && mpirun -np $nbTasks ./solver $alpha $tol $maxit $Ftype $Farg $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName
