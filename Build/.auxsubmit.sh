#! /bin/sh

make && mpirun -np $nbTasks ./solver $alpha $tol $maxit $Ftype $Farg $meshFile $outFFile $outUFile $outUeFile $outEFile
