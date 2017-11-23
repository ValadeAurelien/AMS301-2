#! /bin/sh

nbTasks=2

alpha=1
tol=1e-3
maxit=1e3
Ftype=1
Farg=1

Wdir=benchmark
meshFile=$Wdir/bigcarre.msh
outFFile=$Wdir/solF.msh
outUFile=$Wdir/solNum.msh
outUeFile=$Wir/solExa.msh
outEFile=$Wdir/solErr.msh

make && mpirun -np $nbTasks ./solver $alpha $tol $maxit $Ftype $Farg $meshFile $outFFile $outUFile $outUeFile $outEFile
