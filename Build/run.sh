#! /bin/sh

nbTasks=1

alpha=1
tol=1.5e-6
maxit=1e4
Ftype=1
Farg=1

NAME=test_$nbTasks

Wdir=benchmark
outWdir=$Wdir/$NAME

meshFile=$Wdir/carre_64.msh
outFFileName=$outWdir/solF
outUFileName=$outWdir/solNum
outUeFileName=$outWdir/solExa
outEFileName=$outWdir/solErr

mkdir 2>/dev/null outputs errors $outWdir
rm 2>/dev/null outputs/$NAME errors/$NAME $outWdir/*
cp $meshFile $outWdir

make && mpirun -np $nbTasks --output-filename run_output/$NAME ./solver $alpha $tol $maxit $Ftype $Farg $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName 
