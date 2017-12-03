#! /bin/sh

nbTasks=2

alpha=1
tol=1e-9
maxit=1e5
Ftype=1
Farg=1
solverType=1

NAME=test_$nbTasks

Wdir=benchmark
outWdir=$Wdir/$NAME

meshFile=$Wdir/carre_64.msh
outFFileName=$outWdir/solF
outUFileName=$outWdir/solNum
outUeFileName=$outWdir/solExa
outEFileName=$outWdir/solErr

outRunOpt="" #"--output-filename run_output/$NAME"
mkdir 2>/dev/null outputs errors $outWdir
rm 2>/dev/null outputs/$NAME errors/$NAME $outWdir/*
cp $meshFile $outWdir

make && mpirun -np $nbTasks $outRunOpt ./solver $alpha $tol $maxit $Ftype $Farg $solverType $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName
