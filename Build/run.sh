#! /bin/sh

nbTasks=2
Wdir=benchmark
output_file=log

alpha=0.1
tol=1e-6
maxit=1e5
Ftype=0
Farg=1
meshFile=$Wdir/carre_64.msh
outFFileName=$Wdir/solF
outUFileName=$Wdir/solNum
outUeFileName=$Wir/solExa
outEFileName=$Wdir/solErr

mkdir 2>/dev/null run_output
rm 2>/dev/null -f run_output/$output_file

make && mpirun -np $nbTasks --output-filename run_output/$output_file ./solver $alpha $tol $maxit $Ftype $Farg $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName 
