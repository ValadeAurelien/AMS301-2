#! /bin/sh

nbTasks=2
Wdir=benchmark
output_file=log

alpha=1
tol=1e-3
maxit=1e3
Ftype=1
Farg=1
meshFile=$Wdir/carre.msh
outFFile=$Wdir/solF.msh
outUFile=$Wdir/solNum.msh
outUeFile=$Wir/solExa.msh
outEFile=$Wdir/solErr.msh

mkdir 2>/dev/null run_output
rm 2>/dev/null -f run_output/$output_file

make && mpirun -np $nbTasks --output-filename run_output/$output_file ./solver $alpha $tol $maxit $Ftype $Farg $meshFile $outFFile $outUFile $outUeFile $outEFile 
