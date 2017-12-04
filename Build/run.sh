#! /bin/sh

function printAndDo {
    echo "nbTasks $1
     outRunOpt $2
     alpha $3
     tol $4
     maxit $5
     Ftype $6
     Farg $7
     solverType $8
     meshFile $9
     outFFileName ${10}
     outUFileName ${11}
     outUeFileName ${12}
     outEFileName ${13}" | column -t | sed 's/ /./g'
    echo ARGV : $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13}
    make && mpirun -np $1 $2 ./solver $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13}
}

nbTasks=2

alpha=1
tol=1e-9
maxit=2
Ftype=1
Farg=1
solverType=0

NAME=test_$nbTasks

Wdir=benchmark
outWdir=$Wdir/$NAME

meshFile=$Wdir/carre_64.msh
outFFileName=$outWdir/solF
outUFileName=$outWdir/solNum
outUeFileName=$outWdir/solExa
outEFileName=$outWdir/solErr

outRunOpt="--

	" #"--output-filename run_output/$NAME"
mkdir 2>/dev/null outputs errors $outWdir
rm 2>/dev/null outputs/$NAME errors/$NAME $outWdir/*
cp $meshFile $outWdir

printAndDo $nbTasks $outRunOpt $alpha $tol $maxit $Ftype $Farg $solverType $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName
