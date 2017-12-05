#! /bin/sh

function printAndDo {
    echo "nbTasks $1
     alpha $2
     tol $3
     maxit $4
     Ftype $5
     Farg $6
     solverType $7
     meshFile $8
     outFFileName $9
     outUFileName ${10}
     outUeFileName ${11}
     outEFileName ${12}" | column -t | sed 's/ /./g'
    echo ARGV : $2 $3 $4 $5 $6 $7 $8 ${9} ${10} ${11} ${12}
    make && (mpirun -np $1 ./solver $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12}) #> log
}

nbTasks=16

alpha=1
tol=1e-9
maxit=1
Ftype=0
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

mkdir 2>/dev/null outputs errors $outWdir
rm 2>/dev/null outputs/$NAME errors/$NAME $outWdir/*
cp $meshFile $outWdir

printAndDo $nbTasks $outRunOpt $alpha $tol $maxit $Ftype $Farg $solverType $meshFile $outFFileName $outUFileName $outUeFileName $outEFileName
