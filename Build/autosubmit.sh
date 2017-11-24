#! /bin/sh

nbTasks=16

alpha=1
tol=1e-3
maxit=1e3
Ftype=1
Farg=1

NAME=test_$nbTasks

Wdir=benchmark
outWdir=$Wdir/$NAME

meshFile=$Wdir/bigcarre_64.msh
outFFileName=$outWdir/solF
outUFileName=$outWdir/solNum
outUeFileName=$outWdir/solExa
outEFileName=$outWdir/solErr


ARGS_LIST="-v nbTasks=$nbTasks -v alpha=$alpha -v tol=$tol -v maxit=$maxit -v Ftype=$Ftype -v Farg=$Farg -v Wdir=$Wdir -v meshFile=$meshFile -v outFFileName=$outFFileName -v outUFileName=$outUFileName -v outUeFileName=$outUeFileName -v outEFileName=$outEFileName"

mkdir 2>/dev/null outputs errors $outWdir
rm 2>/dev/null outputs/$NAME errors/$NAME $outWdir/*
qsub -N $NAME -o outputs/$NAME -e errors/$NAME -cwd -pe orte $nbTasks $ARGS_LIST ./.auxsubmit.sh
