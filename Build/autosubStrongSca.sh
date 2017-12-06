#! /bin/sh

nbTasksMin=0
nbTasksStep=4
nbTasksMax=64

alpha=1
tol=1e-6
maxit=1e6
Ftype=1
Farg=1
solverType=0
saveMshs=0

NAME=test_$nbTasks

Wdir=benchmark
outWdir=$Wdir/$NAME

meshFile=$Wdir/bigcarre_64.msh
outFFileName=$outWdir/solF
outUFileName=$outWdir/solNum
outUeFileName=$outWdir/solExa
outEFileName=$outWdir/solErr


ARGS_LIST="-v nbTasksMin=$nbTasksMin -v nbTasksStep=$nbTasksStep -v nbTasksMax=$nbTasksMax -v alpha=$alpha -v tol=$tol -v maxit=$maxit -v Ftype=$Ftype -v Farg=$Farg -v solverType=$solverType -v saveMshs=$saveMshs -v meshFile=$meshFile -v outFFileName=$outFFileName -v outUFileName=$outUFileName -v outUeFileName=$outUeFileName -v outEFileName=$outEFileName "

mkdir 2>/dev/null outputs errors $outWdir
rm 2>/dev/null outputs/$NAME errors/$NAME $outWdir/*
cp $meshFile $outWdir

qsub -N $NAME -o outputs/$NAME -e errors/$NAME -cwd -pe orte $nbTasks $ARGS_LIST ./.auxsubmitStrongSca.sh
