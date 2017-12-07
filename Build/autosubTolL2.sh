#! /bin/sh

nbTasks=16

alpha=1
logtolmin=-10
logtolmax=-3
logtolnb=10
maxit=1e6
Ftype=1
Farg=1
solverType=1
saveMshs=0

NAME=tolL2_gc_$nbTasks

Wdir=benchmark
outWdir=$Wdir/$NAME

meshFile=$Wdir/mshs/bigcarre_64.msh
outFFileName=$outWdir/solF
outUFileName=$outWdir/solNum
outUeFileName=$outWdir/solExa
outEFileName=$outWdir/solErr


ARGS_LIST="-v nbTasks=$nbTasks -v alpha=$alpha -v logtolmin=$logtolmin -v logtolmax=$logtolmax -v logtolnb=$logtolnb -v maxit=$maxit -v Ftype=$Ftype -v Farg=$Farg -v solverType=$solverType -v saveMshs=$saveMshs -v meshFile=$meshFile -v outFFileName=$outFFileName -v outUFileName=$outUFileName -v outUeFileName=$outUeFileName -v outEFileName=$outEFileName "

mkdir 2>/dev/null outputs errors $outWdir
rm 2>/dev/null outputs/$NAME errors/$NAME $outWdir/*
cp $meshFile $outWdir

qsub -N $NAME -o outputs/$NAME -e errors/$NAME -cwd -pe orte $nbTasks $ARGS_LIST ./.auxsubmitTolL2.sh
