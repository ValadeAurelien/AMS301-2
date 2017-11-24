#! /bin/sh

nbTasks=16

alpha=1
tol=1e-3
maxit=1e3
Ftype=1
Farg=1

Wdir=benchmark
meshFile=$Wdir/bigcarre.msh
outFFile=$Wdir/solF.msh
outUFile=$Wdir/solNum.msh
outUeFile=$Wdir/solExa.msh
outEFile=$Wdir/solErr.msh

NAME=test_$nbTasks

ARGS_LIST="-v nbTasks=$nbTasks -v alpha=$alpha -v tol=$tol -v maxit=$maxit -v Ftype=$Ftype -v Farg=$Farg -v Wdir=$Wdir -v meshFile=$meshFile -v outFFile=$outFFile -v outUFile=$outUFile -v outUeFile=$outUeFile -v outEFile=$outEFile"

mkdir 2>/dev/null outputs errors 
rm 2>/dev/null outputs/$NAME errors/$NAME
qsub -N $NAME -o outputs/$NAME -e errors/$NAME -cwd -pe orte $nbTasks $ARGS_LIST ./.auxsubmit.sh
