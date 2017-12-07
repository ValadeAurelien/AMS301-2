#! /bin/sh

alpha=1 
tol=1e-9 
maxit=1e5 
Ftype=1 
Farg=1 
solverType=1 
saveMshs=0 
 
NAME=WeakSca_${solverType}
Wdir=benchmark
outWdir=$Wdir/$NAME

outFFileName=$outWdir/solF
outUFileName=$outWdir/solNum
outUeFileName=$outWdir/solExa
outEFileName=$outWdir/solErr


ARGS_LIST="-v Wdir=$Wdir -v alpha=$alpha -v tol=$tol -v maxit=$maxit -v Ftype=$Ftype -v Farg=$Farg -v solverType=$solverType -v saveMshs=$saveMshs -v outFFileName=$outFFileName -v outUFileName=$outUFileName -v outUeFileName=$outUeFileName -v outEFileName=$outEFileName "

mkdir 2>/dev/null outputs errors $outWdir
rm 2>/dev/null outputs/$NAME errors/$NAME $outWdir/*
cp $meshFile $outWdir

qsub -N $NAME -o outputs/$NAME -e errors/$NAME -cwd -pe orte $nbTasksMax $ARGS_LIST ./.auxsubmitStrongSca.sh
