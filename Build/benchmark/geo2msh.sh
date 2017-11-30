if [[ $# -lt 4 ]]
then
    echo ./geo2msh NAME CLMIN CLMAX NBPARTS
    exit
fi

gmsh -2 -clmin $2 -clmax $3 $1 -part $4 -o ${1%.*}_${4}.msh -v 10
