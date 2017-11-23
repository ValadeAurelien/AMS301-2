if [ $# -lt 3 ]
then
    echo Missing args : ./geo2msh NAME CLMIN CLMAX
else
    gmsh -2 -clmin $2 -clmax $3 -part 4 $1 -o ${1%.*}.msh
fi
