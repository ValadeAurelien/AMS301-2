if [[ $# -lt 3 ]]
then
    echo ./geo2msh NAME CLMIN CLMAX
    exit
fi

gmsh -2 -clmin $2 -clmax $3 $1 -o ${1%.*}.msh -v 10
