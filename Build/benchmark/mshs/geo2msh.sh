if [[ $# -lt 6 ]]
then
    echo ./geo2msh NAME CLMIN CLMAX NBPARTSMIN NBPARTSSTEP NBPARTSTEP
    exit
fi

for i in $(seq $4 $5 $6)
do
    echo "--- $i ---"
    echo "gmsh -2 -clmin $2 -clmax $3 $1 -part $i -o ${1%.*}_${i}.msh -v 10"
    gmsh -2 -clmin $2 -clmax $3 $1 -part $i -o ${1%.*}_${i}.msh -v 10
done
