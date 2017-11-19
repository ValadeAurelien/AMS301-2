#include "headers.hpp"

int myRank;
int nbTasks;

int main(int argc, char* argv[]) {

    // 1. Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nbTasks);

    // 2. Read the mesh, and build lists of nodes for MPI exchanges, local numbering
    Mesh m;
    readMsh(m, "benchmark/carre.msh");
    buildListsNodesMPI(m);
    buildLocalNumbering(m);

    // 3. Build problem (fields and system)
    Vector uNum(m.nbOfNodes);
    Vector uExa(m.nbOfNodes);
    Vector f(m.nbOfNodes);
    for (int i = 0; i < m.nbOfNodes; ++i) {
        double x = m.coords(i, 0);
        double y = m.coords(i, 1);
        uNum(i) = 0.;
        uExa(i) = cos(M_PI * x) * cos(2 * M_PI * y);
        f(i) = x*y;
    }

    Problem p;
    double alpha = 1;
    buildLinearSystem(p, m, alpha, f);
    // buildDirichletBC(p,m,uExa); // (Only for the extension of the project)

    // 4. Solve problem
    double tol = 1e-9; // (Currently useless)
    int maxit = 1000;
    jacobi(p.A, p.b, uNum, m, tol, maxit);

    // 5. Compute error and export fields

    double L2_err;
    computeL2Err(L2_err, uNum, uExa, m);

    saveToMsh(uNum, m, "solNum", "benchmark/solNum.msh");
    saveToMsh(uExa, m, "solRef", "benchmark/solExa.msh");
    //saveToMsh(uErr, m, "solErr", "benchmark/solErr.msh");

    // 6. Finilize MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
