#include "headers.hpp"
#include <stdlib.h>

int myRank;
int nbTasks;

int main(int argc, char* argv[]) {

    if (argc<11) {
	cout << "#Missing args, got :";
	for (int i=0; i<argc; ++i ) cout << " " << argv[i];
	cout << endl;
	return EXIT_FAILURE;
    }
    double alpha     = atof(argv[1]);
    double tol       = atof(argv[2]);
    int    maxit     = floor(atof(argv[3]));
    int    Ftype     = atof(argv[4]);
    double Farg      = atof(argv[5]);
    char*  meshFile  = argv[6];
    char*  outFFile  = argv[7];
    char*  outUFile  = argv[8];
    char*  outUeFile = argv[9];
    char*  outEFile  = argv[10];

// 1. Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nbTasks);

    // 2. Read the mesh, and build lists of nodes for MPI exchanges, local numbering
    Mesh m;
    readMsh(m, meshFile);
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
	switch (Ftype) {
	case CSTE : 
	    uExa(i) = alpha*Farg;
	    f(i) = Farg;
	    break;
	case COSCOS :
	    f(i) = cos(M_PI * x) * cos(Farg * M_PI * y);
	    uExa(i) = (alpha + pow(Farg*M_PI, 2)) * f(i);
	    break;
	default :
	    cout << "Type of F not understood..." << endl;
	    return EXIT_FAILURE;
	}
    }

    saveToMsh(uNum, m, "solF", outFFile);

    Problem p;
    buildLinearSystem(p, m, alpha, f);
    // buildDirichletBC(p,m,uExa); // (Only for the extension of the project)

    // 4. Solve problem
    jacobi(p.A, p.b, uNum, m, tol, maxit);

    // 5. Compute error and export fields

    double L2_err;
    computeL2Err(L2_err, uNum, uExa, m);
    Vector uErr = (uNum-uExa).cwiseAbs();

    saveToMsh(uNum, m, "solF", outFFile);
    saveToMsh(uNum, m, "solNum", outUFile);
    saveToMsh(uExa, m, "solRef", outUeFile);
    saveToMsh(uErr, m, "solErr", outEFile);

    // 6. Finilize MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
