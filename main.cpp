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
    double	alpha	   = atof(argv[1]);
    double	tol	   = atof(argv[2]);
    int		maxit	   = floor(atof(argv[3]));
    int		Ftype	   = atof(argv[4]);
    double	Farg	   = atof(argv[5]);
    int		solverType = atoi(argv[6]);
    char*	meshFile   = argv[7];
    char*	outFFile   = argv[8];
    char*	outUFile   = argv[9];
    char*	outUeFile  = argv[10];
    char*	outEFile   = argv[11];

    // 1. Initialize MPI
    clock_t c0 = clock();
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nbTasks);

    // 2. Read the mesh, and build lists of nodes for MPI exchanges, local numbering
    Mesh m;
    readMsh(m, meshFile);
    buildListsNodesMPI(m);
    int totNbOfNodes = m.nbOfNodes;
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
	    uExa(i) = Farg/alpha;
	    f(i) = Farg;
	    break;
	case COSCOS :
	    f(i) = cos(M_PI * x) * cos(Farg * M_PI * y);
	    uExa(i) = f(i)/(alpha + (1 + pow(Farg, 2))*M_PI*M_PI);
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
    switch(solverType) {
    case JACOBI:
	jacobi(p.A, p.b, uNum, m, tol, maxit);
	break;
    case CONJUGATEGRADIENT:
	conjugateGradient(p.A, p.b, uNum, m, tol, maxit);
	break;
    default:
	printf("Wrong solver choice...");
	return -1;
    }
    // 5. Compute error and export fields

    double L2_err;
    computeL2Err(L2_err, uNum, uExa, m, PRINT);
    Vector uErr = (uNum-uExa).cwiseAbs();
    if (!myRank)
	cout << "#time " << ((double) clock() - c0)/CLOCKS_PER_SEC << endl;
    saveToMsh(f, m, "solF", outFFile);
    saveToMsh(uNum, m, "solNum", outUFile);
    saveToMsh(uExa, m, "solRef", outUeFile);
    saveToMsh(uErr, m, "solErr", outEFile);

    // 6. Finilize MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
