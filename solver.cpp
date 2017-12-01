#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Solution of the system Au=b with Jacobi
//================================================================================

void jacobi(SpMatrix& A, Vector& b, Vector& u, Mesh& m, double tol, int maxit)
{
    if(myRank == 0)
	printf("#== jacobi\n");
  
    // Compute the solver matrices
    Vector Mdiag(A.rows()),
	r(A.rows()),
	nul = Matrix::Zero(A.rows(), 1);
    SpMatrix N(A.rows(), A.cols());
    for(int k = 0; k < A.outerSize(); ++k){
	for(SpMatrix::InnerIterator it(A,k); it; ++it){
	    if(it.row() == it.col())
		Mdiag(it.row()) = it.value();
	    else
		N.coeffRef(it.row(), it.col()) = -it.value();
	}
    }
    exchangeAddInterfMPI(Mdiag, m);
  
    // Jacobi solver
    double res_tot;
    int it = 0;
    do {
	// Compute N*u
	Vector Nu = N*u;
	exchangeAddInterfMPI(Nu, m);
	Nu += b;
	// Update field
	for(int n=0; n<m.nbOfNodes; n++)
	    u(n) = 1/Mdiag(n) * Nu(n);

	r = b - A*u;
	if (nbTasks>1)
	    computeL2Err(res_tot, r, nul, m, NO_PRINT);
	else
	    res_tot = sqrt(r.norm());
    
	if((it % 100) == 0){
	    if(myRank == 0){
		printf("it %6.d res %.15f\n", it, res_tot);
	    }
	}
	it++;
    } while (res_tot > tol && it < maxit);
  
    if(myRank == 0){
	printf("\r#   -> final iteration: %i (prescribed max: %i)\n", it, maxit);
	printf("#   -> final residual: %e (prescribed tol: %e)\n", res_tot, tol);
    }
}

void conjugateGradient(SpMatrix& A, Vector& b, Vector& u, Mesh& m, double tol, int maxit)
{
    if(myRank == 0)
	printf("#== Conjugate Gradient\n");
  
    // Compute the solver matrices
    Vector r(A.rows()),
	rold(A.rows()),
	p(A.rows()),
	Au(A.rows()),
	Ap(A.rows()),
	nul = Matrix::Zero(A.rows(), 1);
    Au = A*u;
    exchangeAddInterfMPI(Au, m);
    r = b - Au; p = r; rold = r;
    // Jacobi solver
    double res2 = 0,
	res2old = 0,
	res_tot,
	alpha, beta;
    int it = 0;
    do {
	res2 = r.dot(r);
	
	Ap = A*p;
	exchangeAddInterfMPI(Ap, m);

	alpha =  res2 / ( p.transpose()*Ap );
	u += alpha * p;
	
	rold = r; res2old = res2;
	r -= alpha*Ap;
	
	res2 = r.dot(r);
	beta = res2 / res2old;
	p = r + beta * p;

	MPI_Barrier(MPI_COMM_WORLD);
	if (nbTasks>1)
	    computeL2Err(res_tot, r, nul, m, NO_PRINT);
	else
	    res_tot = sqrt(res2);
	
	if((it % 10) == 0){
	    if(myRank == 0){
		printf("it %6.d res %.15f\n", it, res_tot);
	    }
	}
	it++;
    } while (res_tot > tol && it < maxit);
  
    if(myRank == 0){!
	    printf("\r#   -> final iteration: %i (prescribed max: %i)\n", it, maxit);
	printf("#   -> final residual: %e (prescribed tol: %e)\n", sqrt(res2old), tol);
    }
}

