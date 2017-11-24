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
      r(A.rows());
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
  double res2 = 0,
      tol2 = tol*tol;
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
    res2 = pow(r.norm(), 2);
    

    if((it % 10) == 0){
	if(myRank == 0){
	    cout << it << " " << res2 << endl;
	}
    }

    if (myRank)
	MPI_Reduce(&res2, &res2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    else
	MPI_Reduce(MPI_IN_PLACE, &res2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    it++;
  } while (res2 > tol2 && it < maxit);
  
  if(myRank == 0){
    printf("\r   -> final iteration: %i (prescribed max: %i)\n", it, maxit);
    printf("   -> final residual: %e (prescribed tol: %e)\n", sqrt(res2), tol);
  }
}

