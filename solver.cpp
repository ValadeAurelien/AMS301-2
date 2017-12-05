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
    // if (!myRank) {
    // 	Matrix NN = N;
    // 	for(int n=0; n<m.nbOfNodes; n++) 
    // 	    for(int nn=0; nn<m.nbOfNodes; nn++)
    // 		NN(n, nn) = NN(n, nn) / Mdiag(n);
    // 	cout << NN.eigenvalues() << endl;
    // }
    // Jacobi solver
    double res_tot = 0;
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
        /*if(myRank == 1)
            cout << "norme = " << r.norm() << endl;
        sleep(myRank*3);*/
	computeL2Norm(res_tot, r, m, NO_PRINT);

	if((it % 100) == 0){
	    // printf("%i %i %f\n", myRank, it, r.norm());
	    // if (myRank==0) printf("-1 %i %f\n", it, res_tot);
	    if(myRank == 0){
	    	printf("it %6.d res %.15f\n", it, res_tot); 
	    }
	}
	it++;
    } while (res_tot > tol && it < maxit);
  
    if(myRank == 0) {
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
	p(A.rows()),
	Au(A.rows()),
	Ap(A.rows());
    const Vector nul = Matrix::Zero(A.rows(), 1);

    double res = 0,
	res2 = 0,
	res2_old = 0,
	res_tot = 0,
	alpha, beta;
    int it = 0;
    Au = A*u;
    exchangeAddInterfMPI(Au, m);
    r = b - Au; p = r; 
    do {
	res2 = r.dot(r);
	
	Ap = A*p;
	exchangeAddInterfMPI(Ap, m);		

	alpha =  res2 / ( p.transpose()*Ap );
	u += alpha * p;
	
	// cout << "a " << myRank << " " << alpha << endl;
	computeL2Norm(res_tot, u, m, NO_PRINT);
	cout << "u " << myRank << " " << res_tot << endl;
	res2_old = res2;
	r -= alpha*Ap;
	res2 = r.dot(r);
	
	beta = res2 / res2_old;
	p = r + beta * p;

	computeL2Norm(res_tot, r, m, NO_PRINT);
	if((it % 1) == 0) {
	    // printf("loc %i %i %f\n", myRank, it, r.norm());
	    if (myRank==0) printf("glob -1 %i %f\n", it, res_tot);
	    // if(myRank == 0){
	    // 	printf("it %6.d res %.15f\n", it, res_tot); 
	    // }
	}
	it++;
    } while (res_tot > tol && it < maxit);
  
    if(myRank == 0){!
	    printf("\r#   -> final iteration: %i (prescribed max: %i)\n", it, maxit);
	printf("#   -> final residual: %e (prescribed tol: %e)\n", sqrt(res2_old), tol);
    }
}

    // Au = A*u;
    // exchangeAddInterfMPI(Au, m);
    // r = b - Au; p = r; 
    // computeL2Norm(res_tot, r, m, NO_PRINT);
    // do {
    // 	res2_tot = pow(res_tot, 2);
	
    // 	Ap = A*p;
    // 	exchangeAddInterfMPI(Ap, m);		

    // 	alpha =  res2_tot / ( p.transpose()*Ap );
    // 	u += alpha * p;
	
    // 	// cout << "a " << myRank << " " << alpha << endl;
	
    // 	res2_tot_old = res2_tot;
    // 	r -= alpha*Ap;
    // 	computeL2Norm(res_tot, r, m, NO_PRINT);
    // 	res2_tot = pow(res_tot, 2);
	
    // 	beta = res2_tot / res2_tot_old;
    // 	p = r + beta * p;
	
    // 	if((it % 1) == 0){
    // 	    // printf("loc %i %i %f\n", myRank, it, r.norm());
    // 	    if (myRank==0) printf("glob -1 %i %f\n", it, res_tot);
    // 	    // if(myRank == 0){
    // 	    // 	printf("it %6.d res %.15f\n", it, res_tot); 
    // 	    // }
    // 	}
    // 	it++;
    // } while (res_tot > tol && it < maxit);
