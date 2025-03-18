/*
Copyright (c) 2022, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef EIGENVALUE_SOLVER_INCLUDE
#define EIGENVALUE_SOLVER_INCLUDE

#include <ARPACK++/arrgsym.h>
#include <Misha/LinearSolvers.h>
#include <Misha/SparseMatrix.h>
#pragma comment( lib , "arpack_x64.lib" )


//Solve Ax = lambda*Bx. A,B must be symmetric, A positive-semidefinite and B positive-definite
class SparseEigenProblem
{
	typedef Eigen::SparseMatrix< double > E_SparseMatrix; 
	typedef Eigen::Triplet< double > E_Triplet;
	typedef Eigen::VectorXd E_Vector;
#ifdef EIGEN_USE_MKL_ALL
	typedef Eigen::PardisoLDLT< E_SparseMatrix > E_Cholesky;
#else // !EIGEN_USE_MKL_ALL
	typedef Eigen::SimplicialLDLT< E_SparseMatrix > E_Cholesky;
#endif // EIGEN_USE_MKL_ALL
	typedef Eigen::MatrixXd E_Matrix;
	bool _computeEigenvectors( ARrcSymGenEig< double > & prob , const int numEigenvectors , std::vector< std::pair< double , std::vector< double > > >& spectrum , bool verbose );
public:
	template< class Real >
	SparseEigenProblem( const SparseMatrix< Real , int >& a , const SparseMatrix< Real , int >& b )
	{
		dimension = b.rows;
		A = E_SparseMatrix( dimension , dimension );
		B = E_SparseMatrix( dimension , dimension );
		std::vector< E_Triplet > triplets_A , triplets_B;
		triplets_A.reserve( a.Entries() ) , triplets_B.reserve( b.Entries() );
		for( int i=0 ; i<a.rows ; i++ ) for( int j=0 ; j<a.rowSizes[i] ; j++ ) triplets_A.push_back( E_Triplet( i , a[i][j].N , a[i][j].Value ) );
		for( int i=0 ; i<b.rows ; i++ ) for( int j=0 ; j<b.rowSizes[i] ; j++ ) triplets_B.push_back( E_Triplet( i , b[i][j].N , b[i][j].Value ) );
		A.setFromTriplets( triplets_A.begin() , triplets_A.end() );
		B.setFromTriplets( triplets_B.begin() , triplets_B.end() );
	}
	int dimension;
	E_SparseMatrix A;
	E_SparseMatrix B;
	/*
	mode:
	"LM" = Largest Magnitude (default)
	"SM" = Smallest Magnitude
	"LA" = Largest  Algebraic Value
	"SA" = Smallest Algebraic Value
	"BE" = Both End
	*/
	bool computePartialSpectrum_RegularMode(const int numEigenvectors, std::vector< std::pair< double , std::vector< double > > >& spectrum , char * mode = "LM" ,                               bool verbose=false );
	bool computePartialSpectrum_ShiftedMode(const int numEigenvectors, std::vector< std::pair< double , std::vector< double > > >& spectrum , char * mode = "LM" , double spectral_shift = 0.0 , bool verbose=false );
};

bool SparseEigenProblem::_computeEigenvectors( ARrcSymGenEig<double> & prob , const int numEigenvectors , std::vector< std::pair< double , std::vector< double > > >& spectrum , bool verbose )
{
	// Finding eigenvalues and eigenvectors.
	prob.FindEigenvectors();
	// Printing solution.
	int nconv;
	nconv = prob.ConvergedEigenvalues();
	if( verbose )
	{
		std::cout << "Real symmetric eigenvalue problem: A*x - B*x*lambda" << std::endl;
		std::cout << "Dimension of the system            : " << dimension << std::endl;
		std::cout << "Number of 'requested' eigenvalues  : " << numEigenvectors << std::endl;
		std::cout << "Number of 'converged' eigenvalues  : " << nconv << std::endl;
		std::cout << "Number of Arnoldi vectors generated: " << prob.GetNcv() << std::endl;
		std::cout << "Number of iterations taken         : " << prob.GetIter() << std::endl;
		std::cout << std::endl;
	}

	if( !prob.EigenvaluesFound() ) return false;
	if( nconv < numEigenvectors ) return false;

	spectrum.resize( numEigenvectors );
	for( int j=0 ; j<numEigenvectors ; j++ )
	{
		spectrum[j].first = prob.Eigenvalue(j);
		spectrum[j].second.resize( dimension );
		for( int i=0 ; i<dimension ; i++ ) spectrum[j].second[i] = prob.Eigenvector( j , i );
	}
	return true;
}

bool SparseEigenProblem::computePartialSpectrum_RegularMode( const int numEigenvectors , std::vector< std::pair< double , std::vector< double > > >& spectrum , char * mode , bool verbose )
{
	if( verbose ) printf( "Computing Cholesky(B)\n" );
	E_Cholesky B_Cholesky(B);
	if( verbose ) printf( "Solving Eigenvalue Problem\n" );
	ARrcSymGenEig< double > prob( dimension , numEigenvectors , mode );
	E_Vector E_getVector(dimension);
	E_Vector E_putVector(dimension);
	double* A_getVector;
	double* A_putVector;
	// Finding an Arnoldi basis.
	int iter = 0;
	while( !prob.ArnoldiBasisFound() )
	{
		iter++;
		if( verbose ) printf( "Arnoldi Iter %07d \r" , iter );
		prob.TakeStep();
		int probIdo = prob.GetIdo();
		if ((probIdo == 1) || (probIdo == -1) || (probIdo == 2)){
			A_getVector = prob.GetVector();
			A_putVector = prob.PutVector();
			for (int i = 0; i < dimension; i++) E_getVector[i] = A_getVector[i];
			if ((probIdo == 1) || (probIdo == -1)) {
				//w <-inv(B)*A*v, v<-A*v
				E_getVector = A*E_getVector;
				E_putVector = B_Cholesky.solve(E_getVector);
				for (int i = 0; i < dimension; i++) A_getVector[i] = E_getVector[i];
			}
			else if (probIdo == 2) {
				//w <- B*v.
				E_putVector = B*E_getVector;
			}
			for (int i = 0; i < dimension; i++) A_putVector[i] = E_putVector[i];
		}
	}
	return _computeEigenvectors( prob , numEigenvectors , spectrum , verbose );
}
bool SparseEigenProblem::computePartialSpectrum_ShiftedMode( const int numEigenvectors , std::vector< std::pair< double , std::vector< double > > >& spectrum , char * mode , double spectral_shift , bool verbose )
{
	if( verbose ) printf( "Computing Cholesky( A - B * %g )\n" , spectral_shift );
	E_Cholesky A_minus_shift_B_Cholesky(A - spectral_shift*B);
	if( verbose ) printf( "Solving Eigenvalue Problem\n" );
	ARrcSymGenEig<double> prob('S', dimension, numEigenvectors, spectral_shift, mode);
	
	E_Vector E_getVector( dimension );
	E_Vector E_putVector( dimension );

	double* A_getVector;
	double* A_putVector;
	// Finding an Arnoldi basis.
	int iter = 0;
	while (!prob.ArnoldiBasisFound()){
		iter++;
		if( verbose ) printf( "Arnoldi Iter %07d \r" , iter );
		prob.TakeStep();
		int probIdo = prob.GetIdo();
		if ((probIdo == -1) || (probIdo == 1) || (probIdo == 2)){
			A_putVector = prob.PutVector();
			if (probIdo == -1){
					//w <-inv(A-shift*B)*B*v
					A_getVector = prob.GetVector();
					for (int i = 0; i < dimension; i++) E_getVector[i] = A_getVector[i];
					E_getVector = B*E_getVector;
					E_putVector = A_minus_shift_B_Cholesky.solve(E_getVector);
			}
			else if (probIdo == 1){
					//w <-inv(A-shift*B)*u
					A_getVector = prob.GetProd();
					for (int i = 0; i < dimension; i++) E_getVector[i] = A_getVector[i];
					E_putVector = A_minus_shift_B_Cholesky.solve(E_getVector);
			}
			else{//(probIdo == 2)
				//w <- B*v.
				A_getVector = prob.GetVector();
				for (int i = 0; i < dimension; i++) E_getVector[i] = A_getVector[i];
				E_putVector = B*E_getVector;
			}
			for (int i = 0; i < dimension; i++) A_putVector[i] = E_putVector[i];
		}
	}
	return _computeEigenvectors( prob , numEigenvectors , spectrum , verbose );
}
#endif// EIGENVALUE_SOLVER_INCLUDE