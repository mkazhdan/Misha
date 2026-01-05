/*
Copyright (c) 2025, Michael Kazhdan and Matthew Bolitho
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

#ifndef MG_INCLUDED
#define MG_INCLUDED

#include <Eigen/Sparse>
#include "MultiThreading.h"

namespace MishaK
{
#if 1 // NEW_CODE
	struct JacobiRelaxer
	{
		JacobiRelaxer( const Eigen::SparseMatrix< double > &M );
		void operator()( const Eigen::VectorXd &b , Eigen::VectorXd &x ) const;
		void operator()( const Eigen::MatrixXd &B , Eigen::MatrixXd &X ) const;
	protected:

		const Eigen::SparseMatrix< double > &_M;
		Eigen::VectorXd _D;
	};
#endif // NEW_CODE

	struct GSRelaxer
	{
		GSRelaxer( const Eigen::SparseMatrix< double > &M );
		void operator()( const Eigen::VectorXd &b , Eigen::VectorXd &x ) const;
		void operator()( const Eigen::VectorXd &b , Eigen::VectorXd &x , const std::vector< std::vector< unsigned int > > &mcIndices ) const;
#if 1 // NEW_CODE
		void operator()( const Eigen::MatrixXd &B , Eigen::MatrixXd &X ) const;
		void operator()( const Eigen::MatrixXd &B , Eigen::MatrixXd &X , const std::vector< std::vector< unsigned int > > &mcIndices ) const;
#endif // NEW_CODE
	protected:

		const Eigen::SparseMatrix< double > &_M;
		Eigen::VectorXd _D;
	};

#if 1 // NEW_CODE
	template< typename Solver , bool Jacobi=false >
#else // !NEW_CODE
	template< typename Solver >
#endif // NEW_CODE
	struct MGSolver
	{
		template< typename ProlongationFunctor /* = std::function< Eigen::SparseMatrix< double > ( unsigned int ) > */ >
		MGSolver( const Eigen::SparseMatrix< double > &M , ProlongationFunctor Ps , size_t pNum );
		~MGSolver( void );

#if 1 // NEW_CODE
		void vCycle( const Eigen::VectorXd & b , Eigen::VectorXd & x , unsigned int iters );

		void vCycle( const Eigen::MatrixXd & B , Eigen::MatrixXd & X , unsigned int iters );

		template< typename MCIndicesFunctor /* = std::function< const std::vector< std::vector< unsigned int > > & ( unsigned int ) > */ >
		void vCycle( const Eigen::VectorXd &b , Eigen::VectorXd &x , unsigned int iters , MCIndicesFunctor MCI );

		template< typename MCIndicesFunctor /* = std::function< const std::vector< std::vector< unsigned int > > & ( unsigned int ) > */ >
		void vCycle( const Eigen::MatrixXd & B , Eigen::MatrixXd & X , unsigned int iters , MCIndicesFunctor MCI );
#else // !NEW_CODE
		void vCycle( const Eigen::VectorXd &b , Eigen::VectorXd &x , unsigned int iters );

		template< typename MCIndicesFunctor /* = std::function< const std::vector< std::vector< unsigned int > > & ( unsigned int ) > */ >
		void vCycle( const Eigen::VectorXd &b , Eigen::VectorXd &x , unsigned int iters , MCIndicesFunctor MCI );
#endif // NEW_CODE

	protected:
#if 1 // NEW_CODE
		template< typename T >
		void _vCycle( const T & b , T & x , unsigned int iters );

		template< typename T , typename MCIndicesFunctor /* = std::function< const std::vector< std::vector< unsigned int > > & ( unsigned int ) > */ >
		void _vCycle( const T & b , T & x , unsigned int iters , MCIndicesFunctor MCI );
		Eigen::SparseMatrix< double > _M;
		std::vector< std::conditional_t< Jacobi , JacobiRelaxer , GSRelaxer > > _relaxers;
#else // !NEW_CODE
		const Eigen::SparseMatrix< double > &_M;
		std::vector< GSRelaxer > _gsRelaxers;
		std::vector< Eigen::VectorXd > _Xs , _Bs;
#endif // NEW_CODE
		std::vector< Eigen::SparseMatrix< double > > _Ps , _Rs , _Ms;
		Solver *_solver;
	};
#include "MG.inl"
}
#endif // MG_INCLUDED