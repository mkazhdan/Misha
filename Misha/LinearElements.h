/*
Copyright (c) 2026, Michael Kazhdan
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

#ifndef LINEAR_ELEMENTS_INCLUDED
#define LINEAR_ELEMENTS_INCLUDED

#include <vector>

#include <Eigen/Sparse>

#include <Misha/Geometry.h>
#include <Misha/Exceptions.h>
#include <Misha/Miscellany.h>
#include <Misha/Polynomial.h>
#include <Misha/SimplexBasis.h>
#include <Misha/FEM.h>

namespace MishaK
{
	template< unsigned int K >
	struct LinearElements
	{
		static SquareMatrix< double , K > J( const SquareMatrix< double , K > & g );
		static Polynomial::Polynomial< K , 1 , double > Interpolant( const Point< double , K > p[K+1] , const double v[K+1] );
		static Point< double , K > Corner( unsigned int k );

		struct Hat { static Polynomial::Polynomial< K , 1 , double > Element( unsigned int k ); };
		struct CrouzeixRaviart { static Polynomial::Polynomial< K , 1 , double > Element( unsigned int k ) ; static SimplexIndex< K-1 > Boundary( unsigned int k ); };

	protected:
		template< typename EType >
		struct _Elements
		{
			_Elements( void );
			const Polynomial::Polynomial< K , 1 , double > operator[]( unsigned int idx ) const { return _elements[idx]; }
			SquareMatrix< double , K+1 > mass( void ) const;
			SquareMatrix< SquareMatrix< double , K > , K+1 > stiffness( void ) const;
			Point< Point< double , K > , K+1 , double > differential( void ) const;

			static SquareMatrix< double , K+1 > Mass( void );
			static SquareMatrix< SquareMatrix< double , K > , K+1 > Stiffness( void );
			static Point< Point< double , K > , K+1 , double > Differential( void );
		protected:
			Polynomial::Polynomial< K , 1 , double > _elements[K+1];
		};

		template< bool InverseMetricTensor >
		static SquareMatrix< double , K > _J( const SquareMatrix< double , K > & g );

	public:
		using             HatElements = _Elements< Hat >;
		using CrouzeixRaviartElements = _Elements< CrouzeixRaviart >;

		template< bool UseHat >
		struct Mesh
		{
			using ElementIndex = std::conditional_t< UseHat , size_t , std::pair< size_t , bool > >;
			using ElementSimplex = std::conditional_t< UseHat , SimplexIndex< 0 > , SimplexIndex< K-1 > >;

			// Constructor:
			// [NOTE] The object accesses the reference to the simplices throughout its scope
			Mesh( const std::vector< SimplexIndex< K > > & simplices , size_t vNum );

			// Indexing functionality
			size_t elementNum( void ) const;

			ElementSimplex elementSimplex( size_t i ) const;

			std::optional< ElementIndex > operator()( size_t s , unsigned int k ) const;


			// Scalar system matrices
			template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
			Eigen::SparseMatrix< double > mass( MetricFunctor && metric ) const;

			template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
			Eigen::SparseMatrix< double > stiffness( MetricFunctor && metric ) const;


			// Derivative system matrices
			Eigen::SparseMatrix< double > differential( void ) const;

			template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
			Eigen::SparseMatrix< double > differentialMass( MetricFunctor && metric ) const;

			template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
			Eigen::SparseMatrix< double > differentialJ( MetricFunctor && metric ) const;


			// Sanity checking
			template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
			void sanityCheck( MetricFunctor && metric ) const;

		protected:
			using _FaceIndex = MultiIndex< K , size_t , true >;

			const std::vector< SimplexIndex< K > > & _simplices;
			std::map< _FaceIndex , size_t > _simplexFaceMap;
			std::vector< _FaceIndex > _simplexFaces;
			size_t _vNum;

			size_t _index( size_t s , unsigned int k ) const;
		};
	};
#include "LinearElements.inl"
}
#endif // LINEAR_ELEMENTS_INCLUDED