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
#include <Misha/SubSimplexMesh.h>
#include <Misha/RightTriangleQuadrature.h>

namespace MishaK
{
	template< unsigned int K >
	struct LinearElements
	{
		static SquareMatrix< double , K > J( const SquareMatrix< double , K > & g );
		static Polynomial::Polynomial< K , 1 , double > Interpolant( const Point< double , K > p[K+1] , const double v[K+1] );
		static Point< double , K > Corner( unsigned int k );

		enum struct ScalarElementType
		{
			HAT ,
			CROUZEIX_RAVIART
		};
		static inline const std::vector ScalarElementTypeNames = { "Hat" , "Crouzeix-Raviart" };
		template< ScalarElementType Type > struct Mesh;

		struct Hat             { static Polynomial::Polynomial< K , 1 , double > Element( unsigned int k ); };
		struct CrouzeixRaviart { static Polynomial::Polynomial< K , 1 , double > Element( unsigned int k ); };

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

		template< enum ScalarElementType Type >
		struct Mesh : public std::conditional_t< Type==ScalarElementType::HAT , SubSimplexMesh< K , 0 > , std::conditional_t< Type==ScalarElementType::CROUZEIX_RAVIART , SubSimplexMesh< K , K-1 > , std::nullptr_t > >
		{
			static constexpr unsigned int SubK( void )
			{
				if      constexpr( Type==ScalarElementType::HAT ) return 0;
				else if constexpr( Type==ScalarElementType::CROUZEIX_RAVIART ) return K-1;
				else static_assert( false , "[ERROR] Unrecognized Type" );
			}
			using Elements = std::conditional_t< Type==ScalarElementType::HAT , HatElements , std::conditional_t< Type==ScalarElementType::CROUZEIX_RAVIART , CrouzeixRaviartElements , std::nullptr_t > >;
			using ElementIndex = std::conditional_t< Type==ScalarElementType::HAT , size_t , std::conditional_t< Type==ScalarElementType::CROUZEIX_RAVIART , typename SubSimplexMesh< K , SubK() >::SignedIndex , std::nullptr_t > >;
			using ElementSimplex = SimplexIndex< SubK() >;

			// Constructor:
			// [NOTE] The object accesses the reference to the simplices throughout its scope
			Mesh( const std::vector< SimplexIndex< K > > & simplices , size_t vNum );

			// Indexing functionality
			size_t elementNum( void ) const;

			ElementSimplex elementSimplex( size_t i ) const;

			ElementIndex operator()( size_t s , unsigned int k ) const;


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

			template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
			Eigen::SparseMatrix< double > elementToVertex( MetricFunctor && metric ) const;

			template< typename SampleFunctor /* = std::function< std::pair< size_t , Point< double , K > ( size_t ) > */ >
			Eigen::SparseMatrix< double > evaluation( size_t sampleNum , SampleFunctor && sampleFunctor ) const;

			// Sanity checking
			template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
			void sanityCheck( MetricFunctor && metric ) const;

		protected:
			using _SubSimplexMesh = SubSimplexMesh< K , SubK() >;
			using _SubSimplexIndex = _SubSimplexMesh::_SubSimplexIndex;
			using _SubSimplexMesh::_simplices;
			using _SubSimplexMesh::_subSimplexIndices;
			size_t _vNum;

			size_t _index( size_t s , unsigned int k ) const;
		};
	};
#include "LinearElements.inl"
}
#endif // LINEAR_ELEMENTS_INCLUDED