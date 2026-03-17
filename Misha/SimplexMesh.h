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

#ifndef SIMPLEX_MESH_INCLUDED
#define SIMPLEX_MESH_INCLUDED
#include <iostream>
#include <vector>
#include <functional>
#include <Eigen/Sparse>
#include "Geometry.h"
#include "SimplexBasis.h"
#include "MultiThreading.h"
#include "Atomic.h"

namespace MishaK
{
	template< unsigned int ... > struct SimplexMesh;

	template< unsigned int Dim >
	struct SimplexMesh< Dim >
	{
		struct Sample
		{
			unsigned int sIdx;
			Point< double , Dim+1 > bcCoordinates;
		};
	};

	template< unsigned int Dim , unsigned int Degree >
	struct SimplexMesh< Dim , Degree >
	{
#if 1 // NEW_CODE
		using NodeMultiIndex = MultiIndex< Degree , unsigned int , false >;
		using FaceMultiIndex = MultiIndex< Dim , unsigned int , false >;
		using NodeMultiIndexMap = typename MultiIndex< Degree , unsigned int , false >::template map< unsigned int >;
		using FaceMultiIndexMap = typename MultiIndex< Dim , unsigned int , false >::template map< unsigned int >;
#else // !NEW_CODE
		typedef MultiIndex< Degree , unsigned int , false > NodeMultiIndex;
		typedef MultiIndex< Dim , unsigned int , false > FaceMultiIndex;
#endif // NEW_CODE
		static const unsigned int NodesPerSimplex = SimplexElements< Dim , Degree >::NodeNum;

		SimplexMesh( void ){}

		template< unsigned int EmbeddingDimension , typename Index , typename VertexFunctor /* = std::function< Point< double , EmbeddingDimension > ( size_t ) > */ >
		static SimplexMesh Init( const std::vector< SimplexIndex< Dim , Index > > &simplices , VertexFunctor && vFunction );

		template< typename Index , typename MetricFunctor /* = std::function< SquareMatrix< double , Dim > ( size_t ) > */ >
		static SimplexMesh Init( const std::vector< SimplexIndex< Dim , Index > > &simplices , MetricFunctor && gFunction );

		size_t simplices( void ) const { return _simplices.size(); }
		SimplexIndex< Dim , unsigned int > simplex( unsigned int idx ) const{ return _simplices[idx]; }
		SquareMatrix< double , Dim > metric( unsigned int idx ) const { return _g[idx]; }
		size_t nodes( void ) const { return _nodeMap.size(); }
		Eigen::SparseMatrix< double > mass( void ) const;
		Eigen::SparseMatrix< double > stiffness( void ) const;
		Eigen::SparseMatrix< double > bistiffness( void ) const;
		Eigen::SparseMatrix< double > system( std::function< SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum > ( const SquareMatrix< double , Dim > & ) > m2s ) const;
		Eigen::SparseMatrix< double > crossFaceGradientEnergy( void ) const;
		template< typename UseFaceFunctor >
		Eigen::SparseMatrix< double > crossFaceGradientEnergy( const UseFaceFunctor &useFaceFunctor ) const;
		template< unsigned int FaceDim >
		Eigen::SparseMatrix< double > _crossFaceGradientEnergy( void ) const;
		Eigen::SparseMatrix< double > evaluationMatrix( const std::vector< typename SimplexMesh< Dim >::Sample > &samples ) const;
		double evaluate( const Eigen::VectorXd &coefficients , typename SimplexMesh< Dim >::Sample sample ) const;
		Polynomial::Polynomial< Dim , Degree , double > evaluate( const Eigen::VectorXd &coefficients , unsigned int simplexIndex ) const;

		SimplexMesh( SimplexMesh &&sm ){ std::swap( _simplices , sm._simplices ) , std::swap( _g , sm._g ) , std::swap( _nodeMap , sm._nodeMap ); }
		SimplexMesh & operator = ( SimplexMesh &&sm ){ std::swap( _simplices , sm._simplices ) , std::swap( _g , sm._g ) , std::swap( _nodeMap , sm._nodeMap ) ; return *this; }

		NodeMultiIndex nodeMultiIndex( unsigned int s , unsigned int n ) const;
		unsigned int nodeIndex( const NodeMultiIndex & multiIndex ) const;
		unsigned int nodeIndex( unsigned int s , unsigned int n ) const { return _localToGlobalNodeIndex.size() ? _localToGlobalNodeIndex[ s*NodesPerSimplex+n ] : nodeIndex( nodeMultiIndex( s , n ) ); }
#if 1 // NEW_CODE
		typename NodeMultiIndexMap::const_iterator cbegin( void ) const { return _nodeMap.cbegin(); }
		typename NodeMultiIndexMap::const_iterator cend  ( void ) const { return _nodeMap.cend  (); }
		const NodeMultiIndexMap &nodeMap( void ) const { return _nodeMap; }
#else // !NEW_CODE
		typename NodeMultiIndex::map::const_iterator cbegin( void ) const { return _nodeMap.cbegin(); }
		typename NodeMultiIndex::map::const_iterator cend  ( void ) const { return _nodeMap.cend  (); }
		const typename NodeMultiIndex::map &nodeMap( void ) const { return _nodeMap; }
#endif // NEW_CODE

		double volume( void ) const;
		void makeUnitVolume( void );
		void hashLocalToGlobalNodeIndex( void );

		template< typename MetricFunctor /* = std::function< SquareMatrix< double , Dim > ( size_t ) > */ >
		void updateMetric( MetricFunctor && gFunction );

		template< unsigned int EmbeddingDimension , typename VertexFunctor /* = std::function< Point< double , EmbeddingDimension > ( size_t ) > */ >
		void updateMetricFromPositions( VertexFunctor && vFunction );

#if 1 // NEW_CODE
		template< typename Index >
		static Eigen::SparseMatrix< double > Prolongation( const std::vector< SimplexIndex< Dim , Index > > & simplices );
#endif // NEW_CODE

	protected:
		template< unsigned int EmbeddingDimension , typename Index , typename VertexFunctor /* = std::function< Point< double , EmbeddingDimension > ( size_t ) > */ >
		void _initFromPositions( const std::vector< SimplexIndex< Dim , Index > > &simplices , VertexFunctor && vFunction );

		template< typename Index , typename MetricFunctor /* = std::function< SquareMatrix< double , Dim > ( size_t ) > */ >
		void _initFromMetric( const std::vector< SimplexIndex< Dim , Index > > &simplices , MetricFunctor && gFunction );

#if 1 // NEW_CODE
		template< typename Index >
		static NodeMultiIndex _NodeMultiIndex( unsigned int s , unsigned int n , const std::vector< SimplexIndex< Dim , Index > > & simplices );

		template< typename Index >
		static NodeMultiIndexMap _NodeMap( const std::vector< SimplexIndex< Dim , Index > > &simplices );
#endif // NEW_CODE

		std::vector< SimplexIndex< Dim , unsigned int > > _simplices;
#if 1 // NEW_CODE
		NodeMultiIndexMap _nodeMap;
#else // !NEW_CODE
		typename NodeMultiIndex::map _nodeMap;
#endif // NEW_CODE
		std::vector< SquareMatrix< double , Dim > > _g;
		std::vector< unsigned int > _localToGlobalNodeIndex;

#if 1 // NEW_CODE
		template< unsigned int ... > friend struct SimplexMesh;
#endif // NEW_CODE
	};

#include "SimplexMesh.inl"
}
#endif // SIMPLEX_MESH_INCLUDED