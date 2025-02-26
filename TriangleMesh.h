/*
Copyright (c) 2019, Michael Kazhdan
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

#ifndef TRIANGLE_MESH_INCLUDE
#define TRIANGLE_MESH_INCLUDE

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <Misha/Graph.h>
#include <Misha/SparseMatrix.h>
#include "MultiThreading.h"


namespace MishaK
{
	template< typename Real , unsigned int Dim , typename ... OtherPoints > Real Volume( Point< Real , Dim > point , OtherPoints ... otherPoints );
	template< typename Real , unsigned int Dim > SquareMatrix< Real , 3 > TriangleMassMatrix     ( const Point< Real , Dim > vertices[] );
	template< typename Real , unsigned int Dim > SquareMatrix< Real , 3 > TriangleStiffnessMatrix( const Point< Real , Dim > vertices[] );
	template< typename Real > SquareMatrix< Real , 3 > TriangleTutteMatrix( void );

	template< typename Real > void MakeLumped( SparseMatrix< Real , int > &m );
	template< typename Real > void ReorderMatrixEntries( SparseMatrix< Real , int > &M );

	template< typename VertexKey >
	struct TriangleMesh
	{

#if defined( DEBUG_GRAPH ) && defined( DEBUG_INDEX )
		typedef NumberWrapper< size_t , TriangleMesh , 0 > HalfEdgeIndex;
		typedef NumberWrapper< size_t , TriangleMesh , 1 > EdgeIndex;
#else // !DEBUG_GRAPH || !DEBUG_INDEX
		typedef size_t HalfEdgeIndex;
		typedef size_t EdgeIndex;
#endif // DEBUG_GRAPH && DEBUG_INDEX
		typedef size_t TriangleIndex;

		/** A class representing a triangle */
		struct Triangle
		{
			VertexKey &operator[]( unsigned int idx ){ return _vertices[idx%3]; }
			const VertexKey &operator[]( unsigned int idx ) const { return _vertices[idx%3]; }
		public:
			VertexKey _vertices[3];
		};

		// A representation of a half-edge in the triangle mesh (assumes there is at most one directed edge between two vertices)
		struct HalfEdge
		{
			VertexKey v1 , v2;
			HalfEdge( void ){}
			HalfEdge( VertexKey _v1 , VertexKey _v2 ) : v1(_v1) , v2(_v2){}
			HalfEdge opposite( void ) const { return HalfEdge( v2 , v1 ); }
			VertexKey opposite( VertexKey v ) const;

			bool operator==( HalfEdge he ) const { return v1==he.v1 && v2==he.v2; }
			struct Hash{ size_t operator()( HalfEdge key ) const { return std::hash< VertexKey >{}( key.v1 ) ^ std::hash< VertexKey >{}( key.v2 ); } };

			// Information about the corner opposite a half-edge
			struct OppositeCornerInfo
			{
				OppositeCornerInfo( void ){}
				OppositeCornerInfo( TriangleIndex t , char c ) : triangle(t) , corner(c){}
				TriangleIndex triangle;
				char corner;
			};
		};

		// The triangles in the mesh
		std::vector< Triangle > triangles;

		void reset( void );
		bool isWaterTight( void ) const;
		bool isBoundaryHalfEdge( HalfEdgeIndex eIndex ) const { return oppositeHalfEdgeIndex(eIndex)==-1; }
		bool isBoundaryVertex( VertexKey vKey ) const;

		size_t vertices( void ) const;
		size_t halfEdges( void ) const;
		size_t boundaryHalfEdges( void ) const;
		size_t edges( void ) const;
		int eulerCharacteristic( void ) const;
		size_t genus( void ) const;

		typename HalfEdge::OppositeCornerInfo oppositeCornerInfo( HalfEdgeIndex eIndex ) const { return _halfEdgeData[eIndex].oppositeCornerInfo; }
		typename HalfEdge::OppositeCornerInfo oppositeCornerInfo( HalfEdge he ) const { return _halfEdgeData[ halfEdgeIndex(he) ].oppositeCornerInfo; }

		Graph<     VertexKey , EdgeIndex >     graph( void ) const;
		Graph< TriangleIndex , EdgeIndex > dualGraph( void ) const;

		template< typename Real > WeightedGraph<     VertexKey , EdgeIndex , Real >     graph( const std::vector< Point3D< Real > > &vertices ) const;
		template< typename Real > WeightedGraph< TriangleIndex , EdgeIndex , Real > dualGraph( const std::vector< Point3D< Real > > &vertices ) const;

		HalfEdgeIndex halfEdgeIndex( HalfEdge he ) const;

		HalfEdge halfEdge( HalfEdgeIndex eIndex ) const;
#if defined( DEBUG_GRAPH ) && defined( DEBUG_INDEX )
		HalfEdge halfEdge(     EdgeIndex eIndex ) const { return halfEdge( HalfEdgeIndex( (size_t)eIndex ) ); }
#endif // DEBUG_GRAPH && DEBUG_INDEX

		EdgeIndex edgeIndex( HalfEdgeIndex eIndex ) const;
		EdgeIndex edgeIndex( HalfEdge he ) const { return edgeIndex( halfEdgeIndex( he ) ); }

		HalfEdgeIndex     nextHalfEdgeIndex( HalfEdgeIndex eIndex ) const { return _halfEdgeData[eIndex].nextHalfEdgeIndex; }
		HalfEdgeIndex previousHalfEdgeIndex( HalfEdgeIndex eIndex ) const { return _halfEdgeData[eIndex].previousHalfEdgeIndex; }
		HalfEdgeIndex oppositeHalfEdgeIndex( HalfEdgeIndex eIndex ) const { return _halfEdgeData[eIndex].oppositeHalfEdgeIndex; }

		HalfEdgeIndex     nextHalfEdgeIndex( HalfEdge he ) const { return _halfEdgeData[ halfEdgeIndex(he) ].nextHalfEdgeIndex; }
		HalfEdgeIndex previousHalfEdgeIndex( HalfEdge he ) const { return _halfEdgeData[ halfEdgeIndex(he) ].previousHalfEdgeIndex; }
		HalfEdgeIndex oppositeHalfEdgeIndex( HalfEdge he ) const { return _halfEdgeData[ halfEdgeIndex(he) ].oppositeHalfEdgeIndex; }

		std::vector< HalfEdgeIndex > halfEdgeOneRing( VertexKey v ) const;
		std::vector< VertexKey > vertexOneRing( VertexKey v ) const;
		std::vector< HalfEdgeIndex > face( TriangleIndex ) const;

		std::vector< HalfEdgeIndex > orient( const std::vector< EdgeIndex > &eSegment ) const;

		std::vector< HalfEdgeIndex > boundary( HalfEdgeIndex eIndex ) const;
		std::vector< std::vector< HalfEdgeIndex > > boundaries( void ) const;

		TriangleMesh split( const std::vector< HalfEdgeIndex > &path , const std::vector< VertexKey > &newVertexKeys ) const;

#if 0
		template< typename Real , typename TriangleMatrixFunctor >
		SparseMatrix< Real , int > matrix( TriangleMatrixFunctor F ) const;
#else
		template< typename Real , unsigned int Dim >
		SparseMatrix< Real , int > matrix( const std::vector< Point< Real , Dim > > &vertices , SquareMatrix< Real , 3 >(*F)( const Point< Real , Dim > [] ) ) const;
#endif
		template< typename Real , unsigned int Dim >
		SparseMatrix< Real , int > massMatrix( const std::vector< Point< Real , Dim > > &vertices , bool lump=false ) const;
		template< typename Real , unsigned int Dim >
		SparseMatrix< Real , int > stiffnessMatrix( const std::vector< Point< Real , Dim > > &vertices ) const;
		template< typename Real >
		SparseMatrix< Real , int > tutteLaplacian( void ) const;

	protected:
		struct _HalfEdgeData
		{
			HalfEdgeIndex oppositeHalfEdgeIndex , previousHalfEdgeIndex , nextHalfEdgeIndex;
			HalfEdge halfEdge;
			typename HalfEdge::OppositeCornerInfo oppositeCornerInfo;
			_HalfEdgeData( void ){}
			_HalfEdgeData( HalfEdge k , typename HalfEdge::OppositeCornerInfo o ) : halfEdge(k) , oppositeCornerInfo(o) , oppositeHalfEdgeIndex(-1) , previousHalfEdgeIndex(-1) , nextHalfEdgeIndex(-1){}
		};

#if defined( DEBUG_GRAPH ) && defined( DEBUG_INDEX )
		VectorWrapper< _HalfEdgeData , HalfEdgeIndex > _halfEdgeData;
#else // !DEBUG_GRAPH || !DEBUG_INDEX
		std::vector< _HalfEdgeData > _halfEdgeData;
#endif // DEBUG_GRAPH && DEBUG_INDEX
		std::unordered_map< HalfEdge , HalfEdgeIndex , typename HalfEdge::Hash > _halfEdgeToHalfEdgeIndex;
		std::unordered_map< VertexKey , HalfEdgeIndex > _vertexToHalfEdgeIndex;
		std::vector< HalfEdgeIndex > _triangleToHalfEdgeIndex;

		HalfEdgeIndex     _nextOutgoingHalfEdgeIndex( HalfEdgeIndex e ) const;
		HalfEdgeIndex _previousOutgoingHalfEdgeIndex( HalfEdgeIndex e ) const;
		HalfEdgeIndex _nextBoundaryHalfEdgeIndex( HalfEdgeIndex e ) const;
		HalfEdgeIndex _previousBoundaryHalfEdgeIndex( HalfEdgeIndex e ) const;
	};
#include "TriangleMesh.inl"
}
#endif // TRIANGLE_MESH_INCLUDE