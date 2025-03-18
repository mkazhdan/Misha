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

#ifndef GRAPH_INCLUDED
#define GRAPH_INCLUDED

#define DEBUG_GRAPH
#define NEW_GRAPH

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <limits>
#include <type_traits>
#include <functional>
#include <Misha/UnionFind.h>
#ifdef DEBUG_GRAPH
#include <Misha/NumberWrapper.h>
#endif // DEBUG_GRAPH

namespace MishaK
{
	// An edge
	template< typename VertexKey , typename EdgeKey >
	struct Edge
	{
		VertexKey v1 , v2;
		EdgeKey e;
		VertexKey opposite( VertexKey vKey ) const
		{
			if( v1==vKey ) return v2;
			if( v2==vKey ) return v1;
			MK_THROW( "vertex not in edge" );
		}
	};
	// A weighted edge
	template< typename VertexKey , typename EdgeKey , typename Real >
	struct WeightedEdge : public Edge< VertexKey , EdgeKey >{ Real weight; };

	template< typename VertexKey , typename EdgeKey >
	struct Graph
	{
#ifdef DEBUG_GRAPH
		typedef NumberWrapper< size_t , Graph , 0 > VertexIndex;
		typedef NumberWrapper< size_t , Graph , 1 > EdgeIndex;
#else // !DEBUG_GRAPH
		typedef size_t EdgeIndex;
		typedef size_t VertexIndex;
#endif // DEBUG_GRAPH

		/** Constructs a graph from a list of edges */
		Graph( const std::vector<         Edge< VertexKey , EdgeKey > > &edges );

		/** Constructs a graph from a list of weighted edges */
		template< typename Real >
		Graph( const std::vector< WeightedEdge< VertexKey , EdgeKey , Real > > &edges );

		/** The number of vertices in the graph */
		size_t vertices( void ) const { return _indexToVertex.size(); }

		/** The number of edges in the graph */
		size_t edges( void ) const { return _indexToEdge.size(); }

		/** The number of edges incident on vertex vKey */
		size_t edges( VertexKey vKey ) const { return _neighborEdges[ _vertexIndex(vKey) ].size(); }

		/** Returns the vertex key associated to the idx-th vertex in the graph */
		VertexKey vertexKey( VertexIndex idx ) const;

		/** Returns the edge key associated to the idx-th edge in the graph */
		EdgeKey edgeKey( EdgeIndex idx ) const;

		/** The idx-th edge incident on vertex vKey */
		EdgeKey operator() ( VertexKey vKey , size_t idx ) const { return edge( _neighborEdges[ _vertexIndex( vKey ) ][idx] ); }

		/** The opposite vertex of the idx-th edge incident on the vertex vKey */
		VertexKey neighbor( VertexKey vKey , size_t idx ) const { VertexIndex v = _vertexIndex(vKey) ; return vertexKey( _edges[ _neighborEdges[v][idx] ].opposite(v) ); }

		/** Removes an edge */
		void remove( EdgeKey eKey );

		/** Removes all edges whose incident vertices have valence one */
		template< typename RemovableNodeFunctor >
		void removeDanglingNodes( RemovableNodeFunctor Removable );
		void removeDanglingNodes( void ){ return removeDanglingNodes( []( VertexKey ){ return true; } ); }

		/** Groups together contiguous edges whose incident vertices have valence two */
		std::vector< std::vector< EdgeKey > > segments( void ) const;

		/** Returns the edge associated to a key */
		Edge< VertexKey , EdgeKey > edge( EdgeKey key ) const;

		/** Returns a spanning tree */
		Graph< VertexKey , EdgeKey > spanningTree( bool complement=false ) const;

		/** Returns true if the graph is connected */
		bool isConnected( void ) const;

		/** Returns true if the graph is connected using only the valid edges */
		bool isConnected( std::function< bool ( EdgeKey ) > ValidEdgeFunction ) const;

	protected:
		// Gives the set of undirected edges emanating from a vertex (accessed through indices to vertices, and giving indices to edges)
#ifdef DEBUG_GRAPH
		VectorWrapper< std::vector< EdgeIndex > , VertexIndex > _neighborEdges;
#else // !DEBUG_GRAPH
		std::vector< std::vector< EdgeIndex > > _neighborEdges;
#endif // DEBUG_GRAPH

		// Map from vertices to indices
		std::unordered_map< VertexKey , VertexIndex > _vertexToIndex;

		// Map from edges to indices
		std::unordered_map< EdgeKey , EdgeIndex > _edgeToIndex;

		// Map from indices to vertices
#ifdef DEBUG_GRAPH
		VectorWrapper< VertexKey , VertexIndex > _indexToVertex;
#else // !DEBUG_GRAPH
		std::vector< VertexKey > _indexToVertex;
#endif // DEBUG_GRAPH

		// Map from indices to edges
#ifdef DEBUG_GRAPH
		VectorWrapper< EdgeKey , EdgeIndex > _indexToEdge;
#else // !DEBUG_GRAPH
		std::vector< EdgeKey > _indexToEdge;
#endif // DEBUG_GRAPH

		// The list of edges
#ifdef DEBUG_GRAPH
		VectorWrapper< Edge< VertexIndex , EdgeIndex > , EdgeIndex > _edges;
#else // !DEBUG_GRAPH
		std::vector< Edge< VertexIndex , EdgeIndex > > _edges;
#endif // DEBUG_GRAPH

		Graph( void ){}

		template< typename _Edge > void _set( const std::vector< _Edge > &edges );

		VertexIndex _vertexIndex( VertexKey vKey ) const;

		EdgeIndex _edgeIndex( EdgeKey eKey ) const;

		void _remove( EdgeIndex e );

		void _reset( void );

		std::vector< std::vector< EdgeIndex > > _segments( void ) const;
	};

	template< typename VertexKey , typename EdgeKey , typename Real >
	struct DijkstraVertex
	{
		VertexKey vertex;
		EdgeKey previousEdge;
		Real distanceToSource;
		DijkstraVertex( VertexKey v=-1 , EdgeKey p=-1 , Real d=std::numeric_limits< Real >::infinity() ) : vertex(v) , previousEdge(p) , distanceToSource(d) {}
	};


	template< typename VertexKey , typename EdgeKey , typename Real >
	struct WeightedGraph : public Graph< VertexKey , EdgeKey >
	{
		typedef typename WeightedGraph::EdgeIndex EdgeIndex;
		typedef typename WeightedGraph::VertexIndex VertexIndex;
#ifdef NEW_GRAPH
		struct DijkstraData
		{
			friend struct WeightedGraph;

		protected:
			VertexIndex _source;
#ifdef DEBUG_GRAPH
			VectorWrapper< DijkstraVertex< VertexIndex , EdgeIndex , Real > , VertexIndex > _vertices;
#else // !DEBUG_GRAPH
			std::vector< DijkstraVertex< VertexIndex , EdgeIndex , Real > > _vertices;
#endif // DEBUG_GRAPH
			Real _distanceToSource( VertexIndex target ) const { return _vertices[ target ].distanceToSource; }

			DijkstraData( void ) : _source(-1){}
		};
#endif // NEW_GRAPH

		WeightedGraph( const std::vector< WeightedEdge< VertexKey , EdgeKey , Real > > &edges );
		Real weight( EdgeKey e ) const { return _weights[ _edgeIndex(e) ]; }

#ifdef NEW_GRAPH
		DijkstraData dijkstra( VertexKey source ,                    std::function< bool (EdgeKey) > ValidEdgeFunction=[]( EdgeKey ){ return true; } ) const;
		DijkstraData dijkstra( VertexKey source , VertexKey target , std::function< bool (EdgeKey) > ValidEdgeFunction=[]( EdgeKey ){ return true; } ) const;
#else // !NEW_GRAPH
		std::vector< DijkstraVertex< VertexKey , EdgeKey , Real > > dijkstra( VertexKey source ) const;
		std::vector< DijkstraVertex< VertexKey , EdgeKey , Real > > dijkstra( VertexKey source , VertexKey target ) const;
#endif // NEW_GRAPH
		std::vector< EdgeKey > shortestPath( VertexKey source , VertexKey target , std::function< bool (EdgeKey) > ValidEdgeFunction=[]( EdgeKey ){ return true; } ) const;
#ifdef NEW_GRAPH
		std::vector< EdgeKey > shortestPath( const DijkstraData &dData , VertexKey target ) const;
#endif // NEW_GRAPH

		/** Returns the edges in a shortest cycle in the graph */
		std::vector< EdgeKey > shortestCycle( std::function< bool (EdgeKey) > ValidEdgeFunction=[]( EdgeKey ){ return true; } ) const;

		/** Returns a minimal spanning tree */
		WeightedGraph< VertexKey , EdgeKey , Real > minimalSpanningTree( bool complement=false ) const;


	protected:
#ifdef DEBUG_GRAPH
		VectorWrapper< Real , EdgeIndex > _weights;
#else // !DEBUG_GRAPH
		std::vector< Real > _weights;
#endif // DEBUG_GRAPH
#ifdef NEW_GRAPH
		DijkstraData _dijkstra( VertexIndex source , VertexIndex target , std::function< bool (EdgeKey) > ValidEdgeFunction ) const;
		std::vector< EdgeIndex > _shortestPath( const DijkstraData &dData , VertexIndex target ) const;
#else // !NEW_GRAPH
#ifdef DEBUG_GRAPH
		VectorWrapper< DijkstraVertex< VertexIndex , EdgeIndex , Real > , VertexIndex > _dijkstra( VertexIndex source , VertexIndex target , std::function< bool (EdgeKey) > ValidEdgeFunction ) const;
#else // !DEBUG_GRAPH
		std::vector< DijkstraVertex< VertexIndex , EdgeIndex , Real > > _dijkstra( VertexIndex source , VertexIndex target , std::function< bool (EdgeKey) > ValidEdgeFunction ) const;
#endif // DEBUG_GRAPH
#endif // NEW_GRAPH
	};
#include "Graph.inl"
}
#endif // GRAPH_INCLUDED