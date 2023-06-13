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

///////////
// Graph //
///////////
template< typename VertexKey , typename EdgeKey >
template< typename _Edge >
void Graph< VertexKey , EdgeKey >::_set( const std::vector< _Edge > &edges )
{
	static_assert( std::is_base_of< Edge< VertexKey , EdgeKey > , _Edge >::value , "[ERROR] _Edge must derive from EdgeKey" );

	_neighborEdges.clear();
	_vertexToIndex.clear();
	_edgeToIndex.clear();
	_indexToVertex.clear();
	_indexToEdge.clear();
	_edges.clear();

	for( size_t e=0 ; e<edges.size() ; e++ )
	{
		if( _vertexToIndex.find( edges[e].v1 )==_vertexToIndex.end() ) _vertexToIndex[ edges[e].v1 ] = _vertexToIndex.size();
		if( _vertexToIndex.find( edges[e].v2 )==_vertexToIndex.end() ) _vertexToIndex[ edges[e].v2 ] = _vertexToIndex.size();
		if( _edgeToIndex.find( edges[e].e )==_edgeToIndex.end() ) _edgeToIndex[ edges[e].e ] = e;
		else THROW( "multiple instances of the same edge" );
	}

	_indexToVertex.resize( _vertexToIndex.size() );
	_indexToEdge.resize( _edgeToIndex.size() );
	for( auto iter=_vertexToIndex.begin() ; iter!=_vertexToIndex.end() ; iter++ ) _indexToVertex[ iter->second ] = iter->first;
	for( auto iter=_edgeToIndex.begin() ; iter!=_edgeToIndex.end() ; iter++ ) _indexToEdge[ iter->second ] = iter->first;

	_edges.resize( edges.size() );
	for( size_t e=0 ; e<edges.size() ; e++ )
	{
		_edges[e].v1 = _vertexIndex( edges[e].v1 );
		_edges[e].v2 = _vertexIndex( edges[e].v2 );
		_edges[e].e = _edgeIndex( edges[e].e );
	}

	_neighborEdges.resize( _vertexToIndex.size() );
	for( size_t e=0 ; e<edges.size() ; e++ )
	{
		VertexIndex v1 = _vertexIndex( edges[e].v1 ) , v2 = _vertexIndex( edges[e].v2 );
		_neighborEdges[v1].push_back( e );
		_neighborEdges[v2].push_back( e );
	}
}

template< typename VertexKey , typename EdgeKey >
void Graph< VertexKey , EdgeKey >::_reset( void )
{
	std::unordered_map< EdgeIndex , Edge< VertexIndex , EdgeIndex > > edgeMap;
	std::vector< Edge< VertexKey , EdgeKey > > edges;
	for( VertexIndex v=0 ; v<_neighborEdges.size() ; v++ ) for( size_t i=0 ; i<_neighborEdges[v].size() ; i++ ) edgeMap[ _neighborEdges[v][i] ] = _edges[ _neighborEdges[v][i] ];

	edges.reserve( edgeMap.size() );
	for( auto iter=edgeMap.begin() ; iter!=edgeMap.end() ; iter++ )
	{
		Edge< VertexKey , EdgeKey > edge;
		edge.e = _indexToEdge[ iter->second.e ];
		edge.v1 = _indexToVertex[ iter->second.v1 ];
		edge.v2 = _indexToVertex[ iter->second.v2 ];
		edges.push_back( edge );
	}
	_set( edges );
}

template< typename VertexKey , typename EdgeKey >
Graph< VertexKey , EdgeKey >::Graph( const std::vector< Edge< VertexKey , EdgeKey > > &edges ){ _set( edges ); }

template< typename VertexKey , typename EdgeKey >
template< typename Real >
Graph< VertexKey , EdgeKey >::Graph( const std::vector< WeightedEdge< VertexKey , EdgeKey , Real > > &edges ){ _set( edges ); }

template< typename VertexKey , typename EdgeKey >
VertexKey Graph< VertexKey , EdgeKey >::vertexKey( VertexIndex idx ) const
{
	if( idx>=_indexToVertex.size() ) THROW( "index out of bounds: " , (size_t)idx , " >= " , _indexToVertex.size() );
	return _indexToVertex[idx];
}

template< typename VertexKey , typename EdgeKey >
EdgeKey Graph< VertexKey , EdgeKey >::edgeKey( EdgeIndex idx ) const
{
	if( idx>=_indexToEdge.size() ) THROW( "index out of bounds: " , (size_t)idx , " >= " , _indexToEdge.size() );
	return _indexToEdge[idx];
}

template< typename VertexKey , typename EdgeKey >
Edge< VertexKey , EdgeKey > Graph< VertexKey , EdgeKey >::edge( EdgeKey key ) const
{
	Edge< VertexIndex , EdgeIndex > _e = _edges[ _edgeIndex(key) ];
	Edge< VertexKey , EdgeKey > e;
	e.v1 = vertexKey( _e.v1 );
	e.v2 = vertexKey( _e.v2 );
	e.e = edgeKey( _e.e );
	return e;
}

template< typename VertexKey , typename EdgeKey >
void Graph< VertexKey , EdgeKey >::_remove( EdgeIndex e )
{
	// [WARNING] The edge is still referenced in _edgeToIndex, _indexToEdge, and _edges
	Edge< VertexIndex , EdgeIndex > edge = _edges[e];
	VertexIndex v1 = edge.v1 , v2 = edge.v2;

	for( int j=(int)_neighborEdges[v1].size()-1 ; j>=0 ; j-- ) if( _neighborEdges[v1][j]==e )
	{
		_neighborEdges[v1][j] = _neighborEdges[v1].back();
		_neighborEdges[v1].pop_back();
	}
	for( int j=(int)_neighborEdges[v2].size()-1 ; j>=0 ; j-- ) if( _neighborEdges[v2][j]==e )
	{
		_neighborEdges[v2][j] = _neighborEdges[v2].back();
		_neighborEdges[v2].pop_back();
	}
}

template< typename VertexKey , typename EdgeKey >
void Graph< VertexKey , EdgeKey >::remove( EdgeKey eKey )
{
	// [WARNING] The edge is still referenced in _edgeToIndex, _indexToEdge, and _edges
	_remove( _edgeIndex(eKey) );
}

template< typename VertexKey , typename EdgeKey >
template< typename RemovableNodeFunctor >
void Graph< VertexKey , EdgeKey >::removeDanglingNodes( RemovableNodeFunctor Removable )
{
	bool removeDangling = true;
	while( removeDangling )
	{
		removeDangling = false;
		for( VertexIndex v=0 ; v<_neighborEdges.size() ; v++ ) if( _neighborEdges[v].size()==1 && Removable( vertexKey( v ) ) )
		{
			EdgeIndex e = _neighborEdges[v][0];
			Edge< VertexIndex , EdgeIndex > edge = _edges[e];
			VertexIndex o = edge.opposite( v );
			for( int j=(int)_neighborEdges[o].size()-1 ; j>=0 ; j-- ) if( _edges[ _neighborEdges[o][j] ].opposite(o)==v )
			{
				_neighborEdges[o][j] = _neighborEdges[o].back();
				_neighborEdges[o].pop_back();
			}
			if( _neighborEdges[o].size()==1 ) removeDangling = true;
			_neighborEdges[v].resize(0);
		}
	}
	_reset();
}

template< typename VertexKey , typename EdgeKey >
std::vector< std::vector< typename Graph< VertexKey , EdgeKey >::EdgeIndex > > Graph< VertexKey , EdgeKey >::_segments( void ) const
{
#ifdef DEBUG_GRAPH
	VectorWrapper< bool , VertexIndex > processed( _neighborEdges.size() , false );
#else // !DEBUG_GRAPH
	std::vector< bool > processed( _neighborEdges.size() , false );
#endif // DEBUG_GRAPH
	std::vector< std::vector< EdgeIndex > > segments;
	auto GrowSegment = [&]( VertexIndex vIndex , EdgeIndex eIndex , std::vector< EdgeIndex > &segment )
	{
		segment.resize(0);
		segment.push_back( eIndex );
		VertexIndex nextVIndex = _edges[ eIndex ].opposite( vIndex );
		while( _neighborEdges[nextVIndex].size()==2 )
		{
			if( _neighborEdges[nextVIndex][0]==eIndex )
			{
				vIndex = nextVIndex;
				eIndex = _neighborEdges[nextVIndex][1];
			}
			else
			{
				vIndex = nextVIndex;
				eIndex = _neighborEdges[nextVIndex][0];
			}
			nextVIndex = _edges[ eIndex ].opposite( vIndex );
			segment.push_back(eIndex);
		}
	};

	for( VertexIndex v=0 ; v<_neighborEdges.size() ; v++ )
	{
		if( _neighborEdges[v].size()==2 && !processed[v] )
		{
			std::vector< EdgeIndex > front , back , segment;
			GrowSegment( v , _neighborEdges[v][0] , front );
			GrowSegment( v , _neighborEdges[v][1] , back );
			segment.reserve( back.size() + front.size() - 1 );
			for( int j=(int)back.size()-1 ; j>=0 ; j-- ) processed[ _edges[ back[j] ].v1 ] = processed[ _edges[ back[j] ].v2 ] = true , segment.push_back( back[j] );
			for( size_t j=0 ; j<front.size() ; j++ ) processed[ _edges[ front[j] ].v1 ] = processed[ _edges[ front[j] ].v2 ] = true , segment.push_back( front[j] );
			segments.push_back( segment );
		}
	}
	return segments;
}

template< typename VertexKey , typename EdgeKey >
std::vector< std::vector< EdgeKey > > Graph< VertexKey , EdgeKey >::segments( void ) const
{
	std::vector< std::vector< EdgeIndex > > _segments = this->_segments();
	std::vector< std::vector< EdgeKey > > segments( _segments.size() );
	for( int i=0 ; i<segments.size() ; i++ )
	{
		segments[i].resize( _segments[i].size() );
		for( int j=0 ; j<_segments[i].size() ; j++ ) segments[i][j] = edgeKey( _segments[i][j] );
	}
	return segments;
}

template< typename VertexKey , typename EdgeKey >
typename Graph< VertexKey , EdgeKey >::VertexIndex Graph< VertexKey , EdgeKey >::_vertexIndex( VertexKey vKey ) const
{
	auto iter = _vertexToIndex.find(vKey);
	if( iter==_vertexToIndex.end() ) THROW( "vertex key does not exist" );
	return iter->second;
}

template< typename VertexKey , typename EdgeKey >
typename Graph< VertexKey , EdgeKey >::EdgeIndex Graph< VertexKey , EdgeKey >::_edgeIndex( EdgeKey eKey ) const
{
	auto iter = _edgeToIndex.find(eKey);
	if( iter==_edgeToIndex.end() ) THROW( "edge key does not exist" );
	return iter->second;
}

template< typename VertexKey , typename EdgeKey >
bool Graph< VertexKey , EdgeKey >::isConnected( void ) const
{
	UnionFind unionFind( _vertexToIndex.size() );
	for( EdgeIndex i=0 ; i<_edges.size() ; i++ )
	{
		VertexIndex v1 = _edges[i].v1 , v2 = _edges[i].v2;
		if( unionFind.find(v1)!=unionFind.find(v2) ) unionFind.merge(v1,v2);
	}
	size_t root = unionFind.find( 0 );
	for( size_t v=0 ; v<_vertexToIndex.size() ; v++ ) if( unionFind.find(v)!=root ) return false;
	return true;
}

template< typename VertexKey , typename EdgeKey >
bool Graph< VertexKey , EdgeKey >::isConnected( std::function< bool ( EdgeKey ) > ValidEdgeFunction ) const
{
	UnionFind unionFind( _vertexToIndex.size() );
	for( EdgeIndex i=0 ; i<_edges.size() ; i++ ) if( ValidEdgeFunction( edgeKey( i ) ) )
	{
		VertexIndex v1 = _edges[i].v1 , v2 = _edges[i].v2;
		if( unionFind.find(v1)!=unionFind.find(v2) ) unionFind.merge(v1,v2);
	}
	size_t root = unionFind.find( 0 );
	for( size_t v=0 ; v<_vertexToIndex.size() ; v++ ) if( unionFind.find(v)!=root ) return false;
	return true;
}

template< typename VertexKey , typename EdgeKey >
Graph< VertexKey , EdgeKey > Graph< VertexKey , EdgeKey >::spanningTree( bool complement ) const
{
	std::vector< bool > keptEdges( _edges.size() , false );
	UnionFind unionFind( _vertexToIndex.size() );
	for( size_t i=0 ; i<_edges.size() ; i++ )
	{
		VertexIndex v1 = _edges[i].v1 , v2 = _edges[i].v2;
		if( unionFind.find(v1)!=unionFind.find(v2) ) unionFind.merge(v1,v2) , keptEdges[i] = true;
	}

	size_t count = 0;
	for( int i=0 ; i<_edges.size() ; i++ ) if( keptEdges[i] ) count++;

	if( complement )
	{
		std::vector< Edge< VertexKey , EdgeKey > > stEdges( _edges.size()-count );
		count = 0;
		for( size_t i=0 ; i<_edges.size() ; i++ ) if( !keptEdges[i] )
		{
			stEdges[count].e = _indexToEdge[ _edges[i].e ];
			stEdges[count].v1 = _indexToVertex[ _edges[i].v1 ];
			stEdges[count].v2 = _indexToVertex[ _edges[i].v2 ];
			count++;
		}
		return Graph< VertexKey , EdgeKey >( stEdges );
	}
	else
	{
		std::vector< Edge< VertexKey , EdgeKey > > stEdges( count );
		count = 0;
		for( size_t i=0 ; i<_edges.size() ; i++ ) if( keptEdges[i] )
		{
			stEdges[count].e = _indexToEdge[ _edges[i].e ];
			stEdges[count].v1 = _indexToVertex[ _edges[i].v1 ];
			stEdges[count].v2 = _indexToVertex[ _edges[i].v2 ];
			count++;
		}
		return Graph< VertexKey , EdgeKey >( stEdges );
	}
}

///////////////////
// WeightedGraph //
///////////////////
template< typename VertexKey , typename EdgeKey , typename Real >
WeightedGraph< VertexKey , EdgeKey , Real >::WeightedGraph( const std::vector< WeightedEdge< VertexKey , EdgeKey , Real > > &edges ) : Graph( edges )
{
	_weights.resize( edges.size() );
	for( size_t e=0 ; e<edges.size() ; e++ ) _weights[e] = edges[e].weight;
}

#ifdef NEW_GRAPH
template< typename VertexKey , typename EdgeKey , typename Real >
typename WeightedGraph< VertexKey , EdgeKey , Real >::DijkstraData WeightedGraph< VertexKey , EdgeKey , Real >::_dijkstra( VertexIndex source , VertexIndex target , std::function< bool (EdgeKey) > ValidEdgeFunction ) const
{
	DijkstraData dData;
#ifdef DEBUG_GRAPH
	VectorWrapper< bool , VertexIndex > processedVertices( vertices() , false );
#else // !DEBUG_GRAPH
	std::vector< bool > processedVertices( vertices() , false );
#endif // DEBUG_GRAPH
	dData._source = source;
	dData._vertices.resize( vertices() );
	struct S
	{
		VertexIndex vertex;
		EdgeIndex previousEdge;
		Real dist;
		S( VertexIndex v , EdgeIndex e=-1 , Real d=std::numeric_limits< Real >::infinity() ) : vertex(v) , previousEdge(e) , dist(d) {}
	};

	auto SCompare = []( S s1 , S s2 ){ return s1.dist>s2.dist; };
	std::priority_queue< S , std::vector< S > , decltype( SCompare ) > q( SCompare );
	q.push( S( source , -1 , 0 ) );
	while( q.size() )
	{
		S s = q.top() ; q.pop();
		VertexIndex v = s.vertex;
		if( !processedVertices[v] )
		{
			dData._vertices[v] = DijkstraVertex< VertexIndex , EdgeIndex , Real >( v , s.previousEdge , s.dist );
			if( v==target ) return dData;
			const std::vector< EdgeIndex > &neighbors = _neighborEdges[v];
			for( size_t i=0 ; i<neighbors.size() ; i++ )
			{
				EdgeIndex e = neighbors[i];
				if( ValidEdgeFunction( edgeKey(e) ) )
				{
					VertexIndex o = _edges[e].opposite(v);
					if( !processedVertices[o] ) q.push( S( o , e , s.dist + _weights[e] ) );
				}
			}
			processedVertices[v] = true;
		}
	}
	return dData;
}
#else // !NEW_GRAPH
#ifdef DEBUG_GRAPH
template< typename VertexKey , typename EdgeKey , typename Real >
VectorWrapper< DijkstraVertex< typename Graph< VertexKey , EdgeKey >::VertexIndex , typename Graph< VertexKey , EdgeKey >::EdgeIndex , Real > , typename Graph< VertexKey , EdgeKey >::VertexIndex > WeightedGraph< VertexKey , EdgeKey , Real >::_dijkstra( VertexIndex source , VertexIndex target , std::function< bool (EdgeKey) > ValidEdgeFunction ) const
#else // !DEBUG_GRAPH
template< typename VertexKey , typename EdgeKey , typename Real >
std::vector< DijkstraVertex< typename Graph< VertexKey , EdgeKey >::VertexIndex , typename Graph< VertexKey , EdgeKey >::EdgeIndex , Real > > WeightedGraph< VertexKey , EdgeKey , Real >::_dijkstra( VertexIndex source , VertexIndex target , std::function< bool (EdgeKey) > ValidEdgeFunction ) const
#endif // DEBUG_GRAPH
{
#ifdef DEBUG_GRAPH
	VectorWrapper< bool , VertexIndex > processedVertices( vertices() , false );
	VectorWrapper< DijkstraVertex< VertexIndex , EdgeIndex , Real > , VertexIndex > dVertices( vertices() );
#else // !DEBUG_GRAPH
	std::vector< bool > processedVertices( vertices() , false );
	std::vector< DijkstraVertex< VertexIndex , EdgeIndex , Real > > dVertices( vertices() );
#endif // DEBUG_GRAPH
	struct S
	{
		VertexIndex vertex;
		EdgeIndex previousEdge;
		Real dist;
		S( VertexIndex v , EdgeIndex e=-1 , Real d=std::numeric_limits< Real >::infinity() ) : vertex(v) , previousEdge(e) , dist(d) {}
	};

	auto SCompare = []( S s1 , S s2 ){ return s1.dist>s2.dist; };
	std::priority_queue< S , std::vector< S > , decltype( SCompare ) > q( SCompare );
	q.push( S( source , -1 , 0 ) );
	while( q.size() )
	{
		S s = q.top() ; q.pop();
		VertexIndex v = s.vertex;
		if( !processedVertices[v] )
		{
			dVertices[v] = DijkstraVertex< VertexIndex , EdgeIndex , Real >( v , s.previousEdge , s.dist );
			if( v==target ) return dVertices;
			const std::vector< EdgeIndex > &neighbors = _neighborEdges[v];
			for( size_t i=0 ; i<neighbors.size() ; i++ )
			{
				EdgeIndex e = neighbors[i];
				if( ValidEdgeFunction( edgeKey(e) ) )
				{
					VertexIndex o = _edges[e].opposite(v);
					if( !processedVertices[o] ) q.push( S( o , e , s.dist + _weights[e] ) );
				}
			}
			processedVertices[v] = true;
		}
	}
	return dVertices;
}
#endif // NEW_GRAPH

template< typename VertexKey , typename EdgeKey , typename Real >
std::vector< EdgeKey > WeightedGraph< VertexKey , EdgeKey , Real >::shortestPath( VertexKey source , VertexKey target , std::function< bool (EdgeKey) > ValidEdgeFunction ) const
{
#ifdef NEW_GRAPH
	DijkstraData dData = _dijkstra( _vertexIndex( source ) , _vertexIndex( target ) , ValidEdgeFunction );
	return shortestPath( dData , target );
#else // !NEW_GRAPH
	std::vector< EdgeIndex > _shortest;
	VertexIndex _source = _vertexIndex( source ) , _target( _vertexIndex( target ) );
#ifdef DEBUG_GRAPH
	VectorWrapper< DijkstraVertex< VertexIndex , EdgeIndex , Real > , VertexIndex > dVertices = _dijkstra( _source , _target , ValidEdgeFunction );
#else // !DEBUG_GRAPH
	std::vector< DijkstraVertex< VertexIndex , EdgeIndex , Real > > dVertices = _dijkstra( _source , _target , ValidEdgeFunction );
#endif // DEBUG_GRAPH
	
	VertexIndex v = _target;
	while( dVertices[v].previousEdge!=-1 )
	{
		_shortest.push_back( dVertices[v].previousEdge );
		v = _edges[ dVertices[v].previousEdge ].opposite( v );
	}
	if( v!=_source ) THROW( "could not make it to source" );

	std::vector< EdgeKey > shortest( _shortest.size() );
	for( size_t i=0 ; i<_shortest.size() ; i++ ) shortest[i] = edgeKey( _shortest[ _shortest.size()-1-i ] );

	return shortest;
#endif // NEW_GRAPH
}

#ifdef NEW_GRAPH
template< typename VertexKey , typename EdgeKey , typename Real >
std::vector< typename WeightedGraph< VertexKey , EdgeKey , Real >::EdgeIndex > WeightedGraph< VertexKey , EdgeKey , Real >::_shortestPath( const DijkstraData &dData , VertexIndex target ) const
{
	std::vector< EdgeIndex > shortest;

	VertexIndex v = target;
	while( dData._vertices[v].previousEdge!=-1 )
	{
		shortest.push_back( dData._vertices[v].previousEdge );
		v = _edges[ dData._vertices[v].previousEdge ].opposite( v );
	}
	if( v!=dData._source ) THROW( "could not make it to source" );
	return shortest;
}

template< typename VertexKey , typename EdgeKey , typename Real >
std::vector< EdgeKey > WeightedGraph< VertexKey , EdgeKey , Real >::shortestPath( const DijkstraData &dData , VertexKey target ) const
{
	std::vector< EdgeIndex > _shortest = _shortestPath( dData , _vertexIndex( target ) );
	std::vector< EdgeKey > shortest( _shortest.size() );
	for( size_t i=0 ; i<_shortest.size() ; i++ ) shortest[i] = edgeKey( _shortest[ _shortest.size()-1-i ] );
	return shortest;
}
#endif // NEW_GRAPH


#ifdef NEW_GRAPH
template< typename VertexKey , typename EdgeKey , typename Real >
typename WeightedGraph< VertexKey , EdgeKey , Real >::DijkstraData WeightedGraph< VertexKey , EdgeKey , Real >::dijkstra( VertexKey source , VertexKey target , std::function< bool (EdgeKey) > ValidEdgeFunction ) const
{
	return _dijkstra( _vertexIndex( source ) , _vertexIndex( target ) , ValidEdgeFunction );
}

template< typename VertexKey , typename EdgeKey , typename Real >
typename WeightedGraph< VertexKey , EdgeKey , Real >::DijkstraData WeightedGraph< VertexKey , EdgeKey , Real >::dijkstra( VertexKey source , std::function< bool (EdgeKey) > ValidEdgeFunction ) const
{
	return _dijkstra( _vertexIndex( source ) , -1 , ValidEdgeFunction );
}
#else // !NEW_GRAPH
template< typename VertexKey , typename EdgeKey , typename Real >
std::vector< DijkstraVertex< VertexKey , EdgeKey , Real > > WeightedGraph< VertexKey , EdgeKey , Real >::dijkstra( VertexKey source , VertexKey target ) const
{
#ifdef DEBUG_GRAPH
	VectorWrapper< DijkstraVertex< VertexIndex , EdgeIndex , Real > , VertexIndex > _dVertices = _dijkstra( _vertexIndex( source ) , _vertexIndex( target ) );
#else // !DEBUG_GRAPH
	std::vector< DijkstraVertex< VertexIndex , EdgeIndex , Real > > _dVertices = _dijkstra( _vertexIndex( source ) , _vertexIndex( target ) );
#endif // DEBUG_GRAPH
	std::vector< DijkstraVertex< VertexKey , EdgeKey , Real > > dVertices( _dVertices.size() );
	for( int i=0 ; i<_dVertices.size() ; i++ )
	{
		dVertices[i].distanceToSource = _dVertices[i].distanceToSource;
		if( _dVertices[i].vertex!=-1 ) dVertices[i].vertex = vertex( _dVertices[i].vertex );
		if( _dVertices[i].previousEdge!=-1 ) dVertices[i].previousEdge = edge( _dVertices[i].previousEdge );
	}
	return dVertices;
}

template< typename VertexKey , typename EdgeKey , typename Real >
std::vector< DijkstraVertex< VertexKey , EdgeKey , Real > > WeightedGraph< VertexKey , EdgeKey , Real >::dijkstra( VertexKey source ) const
{
#ifdef DEBUG_GRAPH
	VectorWrapper< DijkstraVertex< VertexIndex , EdgeIndex , Real > , VertexIndex > _dVertices = _dijkstra( _vertexIndex( source ) , -1 );
#else // !DEBUG_GRAPH
	std::vector< DijkstraVertex< VertexIndex , EdgeIndex , Real > > _dVertices = _dijkstra( _vertexIndex( source ) , -1 );
#endif // DEBUG_GRAPH
	std::vector< DijkstraVertex< VertexKey , EdgeKey , Real > > dVertices( _dVertices.size() );
	for( int i=0 ; i<_dVertices.size() ; i++ )
	{
		dVertices[i].distanceToSource = _dVertices[i].distanceToSource;
		if( _dVertices[i].vertex!=-1 ) dVertices[i].vertex = vertex( _dVertices[i].vertex );
		if( _dVertices[i].previousEdge!=-1 ) dVertices[i].previousEdge = edge( _dVertices[i].previousEdge );
	}
	return dVertices;
}
#endif // NEW_GRAPH

template< typename VertexKey , typename EdgeKey , typename Real >
std::vector< EdgeKey > WeightedGraph< VertexKey , EdgeKey , Real >::shortestCycle( std::function< bool (EdgeKey) > ValidEdgeFunction ) const
{
	Real cycleLength = std::numeric_limits< Real >::infinity();
	std::vector< EdgeIndex > _cycle;
	for( EdgeIndex e=0 ; e<edges() ; e++ ) if( ValidEdgeFunction( edgeKey(e) ) )
	{
		WeightedGraph< VertexKey , EdgeKey , Real > subGraph = *this;
		subGraph._remove( e );
#ifdef NEW_GRAPH
		DijkstraData dData = subGraph._dijkstra( _edges[e].v1 , _edges[e].v2 , ValidEdgeFunction );
		Real length = _weights[e]  + dData._distanceToSource( _edges[e].v2 );
		if( length<cycleLength )
		{
			_cycle = _shortestPath( dData , _edges[e].v2 );
			_cycle.push_back( e );
		}
#else // !NEW_GRAPH
#ifdef DEBUG_GRAPH
		VectorWrapper< DijkstraVertex< VertexIndex , EdgeIndex , Real > , VertexIndex > _dVertices = subGraph._dijkstra( _edges[e].v1 , _edges[e].v2 , ValidEdgeFunction );
#else // !DEBUG_GRAPH
		std::vector< DijkstraVertex< VertexIndex , EdgeIndex , Real > > _dVertices = subGraph._dijkstra( _edges[e].v1 , _edges[e].v2 , ValidEdgeFunction );
#endif // DEBUG_GRAPH
		Real length = _weights[e] + _dVertices[ _edges[e].v2 ].distanceToSource;
		if( length<cycleLength )
		{
			cycleLength = length;
			_cycle.resize(0);
			VertexIndex v = _edges[e].v2;
			while( _dVertices[v].previousEdge!=-1 )
			{
				_cycle.push_back( _dVertices[v].previousEdge );
				v = _edges[ _dVertices[v].previousEdge ].opposite(v);
			}
			if( v!=_edges[e].v1 ) THROW( "did not close loop" );
			_cycle.push_back( e );
		}
#endif // NEW_GRAPH
	}
	std::vector< EdgeKey > cycle( _cycle.size() );
	for( int i=0 ; i<cycle.size() ; i++ ) cycle[i] = edgeKey( _cycle[ cycle.size()-1-i ] );
	return cycle;
}

template< typename VertexKey , typename EdgeKey , typename Real >
WeightedGraph< VertexKey , EdgeKey , Real > WeightedGraph< VertexKey , EdgeKey , Real >::minimalSpanningTree( bool complement ) const
{
#ifdef DEBUG_GRAPH
	VectorWrapper< bool , EdgeIndex > keptEdges( edges() , false );
#else // !DEBUG_GRAPH
	std::vector< bool > keptEdges( edges() , false );
#endif // DEBUG_GRAPH
	UnionFind unionFind( vertices() );
	std::vector< EdgeIndex > indices( edges() );
	for( int i=0 ; i<edges() ; i++ ) indices[i] = i;
	std::sort( indices.begin() , indices.end() , [&]( EdgeIndex e1 , EdgeIndex e2 ){ return _weights[e1]<_weights[e2]; } );

	for( size_t i=0 ; i<_edges.size() ; i++ )
	{
		VertexIndex v1 = _edges[ indices[i] ].v1 , v2 = _edges[ indices[i] ].v2;
		if( unionFind.find((size_t)v1)!=unionFind.find((size_t)v2) ) unionFind.merge((size_t)v1,(size_t)v2) , keptEdges[ indices[i] ] = true;
	}

	size_t count = 0;
	for( int i=0 ; i<_edges.size() ; i++ ) if( keptEdges[i] ) count++;
	if( complement )
	{
		std::vector< WeightedEdge< VertexKey , EdgeKey , Real > > mstEdges( _edges.size()-count );
		count = 0;
		for( size_t i=0 ; i<_edges.size() ; i++ ) if( !keptEdges[i] )
		{
			mstEdges[count].e = _indexToEdge[ _edges[i].e ];
			mstEdges[count].v1 = _indexToVertex[ _edges[i].v1 ];
			mstEdges[count].v2 = _indexToVertex[ _edges[i].v2 ];
			mstEdges[count].weight = _weights[i];
			count++;
		}
		return WeightedGraph< VertexKey , EdgeKey , Real >( mstEdges );
	}
	else
	{
		std::vector< WeightedEdge< VertexKey , EdgeKey , Real > > mstEdges( count );
		count = 0;
		for( size_t i=0 ; i<_edges.size() ; i++ ) if( keptEdges[i] )
		{
			mstEdges[count].e = _indexToEdge[ _edges[i].e ];
			mstEdges[count].v1 = _indexToVertex[ _edges[i].v1 ];
			mstEdges[count].v2 = _indexToVertex[ _edges[i].v2 ];
			mstEdges[count].weight = _weights[i];
			count++;
		}
		return WeightedGraph< VertexKey , EdgeKey , Real >( mstEdges );
	}
}