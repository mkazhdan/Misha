/*
Copyright (c) 2019, Michael Kazhdan and Fabian Prada
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

/////////////////////
// System matrices //
/////////////////////
// [WARNING] This should be moved into Geometry.[h/inl]
template< unsigned int I > typename std::enable_if< I!=0 , unsigned int >::type Factorial( void ){ return Factorial< I-1 >() * I; }
template< unsigned int I > typename std::enable_if< I==0 , unsigned int >::type Factorial( void ){ return 1; }

template< typename Real , unsigned int Dim , typename ... OtherPoints >
Real Volume( Point< Real , Dim > point , OtherPoints ... otherPoints )
{
	const size_t NumPoints = sizeof...(OtherPoints)+1;
	SquareMatrix< Real , NumPoints > massMatrix;
	Point< Real , Dim > _otherPoints[] = { otherPoints ... };
	massMatrix(0,0) = Point< Real , Dim >::Dot( point , point );
	for( int i=1 ; i<NumPoints ; i++ ) massMatrix(0,i) = massMatrix(i,0) = Point< Real , Dim >::Dot( point , _otherPoints[i-1] );
	for( int i=1 ; i<NumPoints ; i++ ) for( int j=1 ; j<NumPoints ; j++ ) massMatrix(i,j) = Point< Real , Dim >::Dot( _otherPoints[i-1] , _otherPoints[j-1] );
	return (Real)sqrt( std::max< Real >( 0 , massMatrix.determinant() ) ) / Factorial< NumPoints-1 >();
}

template< typename Real , unsigned int Dim >
SquareMatrix< Real , 3 > TriangleMassMatrix( const Point< Real , Dim > vertices[] )
{
	SquareMatrix< Real , 3 > mass;
	Real area = Volume( vertices[1]-vertices[0] , vertices[2]-vertices[0] );
	mass(1,0) = mass(0,1) = mass(2,0) = mass(0,1) = mass(0,2) = mass(1,2) = mass(2,1) = area / 12;
	mass(0,0) = mass(1,1) = mass(2,2) = area / 6;
	return mass;
}
template< typename Real >
SquareMatrix< Real , 3 > TriangleTutteMatrix( void )
{
	SquareMatrix< Real , 3 > stiffness;
	stiffness(0,0) = stiffness(1,1) = stiffness(2,2) = 1;
	stiffness(0,1) = stiffness(1,0) = stiffness(2,0) = stiffness(0,2) = stiffness(1,2) = stiffness(2,1) = (Real)-0.5;
	return stiffness;
}
template< typename Real , unsigned int Dim >
SquareMatrix< Real , 3 > TriangleStiffnessMatrix( const Point< Real , Dim > vertices[] )
{
	SquareMatrix< Real , 3 > stiffness;
	Real area = Volume( vertices[1]-vertices[0] , vertices[2]-vertices[0] );
	for( int i=0 ; i<3 ; i++ )
	{
		int i0 = (i+1)%3 , i1 = (i+2)%3;
		Real dot = Point< Real , Dim >::Dot( vertices[i0]-vertices[i] , vertices[i1]-vertices[i] );
		stiffness(i0,i1) = stiffness(i1,i0) = - dot / area;
	}
	stiffness(0,0) = - ( stiffness(0,1) + stiffness(0,2) );
	stiffness(1,1) = - ( stiffness(1,2) + stiffness(1,0) );
	stiffness(2,2) = - ( stiffness(2,0) + stiffness(2,1) );
	return stiffness;
}

template< typename Real >
void MakeLumped( SparseMatrix< Real , int > &m )
{
	ThreadPool::ParallelFor
		(
			0 , m.rows ,
			[&]( unsigned int , size_t i )
			{
				Real sum = 0;
				for( int j=0 ; j<m.rowSizes[i] ; j++ )
				{
					sum += m[i][j].Value;
					if( m[i][j].N!=i ) m[i][j].Value = 0;
				}
				for( int j=0 ; j<m.rowSizes[i] ; j++ ) if( m[i][j].N==i ) m[i][j].Value = sum;
			}
		);
}

template< typename Real >
void ReorderMatrixEntries( SparseMatrix< Real , int > &M )
{
	ThreadPool::ParallelFor
		(
			0 , M.rows ,
			[&]( unsigned int , size_t i )
			{
				MatrixEntry< Real , int > *begin = &M[i][0];
				MatrixEntry< Real , int > *end = begin + M.rowSizes[i];
				std::sort( begin , end , []( MatrixEntry< Real , int > &e1 , MatrixEntry< Real , int > &e2 ){ return e1.N<e2.N; } );
			}
		);
}

////////////////////////////
// TriangleMesh::HalfEdge //
////////////////////////////
template< typename VertexKey >
VertexKey TriangleMesh< VertexKey >::HalfEdge::opposite( VertexKey v ) const
{
	if     (  first==v ) return second;
	else if( second==v ) return first;
	else THROW( "could not find vertex in half edge" );
}

//////////////////
// TriangleMesh //
//////////////////
template< typename VertexKey >
void TriangleMesh< VertexKey >::reset( void )
{
	_vertexToHalfEdgeIndex.clear();
	_triangleToHalfEdgeIndex.resize( triangles.size() );
	_halfEdgeToHalfEdgeIndex.reserve( triangles.size()*3 );

	for( TriangleIndex t=0 ; t<triangles.size() ; t++ ) for( int i=0 ; i<3 ; i++ )
	{
		// Create the half through the two other vertices
		HalfEdge he( triangles[t][(i+1)%3] , triangles[t][(i+2)%3] );

		// Check that we haven't added the edge already
		if( _halfEdgeToHalfEdgeIndex.find(he)!=_halfEdgeToHalfEdgeIndex.end() ) ERROR_OUT( "duplicate edge" );

		// Set the indices of the half-edge, vertex, and triangle
		HalfEdgeIndex halfEdgeIndex = (HalfEdgeIndex)( t*3+i );
		_halfEdgeToHalfEdgeIndex[he] = halfEdgeIndex;
		_vertexToHalfEdgeIndex[ triangles[t][i] ] = (HalfEdgeIndex)( t*3+(i+2)%3 );
		_triangleToHalfEdgeIndex[t] = halfEdgeIndex;

		// Create the half-edge data
		_HalfEdgeData halfEdgeData( he , HalfEdge::OppositeCornerInfo( t , i ) );

		halfEdgeData.previousHalfEdgeIndex = (HalfEdgeIndex)( t*3+(i+2)%3 );
		halfEdgeData.nextHalfEdgeIndex     = (HalfEdgeIndex)( t*3+(i+1)%3 );

		_halfEdgeData.push_back( halfEdgeData );
	}

	for( auto iter=_halfEdgeToHalfEdgeIndex.begin() ; iter!=_halfEdgeToHalfEdgeIndex.end() ; iter++ )
	{
		auto oIter = _halfEdgeToHalfEdgeIndex.find( iter->first.opposite() );
		if( oIter!=_halfEdgeToHalfEdgeIndex.end() ) _halfEdgeData[ iter->second ].oppositeHalfEdgeIndex = oIter->second;
	}
}

template< typename VertexKey >
size_t TriangleMesh< VertexKey >::vertices( void ) const { return _vertexToHalfEdgeIndex.size(); }

template< typename VertexKey >
size_t TriangleMesh< VertexKey >::halfEdges( void ) const { return _halfEdgeData.size(); }

template< typename VertexKey >
size_t TriangleMesh< VertexKey >::boundaryHalfEdges( void ) const
{
	size_t sz = 0;
	for( HalfEdgeIndex e=0 ; e<_halfEdgeData.size() ; e++ ) if( _halfEdgeData[e].oppositeHalfEdgeIndex==-1 ) sz++;
	return sz;
}

template< typename VertexKey >
size_t TriangleMesh< VertexKey >::edges( void ) const { return ( halfEdges()+boundaryHalfEdges() ) / 2; }

template< typename VertexKey >
int TriangleMesh< VertexKey >::eulerCharacteristic( void ) const
{
	return (int)( vertices() - edges() + triangles.size() );
}

template< typename VertexKey >
size_t TriangleMesh< VertexKey >::genus( void ) const { return ( 2-eulerCharacteristic() ) / 2; }

template< typename VertexKey >
bool TriangleMesh< VertexKey >::isWaterTight( void ) const
{
	for( auto i=_halfEdgeToHalfEdgeIndex.begin() ; i!=_halfEdgeToHalfEdgeIndex.end() ; i++ ) if( _halfEdgeToHalfEdgeIndex.find( i->first.opposite() )==_halfEdgeToHalfEdgeIndex.end() ) return false;
	return true;
}

template< typename VertexKey >
typename TriangleMesh< VertexKey >::EdgeIndex TriangleMesh< VertexKey >::edgeIndex( HalfEdgeIndex eIndex ) const
{
	if( _halfEdgeData[eIndex].oppositeHalfEdgeIndex==-1 || _halfEdgeData[eIndex].halfEdge.v1<_halfEdgeData[eIndex].halfEdge.v2 ) return (size_t)eIndex;
	else return (size_t)_halfEdgeData[eIndex].oppositeHalfEdgeIndex;
}

template< typename VertexKey >
Graph< VertexKey , typename TriangleMesh< VertexKey >::EdgeIndex > TriangleMesh< VertexKey >::graph( void ) const
{
	std::vector< Edge< VertexKey , EdgeIndex > > edges;
	edges.reserve( _halfEdgeData.size() / 2 );
	for( HalfEdgeIndex i=0 ; i<_halfEdgeData.size() ; i++ ) if( _halfEdgeData[i].oppositeHalfEdgeIndex==-1 || _halfEdgeData[i].halfEdge.v1<_halfEdgeData[i].halfEdge.v2 )
	{
		Edge< VertexKey , EdgeIndex > edge;
		edge.e = (size_t)i;
		edge.v1 = _halfEdgeData[i].halfEdge.v1;
		edge.v2 = _halfEdgeData[i].halfEdge.v2;
		edges.push_back( edge );
	}
	return Graph< VertexKey , EdgeIndex >( edges );
}

template< typename VertexKey >
Graph< typename TriangleMesh< VertexKey >::TriangleIndex , typename TriangleMesh< VertexKey >::EdgeIndex > TriangleMesh< VertexKey >::dualGraph( void ) const
{
	std::vector< Edge< TriangleIndex , EdgeIndex > > edges;
	edges.reserve( _halfEdgeData.size() / 2 );
	for( HalfEdgeIndex i=0 ; i<_halfEdgeData.size() ; i++ ) if( _halfEdgeData[i].oppositeHalfEdgeIndex!=-1 && _halfEdgeData[i].halfEdge.v1<_halfEdgeData[i].halfEdge.v2 )
	{
		Edge< TriangleIndex , EdgeIndex > edge;
		edge.e = (size_t)i;
		edge.v1 = _halfEdgeData[               i                        ].halfEdge.oppositeCornerInfo.triangle;
		edge.v2 = _halfEdgeData[ _halfEdgeData[i].oppositeHalfEdgeIndex ].halfEdge.oppositeCornerInfo.triangle;
		edges.push_back( edge );
	}
	return Graph< TriangleIndex , EdgeIndex >( edges );
}

template< typename VertexKey >
template< typename Real >
WeightedGraph< VertexKey , typename TriangleMesh< VertexKey >::EdgeIndex , Real > TriangleMesh< VertexKey >::graph( const std::vector< Point3D< Real > > &vertices ) const
{
	std::vector< WeightedEdge< VertexKey , EdgeIndex , Real > > edges;
	edges.reserve( _halfEdgeData.size()/ 2  );
	for( HalfEdgeIndex i=0 ; i<_halfEdgeData.size() ; i++ ) if( _halfEdgeData[i].oppositeHalfEdgeIndex==-1 || _halfEdgeData[i].halfEdge.v1<_halfEdgeData[i].halfEdge.v2 )
	{
		WeightedEdge< VertexKey , EdgeIndex , Real > edge;
		edge.e = (size_t)i;
		edge.v1 = _halfEdgeData[i].halfEdge.v1;
		edge.v2 = _halfEdgeData[i].halfEdge.v2;
		edge.weight = (Real)Point3D< Real >::Length( vertices[ (size_t)edge.v1 ] - vertices[ (size_t)edge.v2] );
		edges.push_back( edge );
	}
	return WeightedGraph< VertexKey , EdgeIndex , Real >( edges );
}

template< typename VertexKey >
template< typename Real >
WeightedGraph< typename TriangleMesh< VertexKey >::TriangleIndex , typename TriangleMesh< VertexKey >::EdgeIndex , Real > TriangleMesh< VertexKey >::dualGraph( const std::vector< Point3D< Real > > &vertices ) const
{
	std::vector< WeightedEdge< TriangleIndex , EdgeIndex , Real > > edges;
	edges.reserve( _halfEdgeData.size() / 2 );
	for( HalfEdgeIndex i=0 ; i<_halfEdgeData.size() ; i++ ) if( _halfEdgeData[i].oppositeHalfEdgeIndex!=-1 && _halfEdgeData[i].halfEdge.v1<_halfEdgeData[i].halfEdge.v2 )
	{
		WeightedEdge< TriangleIndex , EdgeIndex , Real > edge;
		edge.e = (size_t)i;
		edge.v1 = _halfEdgeData[               i                        ].oppositeCornerInfo.triangle;
		edge.v2 = _halfEdgeData[ _halfEdgeData[i].oppositeHalfEdgeIndex ].oppositeCornerInfo.triangle;
		std::vector< HalfEdgeIndex > triangle1 = face( edge.v1 );
		std::vector< HalfEdgeIndex > triangle2 = face( edge.v2 );
		Point3D< Real > c1 , c2;
		for( int j=0 ; j<3 ; j++ ) c1 += vertices[ (size_t)_halfEdgeData[ triangle1[j] ].halfEdge.v1 ] , c2 += vertices[ (size_t)_halfEdgeData[ triangle2[j] ].halfEdge.v1 ];
		c1 /= 3 , c2 /= 3;
		edge.weight = (Real)Point3D< Real >::Length( c1 - c2 );
		edges.push_back( edge );
	}
	return WeightedGraph< TriangleIndex , EdgeIndex , Real >( edges );
}

template< typename VertexKey >
std::vector< typename TriangleMesh< VertexKey >::HalfEdgeIndex > TriangleMesh< VertexKey >::face( TriangleIndex t ) const
{
	std::vector< HalfEdgeIndex > f(3);
	f[0] = _triangleToHalfEdgeIndex[t];
	f[1] = _halfEdgeData[ f[0] ].nextHalfEdgeIndex;
	f[2] = _halfEdgeData[ f[1] ].nextHalfEdgeIndex;
	if( _halfEdgeData[ f[2] ].nextHalfEdgeIndex!=f[0] ) ERROR_OUT( "could not form a triangle" );
	return f;
}

template< typename VertexKey >
bool TriangleMesh< VertexKey >::isBoundaryVertex( VertexKey v ) const
{
	HalfEdgeIndex s , e;
	auto iter = _vertexToHalfEdgeIndex.find( v );
	if( iter==_vertexToHalfEdgeIndex.end() ) THROW( "could not find vertex" );
	s = e = iter->second;
	do{ e = _nextOutgoingHalfEdgeIndex( e ); }
	while( e!=-1 && e!=s );
	return e!=s;
}

template< typename VertexKey >
std::vector< typename TriangleMesh< VertexKey >::HalfEdgeIndex > TriangleMesh< VertexKey >::halfEdgeOneRing( VertexKey v ) const
{
	std::vector< HalfEdgeIndex > front;
	HalfEdgeIndex s , e;
	auto iter = _vertexToHalfEdgeIndex.find( v );
	if( iter==_vertexToHalfEdgeIndex.end() ) THROW( "could not find vertex" );
	s = e = iter->second;
	do
	{
		front.push_back( e );
		e = _nextOutgoingHalfEdgeIndex( e );
	}
	while( e!=-1 && e!=s );
	if( e==s ) return front;

	std::vector< HalfEdgeIndex > back , oneRing;
	e = s;
	do
	{
		back.push_back( e );
		e = _previousOutgoingHalfEdgeIndex( e );
	}
	while( e!=-1 );
	oneRing.reserve( back.size() + front.size() - 1 );
	for( int i=1 ; i<back.size() ; i++ ) oneRing.push_back( back[ back.size()-i ] );
	for( int i=0 ; i<front.size() ; i++ ) oneRing.push_back( front[i] );
	return oneRing;
}

template< typename VertexKey >
std::vector< VertexKey > TriangleMesh< VertexKey >::vertexOneRing( VertexKey v ) const
{
	std::vector< VertexKey > front;
	HalfEdgeIndex s , e , p;
	auto iter = _vertexToHalfEdgeIndex.find( v );
	if( iter==_vertexToHalfEdgeIndex.end() ) THROW( "could not find vertex" );
	s = e = iter->second;
	do
	{
		front.push_back( halfEdge( e ).v2 );
		p = e;
		e = _nextOutgoingHalfEdgeIndex( e );
	}
	while( e!=-1 && e!=s );
	if( e==s ) return front;

	front.push_back( halfEdge( _halfEdgeData[p].previousHalfEdgeIndex ).v1 );

	std::vector< VertexKey > back , oneRing;
	e = s;
	do
	{
		back.push_back( halfEdge( e ).v2 );
		e = _previousOutgoingHalfEdgeIndex( e );
	}
	while( e!=-1 );
	oneRing.reserve( back.size() + front.size() );
	for( int i=1 ; i<back.size() ; i++ ) oneRing.push_back( back[ back.size()-i ] );
	for( int i=0 ; i<front.size() ; i++ ) oneRing.push_back( front[i] );
	return oneRing;
}

template< typename VertexKey >
std::vector< typename TriangleMesh< VertexKey >::HalfEdgeIndex > TriangleMesh< VertexKey >::boundary( HalfEdgeIndex eIndex ) const
{
	std::vector< HalfEdgeIndex > boundary;
	if( _halfEdgeData[eIndex].oppositeHalfEdgeIndex!=-1 ) return boundary;

	HalfEdgeIndex e = eIndex;
	do
	{
		boundary.push_back( e );
		e = _nextBoundaryHalfEdgeIndex( e );
	}
	while( e!=eIndex );
	return boundary;
}

template< typename VertexKey >
std::vector< std::vector< typename TriangleMesh< VertexKey >::HalfEdgeIndex > > TriangleMesh< VertexKey >::boundaries( void ) const
{
	std::vector< std::vector< HalfEdgeIndex > > boundaries;
	std::unordered_set< HalfEdgeIndex > boundaryHalfEdgeSet;
	for( HalfEdgeIndex e=0 ; e<halfEdges() ; e++ ) if( _halfEdgeData[e].oppositeHalfEdgeIndex==-1 ) boundaryHalfEdgeSet.insert( e );
	while( !boundaryHalfEdgeSet.empty() )
	{
		auto iter = boundaryHalfEdgeSet.begin();
		boundaries.push_back( boundary( *iter ) );
		for( size_t i=0 ; i<boundaries.back().size() ; i++ ) boundaryHalfEdgeSet.erase( boundaries.back()[i] );
	}
	return boundaries;
}

template< typename VertexKey >
typename TriangleMesh< VertexKey >::HalfEdgeIndex TriangleMesh< VertexKey >::halfEdgeIndex( HalfEdge he ) const
{
	auto iter = _halfEdgeToHalfEdgeIndex.find(he);
	if( iter==_halfEdgeToHalfEdgeIndex.end() ) THROW( "edge is not in mesh" );
	else return iter->second;
}
template< typename VertexKey >
typename TriangleMesh< VertexKey >::HalfEdge TriangleMesh< VertexKey >::halfEdge( HalfEdgeIndex eIndex ) const
{
	if( eIndex>=_halfEdgeData.size() ) THROW( "edge index out of bounds: " , (size_t)eIndex , " >= " , _halfEdgeData.size() );
	return _halfEdgeData[eIndex].halfEdge; 
}

template< typename VertexKey >
std::vector< typename TriangleMesh< VertexKey >::HalfEdgeIndex > TriangleMesh< VertexKey >::orient( const std::vector< EdgeIndex > &eSegment ) const
{
	std::vector< HalfEdgeIndex > heSegment( eSegment.size() );
	if( eSegment.size()==1 ) heSegment[0] = (size_t)eSegment[0];
	else
	{
		VertexKey v;
		HalfEdge he0 = halfEdge( eSegment[0] ) , he1 = halfEdge( eSegment[1] );
		if     ( he0.v1==he1.v1 || he0.v1==he1.v2 ) v = he0.v2;
		else if( he0.v2==he1.v1 || he0.v2==he1.v2 ) v = he0.v1;
		else ERROR_OUT( "could not link first two edges" );
		for( size_t e=0 ; e<eSegment.size() ; e++ )
		{
			HalfEdge he = halfEdge( eSegment[e] );
			if( he.v1==v )
			{
				heSegment[e] = (size_t)eSegment[e];
				v = he.v2;
			}
			else if( he.v2==v )
			{
				heSegment[e] = _halfEdgeData[ (size_t)eSegment[e] ].oppositeHalfEdgeIndex;
				v = he.v1;
			}
			else THROW( "could not link edge " , (size_t)e ,  ": {" , (size_t)he.v1 , " " , (size_t)he.v2 , "}" );
		}
	}

	return heSegment;
}

template< typename VertexKey >
TriangleMesh< VertexKey > TriangleMesh< VertexKey >::split( const std::vector< HalfEdgeIndex > &path , const std::vector< VertexKey > &newVertexKeys ) const
{
	bool isLoop;
	{
		HalfEdge he0 = halfEdge( path[0] );
		HalfEdge he1 = halfEdge( path.back() );
		isLoop = he0.v1==he1.v2;
		if( !isLoop && ( !isBoundaryVertex( he0.v1 ) || !isBoundaryVertex( he1.v2 ) ) ) THROW( "path must either be a loop or have end-points on the boundary" );
	}
	if( isLoop )
	{
		if( path.size()!=newVertexKeys.size() ) THROW( "number of new vertex keys does not match the number of vertices on the loop: " , newVertexKeys.size() , " != " , path.size() );
	}
	else
	{
		if( path.size()+1!=newVertexKeys.size() ) THROW( "number of new vertex keys does not match the number of vertices on the path: " , newVertexKeys.size() , " != " , path.size() + 1 );
	}

	TriangleMesh tMesh;
	tMesh.triangles = triangles;

	for( size_t i=0 ; i<path.size()-1 ; i++ )
	{
		HalfEdge he0 = halfEdge( path[i+0] );
		HalfEdge he1 = halfEdge( path[i+1] );
		if( he0.v2!=he1.v1 ) THROW( "half edge list is not a path: " , i );
		// [WARNING] This is probably conservative and can be worked around
		if( _halfEdgeData[ path[i] ].oppositeHalfEdgeIndex==-1 ) THROW( "path edges cannot be on the boundary" );
	}
	if( _halfEdgeData[ path.back() ].oppositeHalfEdgeIndex==-1 ) THROW( "path edges cannot be on the boundary" );

	for( size_t i=0 ; i<path.size() ; i++ )
	{
		size_t i0 = (i+path.size()-1)%path.size() , i1 = i;
		HalfEdge he0 = halfEdge( path[i0] );
		HalfEdge he1 = halfEdge( path[i1] );
		HalfEdgeIndex o = _halfEdgeData[ path[i0] ].oppositeHalfEdgeIndex;
		VertexKey v = he1.v1;

		std::vector< HalfEdgeIndex > oneRing = halfEdgeOneRing( v );
		size_t begin=0 , end;
		if( isLoop || i!=0 ) for( begin=0 ; begin<oneRing.size() ; begin++ ) if( oneRing[begin]==o ) break;
		if( begin==oneRing.size() ) THROW( "Could not find begin" );
		for( end=begin ; end<begin+oneRing.size() ; end++ ) if( oneRing[ end%oneRing.size() ]==path[i1] ) break;
		if( end==begin+oneRing.size() ) THROW( "Could not find end" );
		for( size_t j=begin ; j<end ; j++ )
		{
			TriangleIndex t = _halfEdgeData[ oneRing[j%oneRing.size()] ].oppositeCornerInfo.triangle;
			for( int k=0 ; k<3 ; k++ ) if( tMesh.triangles[t][k]==v ) tMesh.triangles[t][k] = newVertexKeys[i];
		}
	}
	if( !isLoop )
	{
		size_t i0 = path.size()-1;
		HalfEdge he0 = halfEdge( path[i0] );
		HalfEdgeIndex o = _halfEdgeData[ path[i0] ].oppositeHalfEdgeIndex;
		VertexKey v = he0.v2;

		std::vector< HalfEdgeIndex > oneRing = halfEdgeOneRing( v );
		size_t begin=0 , end=oneRing.size();
		for( begin=0 ; begin<oneRing.size() ; begin++ ) if( oneRing[begin]==o ) break;
		if( begin==oneRing.size() ) THROW( "Could not find begin" );
		for( size_t j=begin ; j<end ; j++ )
		{
			TriangleIndex t = _halfEdgeData[ oneRing[j%oneRing.size()] ].oppositeCornerInfo.triangle;
			for( int k=0 ; k<3 ; k++ ) if( tMesh.triangles[t][k]==v ) tMesh.triangles[t][k] = newVertexKeys.back();
		}
	}
	tMesh.reset();
	return tMesh;
}

template< typename VertexKey >
typename TriangleMesh< VertexKey >::HalfEdgeIndex TriangleMesh< VertexKey >::_nextOutgoingHalfEdgeIndex( HalfEdgeIndex e ) const { return _halfEdgeData[ _halfEdgeData[e].previousHalfEdgeIndex ].oppositeHalfEdgeIndex; }

template< typename VertexKey >
typename TriangleMesh< VertexKey >::HalfEdgeIndex TriangleMesh< VertexKey >::_previousOutgoingHalfEdgeIndex( HalfEdgeIndex e ) const
{
	if( _halfEdgeData[e].oppositeHalfEdgeIndex==-1 ) return -1;
	else return _halfEdgeData[ _halfEdgeData[e].oppositeHalfEdgeIndex ].nextHalfEdgeIndex;
}

template< typename VertexKey >
typename TriangleMesh< VertexKey >::HalfEdgeIndex TriangleMesh< VertexKey >::_nextBoundaryHalfEdgeIndex( HalfEdgeIndex e ) const
{
	e = _halfEdgeData[e].nextHalfEdgeIndex;
	while( _halfEdgeData[e].oppositeHalfEdgeIndex!=-1 ) e = _halfEdgeData[ _halfEdgeData[e].oppositeHalfEdgeIndex ].nextHalfEdgeIndex;
	return e;
}

template< typename VertexKey >
typename TriangleMesh< VertexKey >::HalfEdgeIndex TriangleMesh< VertexKey >::_previousBoundaryHalfEdgeIndex( HalfEdgeIndex e ) const
{
	e = _halfEdgeData[e].previousHalfEdgeIndex;
	while( _halfEdgeData[e].oppositeHalfEdgeIndex!=-1 ) e = _halfEdgeData[ _halfEdgeData[e].oppositeHalfEdgeIndex ].previousHalfEdgeIndex;
	return e;
}

#if 0
template< typename VertexKey >
template< typename Real , typename TriangleMatrixFunctor >
SparseMatrix< Real , int > TriangleMesh< VertexKey >::matrix( TriangleMatrixFunctor F ) const
#else
template< typename VertexKey >
template< typename Real , unsigned int Dim >
SparseMatrix< Real , int > TriangleMesh< VertexKey >::matrix( const std::vector< Point< Real , Dim > > &vertices , SquareMatrix< Real , 3 >(*F)( const Point< Real , Dim > [] ) ) const
#endif
{
	struct Entry
	{
		Entry( void ) : row(-1) , col(-1) , value(0){}
		Entry( int r , int c , Real v ) : row(r) , col(c) , value(v){}
		int row , col;
		Real value;
	};
	SparseMatrix< Real , int > M;
	M.resize( (int)vertices.size() );
	std::vector< std::vector< Entry > > entries( ThreadPool::NumThreads() );
	ThreadPool::ParallelFor
		(
			0 , triangles.size() ,
			[&]( unsigned int thread , size_t i )
			{
				Point< Real , Dim > v[] = { vertices[ (size_t)triangles[t][0] ] , vertices[ (size_t)triangles[t][1] ] , vertices[ (size_t)triangles[t][2] ] };
				SquareMatrix< Real , 3 > m = F( v );
				for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) entries[thread].push_back( Entry( (int)(size_t)triangles[t][i] , (int)(size_t)triangles[t][j] , m(i,j) ) );
			}
		);

	for( int i=0 ; i<entries.size() ; i++ ) for( int j=0 ; j<entries[i].size() ; j++ )	M.rowSizes[ entries[i][j].row ]++;

	ThreadPool::ParallelFor
		(
			0 , M.rows ,
			[&]( unsigned int , size_t i )
			{
				int rowSize = (int)M.rowSizes[i];
				M.rowSizes[i] = 0;
				M.SetRowSize( i , rowSize );
				M.rowSizes[i] = 0;
			}
		);

	for( int i=0 ; i<entries.size() ; i++ ) for( int j=0 ; j<entries[i].size() ; j++ ) M[ entries[i][j].row ][ M.rowSizes[entries[i][j].row]++ ] = MatrixEntry< Real , int >( entries[i][j].col , entries[i][j].value );

	ThreadPool::ParallelFor
		(
			0 , M.rows ,
			[&]( unsigned int , size_t i )
			{
				std::unordered_map< int , Real > row;
				for( int j=0 ; j<M.rowSizes[i] ; j++ ) row[ M[i][j].N ] += M[i][j].Value;
				M.SetRowSize( i , (int)row.size() );
				int j=0;
				for( std::unordered_map< int , Real >::const_iterator iter=row.begin() ; iter!=row.end() ; iter++ ) M[i][j++] = MatrixEntry< Real , int >( iter->first , iter->second );
			}
		);

	return M;
}

template< typename VertexKey >
template< typename Real , unsigned int Dim >
SparseMatrix< Real , int > TriangleMesh< VertexKey >::massMatrix( const std::vector< Point< Real , Dim > > &vertices , bool lump ) const
{
	SparseMatrix< Real , int > mass = matrix( vertices , TriangleMassMatrix );
	if( lump ) MakeLumped( mass );
	ReorderMatrixEntries( mass );
	return mass;
}

template< typename VertexKey >
template< typename Real, unsigned int Dim >
SparseMatrix< Real , int > TriangleMesh< VertexKey >::stiffnessMatrix( const std::vector< Point< Real , Dim > > &vertices ) const
{
	SparseMatrix< Real , int > stiffness = matrix( vertices , TriangleStiffnessMatrix );
	ReorderMatrixEntries( stiffness );
	return stiffness;
}

template< typename VertexKey >
template< typename Real >
SparseMatrix< Real , int > TriangleMesh< VertexKey >::tutteLaplacian( void ) const
{
	struct Entry
	{
		Entry( void ) : row(-1) , col(-1) , value(0){}
		Entry( int r , int c , Real v ) : row(r) , col(c) , value(v){}
		int row , col;
		Real value;
	};
	SparseMatrix< Real , int > M;
	M.resize( (int)vertices() );
	std::vector< std::vector< Entry > > entries( ThreadPool::NumThreads );
	SquareMatrix< Real , 3 > m = TriangleTutteMatrix< Real >();

	ThreadPool::ParallelFor
		(
			0 , triangles.size() ,
			[&]( unsigned int thread , size_t t )
			{
				for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) entries[thread].push_back( Entry( (int)(size_t)triangles[t][i] , (int)(size_t)triangles[t][j] , m(i,j) ) );
			}
		);

	for( int i=0 ; i<entries.size() ; i++ ) for( int j=0 ; j<entries[i].size() ; j++ )	M.rowSizes[ entries[i][j].row ]++;

	ThreadPool::ParallelFor
		(
			0 , M.rows ,
			[&]( unsigned int , size_t i )
			{
				int rowSize = (int)M.rowSizes[i];
				M.rowSizes[i] = 0;
				M.SetRowSize( i , rowSize );
				M.rowSizes[i] = 0;
			}
		);
	for( int i=0 ; i<entries.size() ; i++ ) for( int j=0 ; j<entries[i].size() ; j++ ) M[ entries[i][j].row ][ M.rowSizes[entries[i][j].row]++ ] = MatrixEntry< Real , int >( entries[i][j].col , entries[i][j].value );

	ThreadPool::ParallelFor
		(
			0 , M.rows ,
			[&]( unsigned int , size_t i )
			{
				std::unordered_map< int , Real > row;
				for( int j=0 ; j<M.rowSizes[i] ; j++ ) row[ M[i][j].N ] += M[i][j].Value;
				M.SetRowSize( i , (int)row.size() );
				int j=0;
				for( std::unordered_map< int , Real >::const_iterator iter=row.begin() ; iter!=row.end() ; iter++ ) M[i][j++] = MatrixEntry< Real , int >( iter->first , iter->second );
			}
		);
	ReorderMatrixEntries( M );
	return M;
}
