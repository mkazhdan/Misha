/*
Copyright (c) 2015, Michael Kazhdan
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

#include <Misha/RegularGrid.h>

//////////////////
// IsoSurface2D //
//////////////////

template< typename Real >
void IsoSurface2D< Real >::Extract( const unsigned int res[2] , const Point< Real , 2 > bBox[2] , std::function< Real ( Point< Real , 2 > ) > vFunction , Real isoValue , std::vector< Point< Real , 2 > >& vertices , std::vector< EdgeIndex >& edges , bool fullCaseTable , int interpolationType )
{
	RegularGrid< Real , 2 > grid;
	grid.resize( res );
	Point< Real , 2 > d = bBox[1]-bBox[0];
	for( unsigned int i=0 ; i<res[0] ; i++ ) for( unsigned int j=0 ; j<res[1] ; j++ )
	{
		Point< Real , 3 > p = bBox[0] + Point< Real , 2 >( (Real)(i)/(res[0]-1) * d[0] , (Real)(j)/(res[1]-1) * d[1] );
		grid(i,j) = vFunction( p );
	}
	Extract( grid , isoValue , vertices , edges , fullCaseTable , interpolationType );
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) for( unsigned int j=0 ; j<2 ; j++ ) vertices[i][j] = bBox[0][j] + d[j]*vertices[i][j] / ( res[j]-1 );
}

template< typename Real >
const std::string IsoSurface2D< Real >::InterpolationNames[] = { "linear" , "quadratic" , "cubic" , "catmull-rom" };

template< typename Real >
void IsoSurface2D< Real >::Extract( const unsigned int res[2] , ConstPointer( Real ) values , Real isoValue , std::vector< Point2D< Real > > &vertices , std::vector< EdgeIndex > &edges , bool fullCaseTable , int interpolationType )
{
#ifdef NEW_ISO_SURFACE_2D
	_Extract( res , values , isoValue , vertices , edges , fullCaseTable , interpolationType );
#else // !NEW_ISO_SURFACE_2D
	std::vector< _Vertex > _vertices;
	_Extract( res , values , isoValue , _vertices , edges , fullCaseTable , interpolationType );

	vertices.resize( _vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = _vertices[i].p;
#endif // NEW_ISO_SURFACE_2D
}

template< typename Real >
void IsoSurface2D< Real >::Extract( const RegularGrid< Real , 2 > &grid , Real isoValue , std::vector< Point2D< Real > >& vertices , std::vector< EdgeIndex >& edges , bool fullCaseTable , int interpolationType )
{
	return Extract( grid.res() , grid() , isoValue , vertices , edges , fullCaseTable , interpolationType );
}

template< typename Real >
#ifdef NEW_ISO_SURFACE_2D
void IsoSurface2D< Real >::_Extract( const RegularGrid< Real , 2 > &grid , Real isoValue , std::vector< Point2D< Real > >& vertices , std::vector< EdgeIndex >& edges , bool fullCaseTable , int interpolationType )
#else // !NEW_ISO_SURFACE_2D
void IsoSurface2D< Real >::_Extract( const RegularGrid< Real , 2 > &grid , Real isoValue , std::vector< _Vertex >& vertices , std::vector< EdgeIndex >& edges , bool fullCaseTable , int interpolationType )
#endif // NEW_ISO_SURFACE_2D
{
	_Extract( grid.res() , grid() , isoValue , vertices , edges , fullCaseTable , interpolationType );
}

template< typename Real >
#ifdef NEW_ISO_SURFACE_2D
void IsoSurface2D< Real >::_Extract( const unsigned int res[2] , ConstPointer( Real ) values , Real isoValue , std::vector< Point2D< Real > >& vertices , std::vector< EdgeIndex >& edges , bool fullCaseTable , int interpolationType )
#else // !NEW_ISO_SURFACE_2D
void IsoSurface2D< Real >::_Extract( const unsigned int res[2] , ConstPointer( Real ) values , Real isoValue , std::vector< _Vertex >& vertices , std::vector< EdgeIndex >& edges , bool fullCaseTable , int interpolationType )
#endif // NEW_ISO_SURFACE_2D
{
#ifdef NEW_ISO_SURFACE_2D
	Pointer( _Vertex ) xIsoVertexMap[2];
	Pointer( _Vertex ) yIsoVertexMap;
	xIsoVertexMap[0] = NewPointer< _Vertex >( res[0] );
	xIsoVertexMap[1] = NewPointer< _Vertex >( res[0] );
	yIsoVertexMap    = NewPointer< _Vertex >( res[0] );
#else // !NEW_ISO_SURFACE_2D
	std::unordered_map< long long , int > xIsoVertexMap[2] , yIsoVertexMap;
#endif // NEW_ISO_SURFACE_2D
	Pointer( unsigned char ) flags[2];
	flags[0] = NewPointer< unsigned char >( res[0] );
	flags[1] = NewPointer< unsigned char >( res[0] );

	if( fullCaseTable ) MarchingSquares::SetFullCaseTable();
	else                MarchingSquares::SetCaseTable();

	_SetFlags    ( res[0] ,     values , isoValue , flags[0] );
	_SetXVertices( res[0] , 0 , values ,            flags[0] , isoValue , interpolationType , xIsoVertexMap[0] , vertices );

	for( int y=0 ; y<(int)res[1]-1 ; y++ )
	{
		int y0 = y&1 , y1 = (y+1)&1;
#ifdef NEW_ISO_SURFACE_2D
#else // !NEW_ISO_SURFACE_2D
		xIsoVertexMap[y1].clear() , yIsoVertexMap.clear();
#endif // NEW_ISO_SURFACE_2D
		_SetFlags    ( res[0] ,       values + (y+1)*res[0] , isoValue , flags[y1] );
		_SetXVertices( res[0] , y+1 , values + (y+1)*res[0] ,            flags[y1] , isoValue , interpolationType , xIsoVertexMap[y1] , vertices );
		_SetYVertices( res[0] , y , y>0 ? values + (y-1)*res[0] : NullPointer< Real >() , values + y*res[0] , values + (y+1)*res[0] , y+1<(int)res[1]-1 ? values + (y+2)*res[0] : NullPointer< Real >() , flags[y0] , flags[y1] , isoValue , interpolationType , yIsoVertexMap , vertices );
		_SetEdges    ( res[0] , y , values + y*res[0] , values+(y+1)*res[0] , isoValue , fullCaseTable , xIsoVertexMap[y0] , xIsoVertexMap[y1] , yIsoVertexMap , edges );
	}

	DeletePointer( flags[0] );
	DeletePointer( flags[1] );
#ifdef NEW_ISO_SURFACE_2D
	DeletePointer( xIsoVertexMap[0] );
	DeletePointer( xIsoVertexMap[1] );
	DeletePointer( yIsoVertexMap );
#endif // NEW_ISO_SURFACE_2D
}

template< typename Real >
Real IsoSurface2D< Real >::_LinearInterpolant( Real x1 , Real x2 , Real isoValue ){ return ( isoValue-x1 ) / ( x2-x1 ); }
template< typename Real >
Real IsoSurface2D< Real >::_QuadraticInterpolant( Real x0 , Real x1 , Real x2 , Real x3 , Real isoValue )
{
	// Adjust so that we are looking for a zero-crossing
	x0 -= isoValue , x1 -= isoValue , x2 -= isoValue , x3 -= isoValue;
	if( !x1 ) return 0;
	if( !x2 ) return 1;

	// Estimate the derivatives at x1 and x2
	Real dx1 = (x2-x0) / 2.f , dx2 = (x3-x1) / 2.f;
	// Solve for the quadratic polynomial:
	//		P(x) = a x^2 + b x + c 
	// such that:
	//		P(0) = x1 , P(1) = x2 , and minimizing || P'(0) - dx1 ||^2 + || P'(1) - dx2 ||^2
	//	=>  c = x1 , a = x2 - x1 - b , and minimizing || b - dx1 ||^2 + || 2*a + b - dx2 ||^2
	//	=>  c = x1 , a = x2 - x1 - b , and minimizing || b - dx1 ||^2 + || 2*x2 - 2*x1 - b - dx2 ||^2
	//	=>  c = x1 , a = x2 - x1 - b , and minimizing || b - dx1 ||^2 + || b - ( 2*x2 - 2*x1 - dx2 ) ||^2
	//	=>  c = x1 , b = ( 2*x2 - 2*x1 - dx2 + dx1 ) / 2 , a = x2 - x1 - b
	//	=>  c = x1 , b = ( x2 - x1 ) - ( dx2 - dx1 ) / 2 , a = ( dx2 - dx1 ) / 2

	double a = (dx2-dx1)/2.f , b = (dx1-dx2)/2.f + x2 - x1 , c = x1;
	if( !a )
	{
		// Solve b * x + c = 0
		return (Real)( -c / b );
	}
	else
	{
		// Solve a x^2 + b x + c = 0
		b /= a , c /= a;
		double disc = b*b - 4.*c;
		if( disc<0 ) ERROR_OUT( "Negative discriminant: " , disc );
		disc = sqrt( disc );
		double r1 = ( - b - disc ) / 2. , r2 = ( - b + disc ) / 2.;
		if( r2<0 || r1>1 ) ERROR_OUT( "Roots out of bounds: " , r1 , " " , r2 );
		if( r2>1 ) return (Real)r1;
		else       return (Real)r2;
	}
}
template< typename Real >
Real IsoSurface2D< Real >::_CubicInterpolant( Real x0 , Real x1 , Real x2 , Real x3 , Real isoValue )
{
	static bool firstTime = true;
	static Polynomial::Polynomial1D< 3 > lagrangePolynomials[4];
	if( firstTime )
	{
		Point< double , 1 > positions[] = { Point< double , 1 >( -1. ) , Point< double , 1 >( 0. ) , Point< double , 1 >( 1. ) , Point< double , 1 >( 2. ) };
		SquareMatrix< double , 4 > Einv = Polynomial::Polynomial1D< 3 >::EvaluationMatrix( positions ).inverse();
		for( int i=0 ; i<4 ; i++ )
		{
			Point< double , 4 > values;
			values[i] = 1;
			Point< double , 4 > coefficients = Einv * values;
#if 1
			for( int j=0 ; j<4 ; j++ ) lagrangePolynomials[i].coefficient( j ) = coefficients[j];
#else
			for( int j=0 ; j<4 ; j++ ) lagrangePolynomials[i][j] = coefficients[j];
#endif
		}
	}

	// Adjust so that we are looking for a zero-crossing
	x0 -= isoValue , x1 -= isoValue , x2 -= isoValue , x3 -= isoValue;
	if( !x1 ) return 0;
	if( !x2 ) return 1;

	Polynomial::Polynomial1D< 3 > p = lagrangePolynomials[0] * x0 + lagrangePolynomials[1] * x1 + lagrangePolynomials[2] * x2 + lagrangePolynomials[3] * x3;
	double roots[3] , _roots[3];
	int rootCount = Polynomial::Roots( p , roots );
	int _rootCount = 0;
	for( int i=0 ; i<rootCount ; i++ ) if( roots[i]>=0 && roots[i]<1 ) _roots[ _rootCount++ ] = roots[i];
	if     ( _rootCount==1 ) return (Real)_roots[0];
	else if( _rootCount==3 ) return (Real)_roots[1];
	else ERROR_OUT( "Unexpected number of roots: " , _rootCount );
	return 0;
}

template< typename Real >
Real IsoSurface2D< Real >::_CatmullRomInterpolant( Real x0 , Real x1 , Real x2 , Real x3 , Real isoValue )
{
	static bool firstTime = true;
	static Polynomial::Polynomial1D< 3 > blendingFunctions[4];
	if( firstTime )
	{
#if 1
		blendingFunctions[0].coefficient( 0 ) =  0.0;
		blendingFunctions[0].coefficient( 1 ) = -0.5;
		blendingFunctions[0].coefficient( 2 ) =  1.0;
		blendingFunctions[0].coefficient( 3 ) = -0.5;
		blendingFunctions[1].coefficient( 0 ) =  1.0;
		blendingFunctions[1].coefficient( 1 ) =  0.0;
		blendingFunctions[1].coefficient( 2 ) = -2.5;
		blendingFunctions[1].coefficient( 3 ) =  1.5;
		blendingFunctions[2].coefficient( 0 ) =  0.0;
		blendingFunctions[2].coefficient( 1 ) =  0.5;
		blendingFunctions[2].coefficient( 2 ) =  2.0;
		blendingFunctions[2].coefficient( 3 ) = -1.5;
		blendingFunctions[3].coefficient( 0 ) =  0.0;
		blendingFunctions[3].coefficient( 1 ) =  0.0;
		blendingFunctions[3].coefficient( 2 ) = -0.5;
		blendingFunctions[3].coefficient( 3 ) =  0.5;
#else
		blendingFunctions[0] = Polynomial::Polynomial1D< 3 >( 0.0 , -0.5 ,  1.0 , -0.5 );
		blendingFunctions[1] = Polynomial::Polynomial1D< 3 >( 1.0 ,  0.0 , -2.5 ,  1.5 );
		blendingFunctions[2] = Polynomial::Polynomial1D< 3 >( 0.0 ,  0.5 ,  2.0 , -1.5 );
		blendingFunctions[3] = Polynomial::Polynomial1D< 3 >( 0.0 ,  0.0 , -0.5 ,  0.5 );
#endif
	}

	// Adjust so that we are looking for a zero-crossing
	x0 -= isoValue , x1 -= isoValue , x2 -= isoValue , x3 -= isoValue;
	if( !x1 ) return 0;
	if( !x2 ) return 1;

	Polynomial::Polynomial1D< 3 > p = blendingFunctions[0] * x0 + blendingFunctions[1] * x1 + blendingFunctions[2] * x2 + blendingFunctions[3] * x3;
	double roots[3] , _roots[3];
	int rootCount = Polynomial::Roots( p , roots );
	int _rootCount = 0;
	for( int i=0 ; i<rootCount ; i++ ) if( roots[i]>=0 && roots[i]<1 ) _roots[ _rootCount++ ] = roots[i];
	if     ( _rootCount==1 ) return (Real)_roots[0];
	else if( _rootCount==3 ) return (Real)_roots[1];
	else
	{
		std::cout << p << std::endl;
		printf( "Values: %g %g %g %g\n" , x0 , x1 , x2 , x3 );
		printf( "Roots:" ) ; for( int i=0 ; i<rootCount ; i++ ) printf( " %g" , roots[i] ) ; printf( "\n" );
		ERROR_OUT( "Unexpected number of roots: " , _rootCount );
	}
	return 0;
}

template< typename Real >
void IsoSurface2D< Real >::_SetFlags( int resX , ConstPointer( Real ) values , Real isoValue , Pointer( unsigned char ) flags )
{
	for( int i=0 ; i<resX ; i++ ) flags[i] = MarchingSquares::ValueLabel( values[i] , isoValue );
}

template< typename Real >
#ifdef NEW_ISO_SURFACE_2D
void IsoSurface2D< Real >::_SetYVertices( int resX , int y , ConstPointer( Real ) values0 , ConstPointer( Real ) values1 , ConstPointer( Real ) values2 , ConstPointer( Real ) values3 , ConstPointer( unsigned char ) flags1 , ConstPointer( unsigned char ) flags2 , Real isoValue , int interpolationType , Pointer( _Vertex ) isoVertexMap , std::vector< Point2D< Real > >& vertices )
#else // !NEW_ISO_SURFACE_2D
void IsoSurface2D< Real >::_SetYVertices( int resX , int y , ConstPointer( Real ) values0 , ConstPointer( Real ) values1 , ConstPointer( Real ) values2 , ConstPointer( Real ) values3 , ConstPointer( unsigned char ) flags1 , ConstPointer( unsigned char ) flags2 , Real isoValue , int interpolationType , std::unordered_map< long long , int >& isoVertexMap , std::vector< _Vertex >& vertices )
#endif // NEW_ISO_SURFACE_2D
{
	for( int x=0 ; x<resX ; x++ )
	{
#ifdef NEW_ISO_SURFACE_2D
		isoVertexMap[x].index = 0;
#endif // NEW_ISO_SURFACE_2D
		if( flags1[x]!=flags2[x] )
		{
			Real v0 = values0 ? values0[x] : values1[x] , v1 = values1[x] , v2 = values2[x] , v3 = values3 ? values3[x] : values2[x];
			Real iso;
			switch( interpolationType )
			{
			case INTERPOLATE_LINEAR:      iso =     _LinearInterpolant(      v1 , v2 ,      isoValue ) ; break;
			case INTERPOLATE_QUADRATIC:   iso =  _QuadraticInterpolant( v0 , v1 , v2 , v3 , isoValue ) ; break;
			case INTERPOLATE_CUBIC:       iso =      _CubicInterpolant( v0 , v1 , v2 , v3 , isoValue ) ; break;
			case INTERPOLATE_CATMULL_ROM: iso = _CatmullRomInterpolant( v0 , v1 , v2 , v3 , isoValue ) ; break;
			default: ERROR_OUT( "Unrecognized interpolation type: " , interpolationType );
			}
			Point2D< Real > p = Point2D< Real >( (Real)x , (Real)y + iso );
#ifdef NEW_ISO_SURFACE_2D
			isoVertexMap[x].index = 1;
			isoVertexMap[x].p = p;
#else // !NEW_ISO_SURFACE_2D
			long long key = x;
			isoVertexMap[key] = (int)vertices.size();
			vertices.push_back( _Vertex( p , 1 , x , y ) );
#endif // NEW_ISO_SURFACE_2D
		}
	}
#ifdef NEW_ISO_SURFACE_2D
	for( int x=0 ; x<resX ; x++ )
		if( isoVertexMap[x].index ){ isoVertexMap[x].index = (int)vertices.size() ; vertices.push_back( isoVertexMap[x].p ); }
		else isoVertexMap[x].index = -1;
#endif // NEW_ISO_SURFACE_2D
}

template< typename Real >
#ifdef NEW_ISO_SURFACE_2D
void IsoSurface2D< Real >::_SetXVertices( int resX , int y , ConstPointer( Real ) values , ConstPointer( unsigned char ) flags , Real isoValue , int interpolationType , Pointer( _Vertex ) xIsoVertexMap , std::vector< Point2D< Real > >& vertices )
#else // !NEW_ISO_SURFACE_2D
void IsoSurface2D< Real >::_SetXVertices( int resX , int y , ConstPointer( Real ) values , ConstPointer( unsigned char ) flags , Real isoValue , int interpolationType , std::unordered_map< long long , int >& xIsoVertexMap , std::vector< _Vertex >& vertices )
#endif // NEW_ISO_SURFACE_2D
{
	for( int x=0 ; x<resX-1 ; x++ )
	{
#ifdef NEW_ISO_SURFACE_2D
		xIsoVertexMap[x].index = 0;
#endif // NEW_ISO_SURFACE_2D
		int idx1 = x , idx2 = x+1;
		if( flags[idx1]!=flags[idx2] )
		{
			Real v0 = idx1>0 ? values[idx1-1] : values[idx1] , v1 = values[idx1] , v2 = values[idx2] , v3 = idx2<resX-1 ? values[idx2+1] : values[idx2];
			Real iso;
			switch( interpolationType )
			{
			case INTERPOLATE_LINEAR:      iso =     _LinearInterpolant(      v1 , v2 ,      isoValue ) ; break;
			case INTERPOLATE_QUADRATIC:   iso =  _QuadraticInterpolant( v0 , v1 , v2 , v3 , isoValue ) ; break;
			case INTERPOLATE_CUBIC:       iso =      _CubicInterpolant( v0 , v1 , v2 , v3 , isoValue ) ; break;
			case INTERPOLATE_CATMULL_ROM: iso = _CatmullRomInterpolant( v0 , v1 , v2 , v3 , isoValue ) ; break;
			default: ERROR_OUT( "Unrecognized interpolation type: " , interpolationType );
			}
			Point2D< Real > p = Point2D< Real >( (Real)x + iso , (Real)y );
#ifdef NEW_ISO_SURFACE_2D
			xIsoVertexMap[x].index = 1;
			xIsoVertexMap[x].p = p;
#else // !NEW_ISO_SURFACE_2D
			long long key = x;
			xIsoVertexMap[key] = (int)vertices.size();
			vertices.push_back( _Vertex( p , 0 , x , y ) );
#endif // NEW_ISO_SURFACE_2D
		}
	}
#ifdef NEW_ISO_SURFACE_2D
	for( int x=0 ; x<resX-1 ; x++ )
		if( xIsoVertexMap[x].index ){ xIsoVertexMap[x].index = (int)vertices.size() ; vertices.push_back( xIsoVertexMap[x].p ); }
		else xIsoVertexMap[x].index = -1;
#endif // NEW_ISO_SURFACE_2D
}

template< typename Real >
#ifdef NEW_ISO_SURFACE_2D
void IsoSurface2D< Real >::_SetEdges( int resX , int y , ConstPointer( Real ) values1 , ConstPointer( Real ) values2 , Real isoValue , bool fullCaseTable , ConstPointer( _Vertex ) xIsoVertexMap1 , ConstPointer( _Vertex ) xIsoVertexMap2 , ConstPointer( _Vertex ) yIsoVertexMap , std::vector< EdgeIndex >& edges )
#else // !NEW_ISO_SURFACE_2D
void IsoSurface2D< Real >::_SetEdges( int resX , int y , ConstPointer( Real ) values1 , ConstPointer( Real ) values2 , Real isoValue , bool fullCaseTable , const std::unordered_map< long long , int >& xIsoVertexMap1 , const std::unordered_map< long long , int >& xIsoVertexMap2 , const std::unordered_map< long long , int >& yIsoVertexMap , std::vector< EdgeIndex >& edges )
#endif // NEW_ISO_SURFACE_2D
{
	for( int x=0 ; x<resX-1 ; x++ )
	{
		Real _values[Square::CORNERS];
		for( int cx=0 ; cx<2 ; cx++ )
		{
			_values[ Square::CornerIndex(cx,0) ] = values1[ x+cx ];
			_values[ Square::CornerIndex(cx,1) ] = values2[ x+cx ];
		}
		int mcIndex = fullCaseTable ? MarchingSquares::GetFullIndex( _values , isoValue ) : MarchingSquares::GetIndex( _values , isoValue );
		const MarchingSquares::FaceEdges &isoEdges = fullCaseTable ? MarchingSquares::fullCaseTable( mcIndex ) : MarchingSquares::caseTable( mcIndex );

		for( int e=0 ; e<isoEdges.count ; e++ )
		{
			const std::pair< int , int > &isoEdge = isoEdges[e];
			EdgeIndex edge;

			for( int v=0 ; v<2 ; v++ )
			{
				int orientation , i1;
				Square::FactorEdgeIndex( v==0 ? isoEdge.first : isoEdge.second , orientation , i1 );
#ifdef NEW_ISO_SURFACE_2D
				int vIndex;
#else // !NEW_ISO_SURFACE_2D
				long long key;
				std::unordered_map< long long , int >::const_iterator iter;
#endif // NEW_ISO_SURFACE_2D
				bool success;
				switch( orientation )
				{
#ifdef NEW_ISO_SURFACE_2D
				case 0:
					if( i1==0 ){ vIndex = xIsoVertexMap1[x].index ; success = vIndex!=-1; }
					else       { vIndex = xIsoVertexMap2[x].index ; success = vIndex!=-1; }
					break;
				case 1:
					vIndex = yIsoVertexMap[x+i1].index ; success = vIndex!=-1;
					break;
#else // !NEW_ISO_SURFACE_2D
				case 0:
					key = (x   );
					if( i1==0 ){ iter = xIsoVertexMap1.find( key ) ; success = iter!=xIsoVertexMap1.end(); }
					else       { iter = xIsoVertexMap2.find( key ) ; success = iter!=xIsoVertexMap2.end(); }
					break;
				case 1:
					key = (x+i1);
					iter = yIsoVertexMap.find( key ) ; success = iter!=yIsoVertexMap.end();
					break;
#endif // NEW_ISO_SURFACE_2D
				}

				if( !success ) ERROR_OUT( "Could not find iso-vertex in map" );
#ifdef NEW_ISO_SURFACE_2D
				edge[v] = vIndex;
#else // !NEW_ISO_SURFACE_2D
				edge[v] = iter->second;
#endif // NEW_ISO_SURFACE_2D
			}
			edges.push_back( edge );
		}
	}
}
