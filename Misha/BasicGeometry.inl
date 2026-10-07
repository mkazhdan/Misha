/*
Copyright (c) 2026, Michael Kazhdan and Hongyi Liu
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

template< HasXForm XForm >
OrientedEmbeddedMesh Implicit( const MishaK::EquationParser::Node & F , unsigned int res , XForm && CubeToWorld , double isoValue )
{
	auto GridToWorld = [&]( MishaK::Point< double , 3 > gCoord ){ return CubeToWorld( gCoord / static_cast< double >( res-1 ) ); };

	MishaK::EquationParser::Node Fx = F.d( 0 ) , Fy = F.d( 1 ) , Fz = F.d( 2 );

	MishaK::RegularGrid< 3 , double > grid;
	grid.resize( res );
	for( unsigned int i=0 ; i<res ; i++ ) for( unsigned int j=0 ; j<res ; j++ ) for( unsigned int k=0 ; k<res ; k++ )
	{
		MishaK::Point< double , 3 > p = GridToWorld( MishaK::Point< double , 3 >( static_cast< double >( i ) , static_cast< double >( j ) , static_cast< double >( k ) ) );
		grid(i,j,k) = F( &p[0] );
	}

	OrientedEmbeddedMesh mesh;
	std::vector< MishaK::Point< double , Dim > > & vertices  = std::get< 0 >( mesh );
	std::vector< MishaK::Point< double , Dim > > & normals   = std::get< 1 >( mesh );
	std::vector< MishaK::SimplexIndex< K > >     & simplices = std::get< 2 >( mesh );

	MishaK::IsoSurface3D< double , unsigned int >::Extract( grid , isoValue , vertices , simplices , true , MishaK::IsoSurface3D< double , unsigned int >::INTERPOLATE_LINEAR , true );
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i] = GridToWorld( vertices[i] );

	normals.resize( vertices.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ )
	{
		normals[i] = MishaK::Point< double , 3 >( Fx( &vertices[i][0] ) ,  Fx( &vertices[i][1] ) ,  Fz( &vertices[i][2] ) );
		normals[i] /= MishaK::Point< double , 3 >::Length( normals[i] );
	}

	return mesh;
}

inline OrientedEmbeddedMesh Implicit( const MishaK::EquationParser::Node & F , unsigned int res , double isoValue )
{
	return Implicit( F , res , []( MishaK::Point< double , Dim > p ){ return p; } , isoValue );
}

template< HasXForm XForm >
OrientedEmbeddedMesh Implicit( std::string eqn_xyz , unsigned int res , XForm && CubeToWorld , double isoValue )
{
	return Implicit( MishaK::EquationParser::Node( eqn_xyz , { "x" , "y" , "z" } ) , res , std::forward< XForm >( CubeToWorld ) , isoValue );
}

inline OrientedEmbeddedMesh ImplicitSphere( unsigned int res , double radius , double scale )
{
	std::string eqn = "x^2 + y^2 + z^2";
	auto CubeToWorld = [&]( MishaK::Point< double , 3 > p )
		{
			p = ( p - MishaK::Point< double , 3 >(0.5,0.5,0.5) ) * 2;
			return p / scale;
		};

	return Implicit( eqn , res , CubeToWorld , 1. );
}

inline OrientedEmbeddedMesh ImplicitTorus( unsigned int res , double radius1 , double radius2 , double scale )
{
	// The square distance from a circle of radius r is:
	//		d(p) = < p - p/|p|*r , p - p/|p|*r >
	//		     = < p , p > + r^2 - 2 < p , p > * r / | p |
	std::string eqn = "y^2 + x^2 + z^2 + " + std::to_string( radius1 ) + "^2 - 2 * ( x^2 + z^2 ) * " + std::to_string( radius1 ) + " / ( x^2 + z^2 )^0.5";

	auto CubeToWorld = [&]( MishaK::Point< double , 3 > p )
		{
			p = ( p - MishaK::Point< double , 3 >(0.5,0.5,0.5) ) * 2 * ( radius1+radius2 );
			return p / scale;
		};

	return Implicit( eqn , res , CubeToWorld , radius2*radius2 );
}

inline OrientedEmbeddedMesh EmbedTorus( const Mesh & _mesh , double radius1 , double radius2 )
{
	OrientedEmbeddedMesh mesh;

	std::get< 2 >( mesh ) = std::get< 1 >( _mesh );

	std::vector< MishaK::Point< double , Dim > > & vertices = std::get< 0 >( mesh );
	std::vector< MishaK::Point< double , Dim > > & normals  = std::get< 1 >( mesh );
	const std::vector< MishaK::Point< double , K > > & _vertices = std::get< 0 >( _mesh );

	Torus torus( radius1 , radius2 );
	vertices.resize( _vertices.size() ) , normals.resize( _vertices.size() );
	for( unsigned int i=0 ; i<_vertices.size() ; i++ )
	{
		vertices[i] = torus( _vertices[i] );
		normals[i] = torus.n( _vertices[i] );
	}

	return mesh;
}

#ifdef USE_Q_HULL
inline OrientedEmbeddedMesh RandomEmbeddedSphere( unsigned int numVertices , double radius , MishaK::Point< double , Dim > scale )
{
	OrientedEmbeddedMesh mesh;

	std::vector< MishaK::Point< double , Dim > > & vertices  = std::get< 0 >( mesh );
	std::vector< MishaK::Point< double , Dim > > & normals   = std::get< 1 >( mesh );
	std::vector< MishaK::SimplexIndex< K > >     & simplices = std::get< 2 >( mesh );

	vertices.resize( numVertices );
	normals.resize( numVertices );
	auto RandomSpherePoint = [&]( void )
		{
			MishaK::Point< double , Dim > p = MishaK::RandomSpherePoint< double , Dim >();
			for( unsigned int d=0 ; d<Dim ; d++ ) p[d] *= scale[d];
			return p;
		};
	for( unsigned int i=0 ; i<numVertices ; i++ ) vertices[i] = RandomSpherePoint();

	std::vector< std::vector< int > > polygons;
	qconvex_misha( &vertices[0][0] , static_cast< unsigned int >( vertices.size() ) , Dim , polygons );
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i] = normals[i] = vertices[i] /= MishaK::Point< double , Dim >::Length( vertices[i] );

	for( unsigned int i=0 ; i<polygons.size() ; i++ )
	{
		std::vector< int > & p = polygons[i];
		size_t sz = p.size();

		MishaK::Point3D< double > c , n;
		for( unsigned int j=0 ; j<sz ; j++ ) c += vertices[ p[j] ];
		c /= static_cast< double >( sz );
		for( unsigned int j=0 ; j<sz ; j++ ) n += MishaK::Point3D< double >::CrossProduct( vertices[ p[(j+0)%sz] ] - c , vertices[ p[(j+1)%sz] ] - c );
		if( MishaK::Point3D< double >::Dot( c , n )<0 ) for( unsigned int j=0 ; j<sz/2 ; j++ ) std::swap( p[j] , p[sz-1-j] );
	}

	for( unsigned int i=0 ; i<polygons.size() ; i++ )
	{
		const std::vector< int > &p = polygons[i];
		if( p.size()<3 ) MK_ERROR_OUT( "Degenerate face" );

		MishaK::SimplexIndex< K > tri;
		tri[0] = p[0];

		for( unsigned int j=1 ; j<p.size()-1 ; j++ )
		{
			tri[1] = p[j+0];
			tri[2] = p[j+1];
			simplices.push_back( tri );
		}
	}

	for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i] *= radius;

	return mesh;
}

inline Mesh RandomTorusWithParameterization( std::vector< MishaK::Simplex< double , 2 , 2 > > & pSimplices , unsigned int numVertices , MishaK::Point< double , K > scale )
{
	Mesh mesh;

	std::vector< MishaK::Point< double , K > > & vertices  = std::get< 0 >( mesh );
	std::vector< MishaK::SimplexIndex< K > >   & simplices = std::get< 1 >( mesh );

	vertices.resize( numVertices );
	std::vector< MishaK::Point< double , 2 > > _pVertices;
	auto RandomSquarePoint = [&]( void )
		{
			MishaK::Point< double , K > p;
			while( true )
			{
				bool inSquare = true;
				for( unsigned int k=0 ; k<K ; k++ )
				{
					p[k] = MishaK::Random< double >() * scale[k];
					if( p[k]<0 || p[k]>1 ) inSquare = false;
				}
				if( inSquare ) return p;
			}
		};
	for( unsigned int i=0 ; i<numVertices ; i++ ) vertices[i] = MishaK::Point2D< double >( MishaK::Random< double >() * scale[0] , MishaK::Random< double >() * scale[1] );

	for( unsigned int i=0 ; i<numVertices ; i++ ) for( int j=-1 ; j<=1 ; j++ ) for( int k=-1 ; k<=1 ; k++ ) _pVertices.push_back( MishaK::Point< double , 2 >( vertices[i][0]+j*scale[0] , vertices[i][1]+k*scale[1] ) );
	std::vector< std::vector< int > > polygons;
	qdelaunay_misha( &_pVertices[0][0] , static_cast< unsigned int >( _pVertices.size() ) , 2 , polygons );

	for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i][0] *= 2.*M_PI / scale[0] , vertices[i][1] *= 2.*M_PI / scale[1];
	for( unsigned int i=0 ; i<_pVertices.size() ; i++ ) _pVertices[i][0] *= 2.*M_PI / scale[0] , _pVertices[i][1] *= 2.*M_PI / scale[1];

	for( unsigned int i=0 ; i<polygons.size() ; i++ )
	{
		std::vector< int > & p = polygons[i];
		size_t sz = p.size();

		MishaK::Point2D< double > c;
		MishaK::Point3D< double > n;
		for( unsigned int j=0 ; j<sz ; j++ ) c += _pVertices[ p[j] ];
		c /= static_cast< double >( sz );
		for( unsigned int j=0 ; j<sz ; j++ )
		{
			MishaK::Point3D< double > v1( _pVertices[ p[(j+0)%sz] ][0] - c[0] , _pVertices[ p[(j+0)%sz] ][1] - c[1] , 0 );
			MishaK::Point3D< double > v2( _pVertices[ p[(j+1)%sz] ][0] - c[0] , _pVertices[ p[(j+1)%sz] ][1] - c[1] , 0 );
			n += MishaK::Point3D< double >::CrossProduct( v1 , v2 );
		}
		if( n[2]<0 ) for( unsigned int j=0 ; j<sz/2 ; j++ ) std::swap( p[j] , p[sz-1-j] );
	}

	for( unsigned int i=0 ; i<polygons.size() ; i++ )
	{
		const std::vector< int > &p = polygons[i];
		if( p.size()<3 ) MK_ERROR_OUT( "Degenerate face: " , p.size() );


		// Will only add a polygon if:
		// 1. It contains some vertex in [0,1]x[0,1]
		// 2. All its vertices are in [-1,1]x[-1,1]
		// [ASSUMING]
		// All the vertices are in [-1,1]x[-1,1], [0,2]x[-1,1], [-1,1]x[0,2], or [0,2][0,2].
		// Note:
		// 1. If it has a vertex not in [-1,1]x[-1,1], there is a left/down-shift of the triangle will have all three vertices in [-1,1]x[-1,1]
		// 2. If it has a vertex in [0,1]x[0,1] and all vertices are in [-1,1]x[-1,1], no left/down-shift will have a vertex in [0,1]x[0,1]

		bool addPolygon;
		{
			int xMin=1 , xMax=-1 , yMin=1 , yMax=-1;
			for( unsigned int j=0 ; j<p.size() ; j++ )
			{
				MishaK::Point< double , 2 > v = _pVertices[ p[j] ];
				int idx[2];
				if     ( v[0]< 0.      ) idx[0] = -1;
				else if( v[0]>=2.*M_PI ) idx[0] =  1;
				else                     idx[0] =  0;

				if     ( v[1]< 0.      ) idx[1] = -1;
				else if( v[1]>=2.*M_PI ) idx[1] =  1;
				else                     idx[1] =  0;

				xMin = std::min< int >( xMin , idx[0] );
				xMax = std::max< int >( xMax , idx[0] );
				yMin = std::min< int >( yMin , idx[1] );
				yMax = std::max< int >( yMax , idx[1] );
			}
			addPolygon = xMax==0 && yMax==0;
			if( addPolygon )
			{
				for( unsigned int j=0 ; j<p.size() ; j++ ) for( unsigned int k=0 ; k<j ; k++ ) if( p[j]/9 == p[k]/9 )
				{
					MK_WARN( "Degenerate polygon edge" );
					addPolygon = false;
				}
			}
		}
		if( addPolygon )
		{
			MishaK::SimplexIndex< 2 > tri;
			MishaK::Simplex< double , 2 , 2 > pSimplex;
			tri[0] = p[0]/9;
			pSimplex[0] = _pVertices[ p[0] ];

			for( unsigned int j=1 ; j<p.size()-1 ; j++ )
			{
				tri[1] = p[j+0]/9 , tri[2] = p[j+1]/9;
				pSimplex[1] = _pVertices[ p[j+0] ] , pSimplex[2] = _pVertices[ p[j+1] ];
				simplices.push_back( tri );
				pSimplices.push_back( pSimplex );
			}
		}
	}

	return mesh;
}

#endif // USE_Q_HULL

inline Mesh RegularTorusWithParameterization( std::vector< MishaK::Simplex< double , 2 , 2 > > & pSimplices , unsigned int res )
{
	Mesh mesh;

	std::vector< MishaK::Point< double , K > > & vertices  = std::get< 0 >( mesh );
	std::vector< MishaK::SimplexIndex< K > >   & simplices = std::get< 1 >( mesh );

	vertices.resize( res*res );
	simplices.resize( res*res*2 );
	pSimplices.resize( res*res*2 );

	// Set the vertices
	for( unsigned int j=0 ; j<res ; j++ ) for( unsigned int i=0 ; i<res ; i++ ) vertices[ j*res+i ] = MishaK::Point< double , K >( static_cast< double >(i)/res , static_cast< double >(j)/res ) * 2. * M_PI;

	// Set the (parameterization) triangles
	for( unsigned int j=0 ; j<res ; j++ ) for( unsigned int i=0 ; i<res ; i++ )
	{
		unsigned int idx = j*res + i;
		double x = static_cast< double >( i ) , y = static_cast< double >(j);

		simplices[2*idx+0][0] = ((j+0)%res)*res + (i+0)%res;
		simplices[2*idx+0][1] = ((j+0)%res)*res + (i+1)%res;
		simplices[2*idx+0][2] = ((j+1)%res)*res + (i+0)%res;

		simplices[2*idx+1][0] = ((j+1)%res)*res + (i+1)%res;
		simplices[2*idx+1][1] = ((j+1)%res)*res + (i+0)%res;
		simplices[2*idx+1][2] = ((j+0)%res)*res + (i+1)%res;

		pSimplices[2*idx+0][0] = MishaK::Point< double , 2 >( x+0. , y+0. ) / res * 2. * M_PI;
		pSimplices[2*idx+0][1] = MishaK::Point< double , 2 >( x+1. , y+0. ) / res * 2. * M_PI;
		pSimplices[2*idx+0][2] = MishaK::Point< double , 2 >( x+0. , y+1. ) / res * 2. * M_PI;

		pSimplices[2*idx+1][0] = MishaK::Point< double , 2 >( x+1. , y+1. ) / res * 2. * M_PI;
		pSimplices[2*idx+1][1] = MishaK::Point< double , 2 >( x+0. , y+1. ) / res * 2. * M_PI;
		pSimplices[2*idx+1][2] = MishaK::Point< double , 2 >( x+1. , y+0. ) / res * 2. * M_PI;
	}

	return mesh;

}

inline Torus::Torus( double r1 , double r2 ) : _r1(r1) , _r2(r2)
{
	if( _r1<r2 ) MK_ERROR_OUT( "First radius needs to exceed second radius" );
}

inline MishaK::Point< double , 3 > Torus::operator()( double s , double t ) const
{
	static const MishaK::Point< double , 3 > X(1,0,0) , Y(0,1,0) , Z(0,0,1);
	MishaK::Point< double , 3 > d = X * sin(s) + Z * cos(s);
	return d * _r1  + ( d * cos(t) + Y * sin(t) ) * _r2;
}

inline MishaK::Point< double , 2 > Torus::operator()( MishaK::Point< double , 3 > p ) const
{
	static const MishaK::Point< double , 3 > X(1,0,0) , Y(0,1,0) , Z(0,0,1);

	double s = atan2( p[0] , p[2] );
	MishaK::Point< double , 3 > d = X * sin(s) + Z * cos(s);
	
	double t = atan2( MishaK::Point< double , 3 >::Dot( p - d*_r1 , Y ) , MishaK::Point< double , 3 >::Dot( p - d*_r1 , d ) );

	return MishaK::Point< double , 2 >( s , t );
}



inline MishaK::Point< double , 3 > Torus::ds( double s , double t ) const
{
	static const MishaK::Point< double , 3 > X(1,0,0) , Y(0,1,0) , Z(0,0,1);
	MishaK::Point< double , 3 > _d = ( cos(s) * X - sin(s) * Z );
	return _d * ( _r1 + cos(t) * _r2 );
}

inline MishaK::Point< double , 3 > Torus::dt( double s , double t ) const
{
	static const MishaK::Point< double , 3 > X(1,0,0) , Y(0,1,0) , Z(0,0,1);
	MishaK::Point< double , 3 > d = X * sin(s) + Z * cos(s);
	return ( -d * sin(t) + Y * cos(t) );
}

inline MishaK::Point< double , 3 > Torus::n( double s , double t ) const
{
	MishaK::Point< double , 3 > n = MishaK::Point3D< double >::CrossProduct( ds( s , t ) , dt( s , t ) );
	return n /= MishaK::Point3D< double >::Length( n );
}

inline MishaK::SquareMatrix< double , 2 > Torus::g( double s , double t ) const
{
	MishaK::Point< double , 3 > _ds = ds(s,t) , _dt=dt(s,t);
	MishaK::SquareMatrix< double , 2 > _g;
	_g(0,0) = MishaK::Point< double , 3 >::Dot( _ds , _ds );
	_g(0,1) = MishaK::Point< double , 3 >::Dot( _ds , _dt );
	_g(1,0) = MishaK::Point< double , 3 >::Dot( _dt , _ds );
	_g(1,1) = MishaK::Point< double , 3 >::Dot( _dt , _dt );
	return _g;
}

inline MishaK::Point< double , 3 > Torus::operator()( MishaK::Point< double , 2 > st ) const { return operator()( st[0] , st[1] ); }
inline MishaK::Point< double , 3 > Torus::ds( MishaK::Point< double , 2 > st ) const { return ds( st[0] , st[1] ); }
inline MishaK::Point< double , 3 > Torus::dt( MishaK::Point< double , 2 > st ) const { return dt( st[0] , st[1] ); }
inline MishaK::Point< double , 3 > Torus::n( MishaK::Point< double , 2 > st ) const { return n( st[0] , st[1] ); }
inline MishaK::SquareMatrix< double , 2 > Torus::g( MishaK::Point< double , 2 > st ) const { return g( st[0] , st[1] ); }



inline OrientedEmbeddedMesh Tetrahedron( void )
{
	OrientedEmbeddedMesh mesh;
	std::vector< MishaK::Point< double , Dim > > & vertices = std::get< 0 >( mesh );
	std::vector< MishaK::SimplexIndex< K > > & simplices = std::get< 2 >( mesh );

	vertices.resize( 4 );
	simplices.resize( 4 );

	// v[] = { (0,a,0) , (cos(0),0,sin(0)) , (cos(2Pi/3),0,sin(2Pi/3)) , (cos(4Pi/3),0,sin(4Pi/3)) }
	// || v[0] - v[1] ||^2 = || v[1] - v[2] ||^2
	// => 1 + a^2 = [ 1 - cos(2Pi/3) ]^2 + [ sin(2Pi/3) ]^2
	//            = 1 - 2 * cos( 2PI/3 ) + 1
	//            = 2 - 2 * cos( 2PI/3 ) );
	// => a = sqrt( 1 - 2 * cos( 2Pi/3 ) )
	double theta = 2 * M_PI / 3.;
	double a = sqrt( 1. - 2. * cos(theta) );

	vertices[0] = MishaK::Point< double , Dim > ( 0 , a , 0 );
	vertices[1] = MishaK::Point< double , Dim > ( cos(0*theta) , 0 , sin(0*theta) );
	vertices[2] = MishaK::Point< double , Dim > ( cos(1*theta) , 0 , sin(1*theta) );
	vertices[3] = MishaK::Point< double , Dim > ( cos(2*theta) , 0 , sin(2*theta) );

	{
		MishaK::Point< double , Dim > c = ( vertices[0] + vertices[1] + vertices[2] + vertices[3] ) / 4.;
		for( unsigned int i=0 ; i<4 ; i++ )
		{
			vertices[i] -= c;
			vertices[i] /= MishaK::Point< double , Dim >::Length( vertices[i] );
		}
	}

	for( unsigned int i=0 ; i<3 ; i++ ) 
	{
		simplices[i][0] = 0;
		simplices[i][1] = (i+1)%3+1;
		simplices[i][2] = (i+0)%3+1;
	}
	simplices[3][0] = 1;
	simplices[3][1] = 2;
	simplices[3][2] = 3;

	return mesh;
}

inline OrientedEmbeddedMesh Octahedron( void )
{
	OrientedEmbeddedMesh mesh;
	std::vector< MishaK::Point< double , Dim > > & vertices = std::get< 0 >( mesh );
	std::vector< MishaK::SimplexIndex< K > > & simplices = std::get< 2 >( mesh );

	simplices.resize( 8 );
	vertices.resize( 6 );
	vertices[0] = MishaK::Point< double , Dim >(  0 ,  0 , -1 );
	vertices[1] = MishaK::Point< double , Dim >(  1 ,  0 ,  0 );
	vertices[2] = MishaK::Point< double , Dim >(  0 ,  1 ,  0 );
	vertices[3] = MishaK::Point< double , Dim >( -1 ,  0 ,  0 );
	vertices[4] = MishaK::Point< double , Dim >(  0 , -1 ,  0 );
	vertices[5] = MishaK::Point< double , Dim >(  0 ,  0 ,  1 );
	for( unsigned int i=0 ; i<4 ; i++ )
	{
		simplices[i][0] = 0;
		simplices[i][1] = (i+2)%4 + 1;
		simplices[i][2] = (i+1)%4 + 1;
		simplices[4+i][0] = 5;
		simplices[4+i][1] = (i+1)%4 + 1;
		simplices[4+i][2] = (i+2)%4 + 1;
	}

	return mesh;
}

inline OrientedEmbeddedMesh Icosahedron( void )
{
	OrientedEmbeddedMesh mesh;
	std::vector< MishaK::Point< double , Dim > > & vertices = std::get< 0 >( mesh );
	std::vector< MishaK::SimplexIndex< K > > & simplices = std::get< 2 >( mesh );

	simplices.resize( 20 );
	vertices.resize( 12 );
	vertices[0] = MishaK::Point< double , Dim >( 0 , -1 , 0 );
	vertices[11] = MishaK::Point< double , Dim >( 0 , 1 , 0 );

	double theta = 2.*M_PI / 5;
	double phi = M_PI / 3.;
	for( unsigned int i=0 ; i<5 ; i++ )
	{
		vertices[1+i] = MishaK::Point< double , Dim >( cos(i*theta)*sin(2*phi) , cos(2*phi) , sin(i*theta)*sin(2*phi) );
		vertices[6+i] = MishaK::Point< double , Dim >( cos(i*theta+theta/2.)*sin(phi) , cos(phi) , sin(i*theta+theta/2.)*sin(phi) );
	}
	for( unsigned int i=0 ; i<5 ; i++ )
	{
		simplices[0+i][0] = 0;
		simplices[0+i][1] = (i+0)%5+1;
		simplices[0+i][2] = (i+1)%5+1;

		simplices[5+i][0] = (i+0)%5+1;
		simplices[5+i][1] = (i+4)%5+6;
		simplices[5+i][2] = (i+0)%5+6;

		simplices[10+i][0] = (i+1)%5+1;
		simplices[10+i][1] = (i+0)%5+1;
		simplices[10+i][2] = (i+0)%5+6;

		simplices[15+i][0] = 11;
		simplices[15+i][1] = (i+1)%5+6;
		simplices[15+i][2] = (i+0)%5+6;
	}

	return mesh;
}
