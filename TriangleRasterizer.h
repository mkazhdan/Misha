#ifndef TRIANGLE_RASTERIZER
#define TRIANGLE_RASTERIZER

#include <Misha/Algebra.h>

namespace Rasterizer
{
	template< typename Real , typename Data >
	struct ClippedTriangle
	{
	protected:
		static const unsigned int _MAX_VERTS = 7;
		Point2D< Real > _v[ _MAX_VERTS ];
		Data _d[ _MAX_VERTS ];
		unsigned int _sz;
	public:
		Point2D< Real > &operator[] ( unsigned int idx ){ return _v[idx]; }
		const Point2D< Real > &operator[] ( unsigned int idx ) const { return _v[idx]; }
		Data &data( unsigned int idx ){ return _d[idx]; }
		const Data &data( unsigned int idx ) const { return _d[idx]; }
		unsigned int size( void ) const { return _sz; }

		ClippedTriangle( void ) : _sz(0){}
		ClippedTriangle( Point2D< Real > v0 , Point2D< Real > v1 , Point2D< Real > v2 , Data d0=0 , Data d1=0 , Data d2=0 ) :_sz(3)
		{
			_v[0] = v0 , _v[1] = v1 , _v[2] = v2;
			_d[0] = d0 , _d[1] = d1 , _d[2] = d2;
		}

		void split( Point2D< Real > normal , Real offset , ClippedTriangle &back , ClippedTriangle &front ) const;
	};

	template< typename Real , typename Data >
	struct TriangleSample : public VectorSpace< Real , TriangleSample< Real , Data > >
	{
	protected:
		Point2D< Real > _v;
		Data _d;
	public:
		Point2D< Real > &operator()( void ){ return _v; }
		const Point2D< Real > &operator()( void ) const { return _v; }
		Data &data( void ){ return _d; }
		const Data &data( void ) const { return _d; }
		TriangleSample( void ){}
		TriangleSample( Point2D< Real > v , Data d=0 ) : _v(v) , _d(d){}

		void Add( const TriangleSample< Real , Data > &e ){ _v += e._v , _d += e._d; }
		void Scale( Real s ){ _v *= s , _d *= s; }
	};

	template< typename Real , typename Data >
	std::ostream &operator << ( std::ostream &os , const ClippedTriangle< Real , Data > &tri )
	{
		for( size_t i=0 ; i<tri.size() ; i++ ) os << tri[i] << "\t" << tri.data(i) << std::endl;
		return os;
	}


	template< typename Real , typename Data >
	std::vector< ClippedTriangle< Real , Data > > RasterizeZ2( Point2D< Real > v0 , Point2D< Real > v1 , Point2D< Real > v2 , Data d0 , Data d1 , Data d2 );

	template< typename Real , typename Data >
#if 1
	std::vector< TriangleSample< int , Data > > SampleZ2( Point2D< Real > v0 , Point2D< Real > v1 , Point2D< Real > v2 , Data d0 , Data d1 , Data d2 , bool interiorOnly=true );
#else
	std::vector< TriangleSample< int , Data > > SampleZ2( Point2D< Real > v0 , Point2D< Real > v1 , Point2D< Real > v2 , Data d0 , Data d1 , Data d2 );
#endif


	////////////////////
	// Implementation //
	////////////////////

	template< class Real , typename Data >
	void ClippedTriangle< Real , Data >::split( Point2D< Real > normal , Real offset , ClippedTriangle &back , ClippedTriangle &front ) const
	{
		Real values[ _MAX_VERTS ];
		bool frontSet = false , backSet = false;

		// Evaluate the plane's function at the vertices and mark if front/back vertices have been found
		for( unsigned int i=0 ; i<size() ; i++ )
		{
			values[i] = Point2D< Real >::Dot( _v[i] , normal ) - offset;
			backSet |= ( values[i]<0 ) , frontSet |= ( values[i]>0 );
		}

		// If no vertices are in front or no vertices are behind, life is easy
		if( !frontSet ){  back = *this ; front = ClippedTriangle() ; return; }
		if( !backSet  ){ front = *this ;  back = ClippedTriangle() ; return; }

		enum{ STATE_BACK , STATE_FRONT , STATE_BOUNDARY };
		auto State = []( Real value )
		{
			if     ( value<0 ) return STATE_BACK;
			else if( value>0 ) return STATE_FRONT;
			else               return STATE_BOUNDARY;
		};
		struct _Vertex
		{
			Point2D< Real > position;
			Data data;
			int state;
			_Vertex( Point2D< Real > p=Point2D< Real >() , Data d=0 , int s=STATE_BOUNDARY ) : position(p) , data(d) , state(s) {}
		};

		_Vertex _vertices[ 2*_MAX_VERTS ];
		unsigned int _vCount=0;
		for( unsigned int i0=0 ; i0<size() ; i0++ )
		{
			int i1 = (i0+1)%size();
			_vertices[_vCount++] = _Vertex( _v[i0] , _d[i0] , State( values[i0] ) );
			// If the (clipped) triangle edge crosses the line, we need to create a new vertex
			if( values[i0]*values[i1]<0 )
			{
				Real t = values[i0] / ( values[i0] - values[i1] );
				_vertices[_vCount++] = _Vertex( _v[i1]*t + _v[i0]*(Real)(1.-t) , _d[i1]*t + _d[i0]*(Real)(1.-t) , STATE_BOUNDARY );
			}
		}

		back = front = ClippedTriangle();
		for( int i=0 ; i<_vCount ; i++ )
		{
			if( _vertices[i].state!=STATE_BACK  ){ front._v[front._sz] = _vertices[i].position , front._d[front._sz] = _vertices[i].data ; front._sz++; }
			if( _vertices[i].state!=STATE_FRONT ){  back._v[ back._sz] = _vertices[i].position ,  back._d[ back._sz] = _vertices[i].data ;  back._sz++; }
		}
		if( front.size()>_MAX_VERTS || back.size()>_MAX_VERTS ) ERROR_OUT( "Uh oh: " , front.size() , " " , back.size() , " " , _MAX_VERTS );
	}

	template< typename Real , typename Data >
	std::vector< ClippedTriangle< Real , Data > > RasterizeZ2( Point2D< Real > v0 , Point2D< Real > v1 , Point2D< Real > v2 , Data d0 , Data d1 , Data d2 )
	{
		std::vector< ClippedTriangle< Real , Data > > fragments;
		ClippedTriangle< Real , Data > t( v0 , v1 , v2 , d0 , d1 , d2 );
		// Get the (integer) bounding box
		int xMin , xMax , yMin , yMax;
		{
			Real _xMin = std::min< Real >( t[0][0] , std::min< Real >( t[1][0] , t[2][0] ) );
			Real _xMax = std::max< Real >( t[0][0] , std::max< Real >( t[1][0] , t[2][0] ) );
			Real _yMin = std::min< Real >( t[0][1] , std::min< Real >( t[1][1] , t[2][1] ) );
			Real _yMax = std::max< Real >( t[0][1] , std::max< Real >( t[1][1] , t[2][1] ) );
			xMin = (int)floor( _xMin ) , xMax = (int)ceil( _xMax ) , yMin = (int)floor( _yMin ) , yMax = (int)ceil( _yMax );
		}
		std::vector< ClippedTriangle< Real , Data > > strips;
		{
			ClippedTriangle< Real , Data > back , front;

			front = t;
			Point2D< Real > normal = Point2D< Real >( 1 , 0 );
			for( int x=xMin+1 ; x<=xMax-1 ; x++ )
			{
				front.split( normal , (Real)x , back , front );
				if( back.size() ) strips.push_back( back );
			}
			if( front.size() ) strips.push_back( front );
		}

		for( int i=0 ; i<strips.size() ; i++ )
		{
			ClippedTriangle< Real , Data > back , front;

			front = strips[i];
			Point2D< Real > normal = Point2D< Real >( 0 , 1 );
			for( int y=yMin+1 ; y<=yMax-1 ; y++ )
			{
				front.split( normal , (Real)y , back , front );
				if( back.size() ) fragments.push_back( back );
			}
			if( front.size() ) fragments.push_back( front );
		}

		return fragments;
	}

	// Computes the sub-set of integeral lattice points interior to the triangle
	template< typename Real , typename Data >
	std::vector< TriangleSample< int , Data > > SampleZ2( Point2D< Real > v0 , Point2D< Real > v1 , Point2D< Real > v2 , Data d0 , Data d1 , Data d2 )
	{
		auto AddYSamples = []( int x , TriangleSample< Real , Data > s0 , TriangleSample< Real , Data > s1 , std::vector< TriangleSample< int , Data > > &samples )
			{
				// Make sure the first point is below the second
				if( s0()[1]>s1()[1] ) std::swap( s0 , s1 );
				int startY = (int)ceil( s0()[1] ) , endY = (int)floor( s1()[1] );

				for( int y=startY ; y<=endY ; y++ )
				{
					Real s = ( (Real)y - s0()[1] ) / ( s1()[1] - s0()[1] );
					samples.push_back( TriangleSample< int , Data >( Point2D< int >(x,y) , s0.data()*(Real)(1.-s) + s1.data()*s ) );
				}
			};

		std::vector< TriangleSample< int , Data > > samples;

		TriangleSample< Real , Data > v[] = { TriangleSample< Real , Data >( v0 , d0 ) , TriangleSample< Real , Data >( v1 , d1 ) , TriangleSample< Real , Data >( v2 , d2 ) };

		// Sort the triangle corners by x-value
		std::sort( v , v+3 , []( TriangleSample< Real , Data > v1 , TriangleSample< Real , Data > v2 ){ return v1()[0] < v2()[0]; } );

		TriangleSample< Real , Data > mid[2];

		// [WARNING] If the x-coordinate of v[1] is an integer, the samples on the vertical line through the x-coordinate of v[1] will be
		// generated twice.

		// If the first triangle is not degenerate
		if( v[0]()[0]!=v[1]()[0] )
		{
			Real x0 = v[0]()[0] , x1 = v[1]()[0];
			TriangleSample< Real , Data > d[2];
			// Get the rates of change of the y-component and the data values
			d[0] = ( v[1] - v[0] ) / ( v[1]()[0] - v[0]()[0] );
			d[1] = ( v[2] - v[0] ) / ( v[2]()[0] - v[0]()[0] );

			// Get the two end-points
			mid[0] = v[0] + d[0] * ( x1-x0 );
			mid[1] = v[0] + d[1] * ( x1-x0 );

			int startX = (int)ceil( x0 ) , endX = (int)floor( x1 );
			for( int x=startX ; x<=endX ; x++ )
			{
				Real s = ( (Real)x - x0 ) / ( x1 - x0 );
				AddYSamples( x , v[0] * ( (Real)1.-s ) + mid[0] * s , v[0] * ( (Real)1.-s ) + mid[1] * s , samples );
			}
		}
		else mid[0] = v[0] , mid[1] = v[1];

		// If the second triangle is not degenerate
		if( v[1]()[0]!=v[2]()[0] )
		{
			Real x0 = v[1]()[0] , x1 = v[2]()[0];
			int startX = (int)ceil( x0 ) , endX = (int)floor( x1 );

			for( int x=startX ; x<=endX ; x++ )
			{
				Real s = ( (Real)x - x0 ) / ( x1 - x0 );
				AddYSamples( x , mid[0] * ( (Real)1.-s ) + v[2] * s , mid[1] * ( (Real)1.-s ) + v[2] * s , samples );
			}
		}
		return samples;
	}
};

#endif // TRIANGLE_RASTERIZER
