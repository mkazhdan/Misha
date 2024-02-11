#include "KDtree.h"
#include "Geometry.h"

template< class Real >
struct KDMesh
{
protected:
	const std::vector< Point3D< Real > >& _v;
	const std::vector< TriangleIndex >& _f;
	KDtree< Real , 3 > _kd;
	mutable std::vector< std::vector< int > > _af;
public:
	KDMesh( const std::vector< Point3D< Real > >& v , const std::vector< TriangleIndex >& f ) : _v(v) , _f(f) , _kd( &v[0][0] , v.size() , sizeof( Point3D< Real > ) ){  }

	void need_adjacent_faces() const
	{
		if( _af.size() == _v.size() ) return;

		_af.clear();
		_af.resize( _v.size() );
		for( size_t i = 0 ; i<_f.size() ; i++ )
			for( int j=0 ; j<3 ; j++ ) 
				_af[ _f[i][j] ].push_back((int) i);
	}
	int closest_vertex_index( Point3D< Real > p ) const
	{
		Real maxDist2 = -1;
		Point3D< Real > *nearest = ( Point3D< Real >* ) _kd.closest_to_pt( &p[0] , maxDist2 );
		if( nearest ) return (int)( ( size_t( nearest ) - size_t( &_v[0][0] ) ) / sizeof( Point3D< Real > ) );
		else return -1;
	}

	// Returns the closest face on the mesh when the nearest vertex is known

	// [IN]  vert: The index of the nearest point to pIn
	// [IN]  pIn:  The point of interest
	// [OUT] pOut: The point on the triangle that is closest
	// [OUT] a:    The barycentric coordinates of the point on the face
	// [Return]    The index of the face
	template< bool TestAdjacentFaces >
	int closest_point_on_face( int vert , Point3D< Real > pIn , Point3D< Real > &pOut , Real a[3] )
	{
		if( vert<0 ) return -1;
		int fIndex = -1;
		Real dst = 0;
		if( TestAdjacentFaces )
		{
#pragma omp critical
			{
				need_adjacent_faces();
			}
		}
#if 1
		std::vector< int > neighborTriangles;
		for( size_t i=0 ; i<_af[vert].size() ; i++ )
		{
			for( int j=0 ; j<3 ; j++ )
			{
				int _vert = _f[ _af[vert][i] ][j];
				for( size_t k=0 ; k<_af[_vert].size() ; k++ ) neighborTriangles.push_back( _af[_vert][k] );
			}
		}
		for( size_t i=0; i<neighborTriangles.size() ; i++ )				// Iterate over all faces adjacent to the nearest vertex
		{
			if( neighborTriangles[i]==fIndex ) continue;
			const TriangleIndex &face = _f[ neighborTriangles[i] ];		// The current face of interest
			BarycentricCoordinates( pIn , _v[ face[0] ] , _v[ face[1] ] , _v[ face[2] ] , a[0] , a[1] , a[2] );
	
			Real sum = 0;
			for( int d=0 ; d<3 ; d++ )
			{
				if( a[d]<0 ) a[d] = 0;
				sum += a[d];
			}
			for( int d=0 ; d<3 ; d++ ) a[d] /= sum;
			Point3D< Real > vertex = _v[face[0]] * a[0] + _v[face[1]] * a[1] +  _v[face[2]] * a[2];
			Real dist = Point3D< Real >::SquareNorm( vertex - pIn );
			if( fIndex<0 || dist<dst )
			{
				fIndex = (int) neighborTriangles[i];
				dst = dist;
				pOut = vertex;
			}
		}
#else
		for( size_t i=0; i<_af[vert].size() ; i++ )			// Iterate over all faces adjacent to the nearest vertex
		{
			const TriangleIndex &face = _f[_af[vert][i]];		// The current face of interest
			BarycentricCoordinates( pIn , _v[ face[0] ] , _v[ face[1] ] , _v[ face[2] ] , a[0] , a[1] , a[2] );
	
			Real sum = 0;
			for( int d=0 ; d<3 ; d++ )
			{
				if( a[d]<0 ) a[d] = 0;
				sum += a[d];
			}
			for( int d=0 ; d<3 ; d++ ) a[d] /= sum;
			Point3D< Real > vertex = _v[face[0]] * a[0] + _v[face[1]] * a[1] +  _v[face[2]] * a[2];
			Real dist = Point3D< Real >::SquareNorm( vertex - pIn );
			if( fIndex<0 || dist<dst )
			{
				fIndex = (int) _af[vert][i];
				dst = dist;
				pOut = vertex;
			}
		}
#endif
		return fIndex;
	}

	// [IN] vert: The index of the nearest point to pIn
	// [IN] pIn:  The point of interest
	// [OUT] a:   The barycentric coordinates of the point on the face
	// [Return]   The index of the face
	template< bool TestAdjacentFaces > int closest_point_on_face( int vert , Point3D< Real > pIn , Real a[3] ){ Point3D< Real > pOut ; return closest_point_on_face< TestAdjacentFaces >( vert , pIn , pOut , a ); }

	// Returns the closest face on the mesh when the nearest vertex is unknown

	// [IN] pIn:   The point of interest
	// [OUT] pOut: The point on the triangle that is closest
	// [OUT] a:    The barycentric coordinates of the point on the face
	// [Return]    The index of the face
	template< bool TestAdjacentFaces > int closest_point_on_face( Point3D< Real > pIn , Point3D< Real > &pOut , Real a[3]	){ return closest_point_on_face< TestAdjacentFaces >( closest_vertex_index( pIn ) , pIn , pOut , a ); }

	// [IN] pIn: The point of interest
	// [OUT] a:  The barycentric coordinates of the point on the face
	// [Return]  The index of the face
	template< bool TestAdjacentFaces > int closest_point_on_face( Point3D< Real > &pIn , Real a[3] ){ Point3D< Real > pOut ; return closest_point_on_face< TestAdjacentFaces >( pIn , pOut , a ); }
};
