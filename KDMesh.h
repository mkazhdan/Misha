/*
Copyright (c) 2025, Michael Kazhdan
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

#ifndef KD_MESH_INCLUDED
#define KD_MESH_INCLUDED

#include <mutex>
#include "KDtree.h"
#include "Geometry.h"

namespace MishaK
{
	template< class Real >
	struct KDMesh
	{
	protected:
		const std::vector< Point3D< Real > >& _v;
		const std::vector< TriangleIndex >& _f;
		SzymonRusinkiewicz::KDtree< Real , 3 > _kd;
		mutable std::vector< std::vector< int > > _af;
	public:
		KDMesh( const std::vector< Point3D< Real > >& v , const std::vector< TriangleIndex >& f ) : _v(v) , _f(f) , _kd( &v[0][0] , v.size() ){  }

		static void BarycentricCoordinates( const Point< Real , 3 >& p , const Point< Real , 3 >& v1 , const Point< Real , 3 >& v2, const Point< Real , 3 >& v3 , Real& a0 , Real& a1 , Real& a2 )
		{
			Point< Real , 3 > p0 =  p - v1;
			Point< Real , 3 > p1 = v2 - v1;
			Point< Real , 3 > p2 = v3 - v1;
			Point< Real , 3 >  n  = Point< Real , 3 >::CrossProduct( p1 , p2 );
			Point< Real , 3 > _v1 = Point< Real , 3 >::CrossProduct( p2 , n  );
			Point< Real , 3 > _v2 = Point< Real , 3 >::CrossProduct( p1 , n  );

			_v1 /= Point< Real , 3 >::Dot( _v1 , p1 );
			_v2 /= Point< Real , 3 >::Dot( _v2 , p2 );

			a1 = Point< Real , 3 >::Dot( _v1 , p0 );
			a2 = Point< Real , 3 >::Dot( _v2 , p0 );
			a0 = Real(1.0) - a1 - a2;
		}

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
				static std::mutex setAtomicMutex;
				std::lock_guard< std::mutex > lock( setAtomicMutex );
				need_adjacent_faces();
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
		template< bool TestAdjacentFaces > int closest_point_on_face( Point3D< Real > pIn , Point3D< Real > &pOut , Real a[3] ){ return closest_point_on_face< TestAdjacentFaces >( closest_vertex_index( pIn ) , pIn , pOut , a ); }

		// [IN] pIn: The point of interest
		// [OUT] a:  The barycentric coordinates of the point on the face
		// [Return]  The index of the face
		template< bool TestAdjacentFaces > int closest_point_on_face( Point3D< Real > &pIn , Real a[3] ){ Point3D< Real > pOut ; return closest_point_on_face< TestAdjacentFaces >( pIn , pOut , a ); }
	};
}
#endif // KD_MESH_INCLUDED