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

#ifndef ORTHOGONAL_FRAME_INCLUDED
#define ORTHOGONAL_FRAME_INCLUDED

#include "Geometry.h"

namespace MishaK
{
	template< unsigned int Dim , typename ... > struct OrthogonalFrame;

	template< unsigned int Dim , typename Real , typename V >
	struct OrthogonalFrame< Dim , Real , V >
	{
		template< typename DotFunctor /* = std::function< Real ( V , V ) > */ , typename RandomVFunctor /* std::function< V (void) > */ , typename ... Vs >
		OrthogonalFrame( DotFunctor && Dot , RandomVFunctor && RandomV , const V v[] , unsigned int vNum )
		{
			static const Real EPS = static_cast< Real >( 1e-8 );
			_init( std::forward< DotFunctor >( Dot ) , std::forward< RandomVFunctor >( RandomV ) , v , vNum , EPS );
		}

#if 0
		template< typename DotFunctor /* = std::function< Real ( V , V ) > */ , typename RandomVFunctor /* std::function< V (void) > */ , typename ... Vs >
		OrthogonalFrame( DotFunctor && Dot , RandomVFunctor && RandomV , Vs ... p )
		{
			static_assert( sizeof...(Vs)<Dim , "[ERROR] Too many vectors" );

			static const Real EPS = static_cast< Real >( 1e-8 );
			V v[] = { p... };
			_init( std::forward< DotFunctor >( Dot ) , std::forward< RandomVFunctor >( RandomV ) , v , sizeof...(Vs) , EPS );
		}
#endif

		const V & operator[] ( unsigned int i ) const { return _v[i]; }

	protected:
		V _v[Dim];

		template< typename DotFunctor /* = std::function< Real ( V , V ) > */ , typename RandomVFunctor /* std::function< V (void) > */ >
		void _init( DotFunctor && Dot , RandomVFunctor && RandomV , const V v[] , unsigned int vNum , double eps )
		{
			static_assert( std::is_convertible_v< DotFunctor , std::function< Real ( V , V ) > > , "[ERROR] DotFunctor poorly formed" );
			static_assert( std::is_convertible_v< RandomVFunctor , std::function< V ( void ) > > , "[ERROR] RandomVFunctor poorly formed" );

			if( vNum>Dim ) MK_ERROR_OUT( "Too many vectors" );

			// Orthogonalize the given vectors
			for( unsigned int i=0 ; i<vNum ; i++ )
			{
				_v[i] = v[i];
				for( unsigned int j=0 ; j<i ; j++ ) _v[i] -= _v[j] * Dot( _v[i] , _v[j] );
				Real l2 = Dot( _v[i] , _v[i] );
				if( !l2 ) MK_THROW( "Degenerate frame" );
				_v[i] /= static_cast< Real >( sqrt( l2 ) );
			}

			// Complete to an orthonormal frame
			for( unsigned i=vNum ; i<Dim ; i++ )
			{
				while( true )
				{
					_v[i] = RandomV();
					for( unsigned int j=0 ; j<i ; j++ ) _v[i] -= _v[j] * Dot( _v[i] , _v[j] );
					Real l2 = Dot( _v[i] , _v[i] );
					if( l2>eps )
					{
						_v[i] /= static_cast< Real >( sqrt( l2 ) );
						break;
					}
				}
			}
		}
	};

	template< unsigned int Dim >
	struct OrthogonalFrame< Dim > : public OrthogonalFrame< Dim , double , Point< double , Dim > >
	{
		using V = Point< double , Dim >;

#if 1
		OrthogonalFrame( const Point< double , Dim > p[] , unsigned int num , bool orthonormal=true ) 
			: OrthogonalFrame< Dim , double , Point< double , Dim > >( []( V v1 , V v2 ){ return V::Dot(v1,v2); } , []( void ){ return RandomSpherePoint< double , Dim >(); } , p , num )
		{
			if( orthonormal )
			{
				SquareMatrix< double , Dim > O;
				for( unsigned int i=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) O(i,j) = OrthogonalFrame< Dim , double , Point< double , Dim > >::_v[i][j];
				if( O.determinant()<0 ) OrthogonalFrame< Dim , double , Point< double , Dim > >::_v[Dim-1] *= -1;
			}
		}

		static SquareMatrix< double  , Dim > RotationMatrix( Point< double , Dim > source , Point< double , Dim > target , double eps )
		{
			if( Point< double , Dim >::SquareNorm(source)<eps ) return SquareMatrix< double , Dim >::Identity();
			if( Point< double , Dim >::SquareNorm(target)<eps ) return SquareMatrix< double , Dim >::Identity();
			source /= Point< double , Dim >::Length( source );
			target /= Point< double , Dim >::Length( target );
			if( Point< double , Dim >::SquareNorm( source-target )<eps  ) return SquareMatrix< double , Dim >::Identity(); 

			Point< double , Dim > v[] = { source , target };
			OrthogonalFrame oFrame( v , 2 );
			double cs = Point3D< double >::Dot( target , oFrame[0] );
			double sn = Point3D< double >::Dot( target , oFrame[1] );

			// Takes:
			//		n[0] ->  cs * n[0] + sn * n[1]
			//		n[1] -> -sn * n[0] + cs * n[1]
			//		n[i] -> n[i] \forall i>1
			SquareMatrix< double , Dim > R =
				OuterProduct( oFrame[0] , oFrame[0] ) * ( cs) + OuterProduct( oFrame[1] , oFrame[0] ) * sn +
				OuterProduct( oFrame[0] , oFrame[1] ) * (-sn) + OuterProduct( oFrame[1] , oFrame[1] ) * cs;
			for( unsigned int i=2 ; i<Dim ; i++ ) R += OuterProduct( oFrame[i] , oFrame[i] );
			return R;
		}
#else
		template< typename ... Points >
		OrthogonalFrame( Points ... p ) 
			: OrthogonalFrame< Dim , double , Point< double , Dim > >( []( V v1 , V v2 ){ return V::Dot(v1,v2); } , []( void ){ return RandomSpherePoint< double , Dim >(); } , p... )
		{}
#endif
	};
}

#endif // ORTHOGONAL_FRAME_INCLUDED