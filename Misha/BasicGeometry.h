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

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <concepts>

#include "Miscellany.h"
#include "Geometry.h"
#include "RegularGrid.h"
#include "IsoSurface3D.h"
#include "AutoDiff/EquationParser.h"

#ifdef USE_Q_HULL
#include <QHull/q_misha.h>
#endif // USE_Q_HULL

namespace MishaK
{
	namespace BasicGeometry
	{
		static const unsigned int K = 2;
		static const unsigned int Dim = K+1;

		using Mesh = std::tuple< std::vector< MishaK::Point< double , K > > , std::vector< MishaK::SimplexIndex< K > > >;
		using OrientedEmbeddedMesh = std::tuple< std::vector< MishaK::Point< double , Dim > > , std::vector< MishaK::Point< double , Dim > > , std::vector< MishaK::SimplexIndex< K > > >;

		template< typename XForm >
		concept HasXForm = requires( const XForm xForm , MishaK::Point< double , Dim > p ){ { xForm(p) } -> std::same_as< MishaK::Point< double , Dim > >; };

		/////////////////
		// Declaration //
		/////////////////

		template< HasXForm XForm >
		OrientedEmbeddedMesh Implicit( const MishaK::EquationParser::Node & F , unsigned int res , XForm && CubeToWorld , double isoValue );

		OrientedEmbeddedMesh Implicit( const MishaK::EquationParser::Node & F , unsigned int res , double isoValue );

		template< HasXForm XForm >
		OrientedEmbeddedMesh Implicit( std::string eqn_xyz , unsigned int res , XForm && CubeToWorld , double isoValue );

		OrientedEmbeddedMesh ImplicitSphere( unsigned int res , double radius , double scale=0.5 );
		OrientedEmbeddedMesh ImplicitTorus( unsigned int res , double radius1 , double radius2 , double scale=0.5 );

		OrientedEmbeddedMesh Tetrahedron( void );
		OrientedEmbeddedMesh Octahedron( void );
		OrientedEmbeddedMesh Iocsahedron( void );

		OrientedEmbeddedMesh EmbedTorus( const Mesh & mesh , double radius1 , double radius2 );

#ifdef USE_Q_HULL
		OrientedEmbeddedMesh RandomEmbeddedSphere( unsigned int numVertices , double radius=1. , MishaK::Point< double , Dim > scale=MishaK::Point< double , Dim >(1.,1.,1.) );
		Mesh RandomTorusWithParameterization( std::vector< MishaK::Simplex< double , 2 , 2 > > & pSimplices , unsigned int numVertices , MishaK::Point< double , K > scale=MishaK::Point< double , K >(1.,1.) );
		OrientedEmbeddedMesh RandomEmbeddedTorusWithParameterization( std::vector< MishaK::Simplex< double , 2 , 2 > > & pSimplices , unsigned int numVertices , double radius1 , double radius2 , MishaK::Point< double , K > scale=MishaK::Point< double , K >(1.,1.) ){ return EmbedTorus( std::forward< Mesh >( RandomTorusWithParameterization( pSimplices , numVertices , scale ) ) , radius1 , radius2 ); }
		Mesh RandomTorus( unsigned int numVertices , MishaK::Point< double , K > scale=MishaK::Point< double , K >(1.,1.) ){ std::vector< MishaK::Simplex< double , 2 , 2 > > pSimplices ; return RandomTorusWithParameterization( pSimplices , numVertices , scale ); }
		OrientedEmbeddedMesh RandomEmbeddedTorus( unsigned int numVertices , double radius1 , double radius2 , MishaK::Point< double , K > scale=MishaK::Point< double , K >(1.,1.) ){ return EmbedTorus( std::forward< Mesh >( RandomTorus( numVertices , scale ) ) , radius1 , radius2 ); }
#endif // USE_Q_HULL

		Mesh RegularTorusWithParameterization( std::vector< MishaK::Simplex< double , 2 , 2 > > & pSimplices , unsigned int res );
		OrientedEmbeddedMesh RegularEmbeddedTorusWithParameterization( std::vector< MishaK::Simplex< double , 2 , 2 > > & pSimplices , unsigned int res , double radius1 , double radius2 ){ return EmbedTorus( std::forward< Mesh >( RegularTorusWithParameterization( pSimplices , res ) ) , radius1 , radius2 ); }
		Mesh RegularTorus( unsigned int res ){ std::vector< MishaK::Simplex< double , 2 , 2 > > pSimplices ; return RegularTorusWithParameterization( pSimplices , res ); }
		OrientedEmbeddedMesh RegularEmbeddedTorus( unsigned int res , double radius1 , double radius2 ){ return EmbedTorus( std::forward< Mesh >( RegularTorus( res ) ) , radius1 , radius2 ); }

		struct Torus
		{
			Torus( double r1 , double r2 );

			MishaK::Point< double , 3 > operator()( double s , double t ) const;
			MishaK::Point< double , 3 > ds( double s , double t ) const;
			MishaK::Point< double , 3 > dt( double s , double t ) const;
			MishaK::Point< double , 3 > n( double s , double t ) const;
			MishaK::SquareMatrix< double , 2 > g( double s , double t ) const;
			MishaK::Point< double , 3 > operator()( MishaK::Point< double , 2 > st ) const;
			MishaK::Point< double , 2 > operator()( MishaK::Point< double , 3 > p ) const;
			MishaK::Point< double , 3 > ds( MishaK::Point< double , 2 > st ) const;
			MishaK::Point< double , 3 > dt( MishaK::Point< double , 2 > st ) const;
			MishaK::Point< double , 3 > n(  MishaK::Point< double , 2 > st ) const;
			MishaK::SquareMatrix< double , 2 > g( MishaK::Point< double , 2 > st ) const;
		protected:
			double _r1 , _r2;
		};
#include "BasicGeometry.inl"
	}
}

