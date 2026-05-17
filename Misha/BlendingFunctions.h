/*
Copyright (c) 2026, Michael Kazhdan
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

#ifndef BLENDING_FUNCTIONS_INCLUDED
#define BLENDING_FUNCTIONS_INCLUDED

#include <functional>
#include "Array.h"
#include "ParameterPack.h"
#include "MultiThreading.h"

namespace MishaK
{
	namespace BlendingFunctions
	{
		enum struct BlendType
		{
			CONSTANT ,
			LINEAR ,
			CATMULL_ROM ,
			UNIFORM_CUBIC ,
			B_SPLINE
		};
		inline const std::vector< std::string > BlendTypeNames = { "constant" , "linear" , "Catmull-Rom" , "uniform cubic" , "B-spline" };


		template< typename BaseBlender >
		struct MultiBlender
		{
			double operator()( const double values[/*BaseBlender::N*/] , double x ) const;

			template< unsigned int Dim >
			double operator()( const double values[/*BaseBlender::N^^Dim*/] , Point< double , Dim > p ) const;

		protected:
			template< unsigned int Dim >
			double _evaluate( const double values[/*BaseBlender::N^^Dim*/] , const double p[/*Dim*/] ) const;
		};

		template< unsigned int Dim , typename Blender >
		double MultiBlendValues( const Blender & blender , const double values[/*N^^Dim*/] , const double p[/*Dim*/] );

		struct ConstantInterpolant : public MultiBlender< ConstantInterpolant >
		{
			friend MultiBlender< ConstantInterpolant >;

			static const unsigned int N = 1;

			ConstantInterpolant( void );

		protected:
			Polynomial::Polynomial1D< N-1 > _blendingFunctions[N];
		};

		struct LinearInterpolant : public MultiBlender< LinearInterpolant >
		{
			friend MultiBlender< LinearInterpolant >;

			static const unsigned int N = 2;

			LinearInterpolant( void );

		protected:
			Polynomial::Polynomial1D< N-1 > _blendingFunctions[N];
		};

		struct CatmullRomInterpolant : public MultiBlender< CatmullRomInterpolant >
		{
			friend MultiBlender< CatmullRomInterpolant >;

			static const unsigned int N = 4;

			CatmullRomInterpolant( void );

		protected:
			Polynomial::Polynomial1D< N-1 > _blendingFunctions[N];
		};

		struct UniformCubicApproximant : public MultiBlender< UniformCubicApproximant >
		{
			friend MultiBlender< UniformCubicApproximant >;

			static const unsigned int N = 4;

			UniformCubicApproximant( void );

		protected:
			Polynomial::Polynomial1D< N-1 > _blendingFunctions[N];
		};

		template< unsigned int Degree >
		struct BSpline : MultiBlender< BSpline< Degree > >
		{
			template< unsigned int _Degree > friend struct BSpline;
			friend MultiBlender< BSpline< Degree > >;

			static const unsigned int N = Degree+1;

			BSpline( void );

		protected:
			Polynomial::Polynomial1D< N-1 > _blendingFunctions[N];

			static Polynomial::Polynomial1D< Degree > _BSplineComponent( unsigned int i );
		};

#include "BlendingFunctions.inl"
	}
}
#endif // BLENDING_FUNCTIONS_INCLUDED
