/*
Copyright (c) 2026, Michael Kazhdan and Matthew Bolitho
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



//////////////////
// MultiBlender //
//////////////////
template< typename BaseBlender >
double MultiBlender< BaseBlender >::operator()( const double values[/*BaseBlender::N*/] , double x ) const
{
	double value = 0;
	for( unsigned int n=0 ; n<BaseBlender::N ; n++ ) value += static_cast< const BaseBlender & >( *this )._blendingFunctions[n](x) * values[n];
	return value;
}


template< typename BaseBlender >
template< unsigned int Dim >
double MultiBlender< BaseBlender >::operator()( const double values[/*BaseBlender::N^^Dim*/] , Point< double , Dim > p ) const
{
	return _evaluate< Dim >( values , &p[0] );
}

template< typename BaseBlender >
template< unsigned int Dim >
double MultiBlender< BaseBlender >::_evaluate( const double values[/*BaseBlender::N^^Dim*/] , const double p[/*Dim*/] ) const
{
	if constexpr( Dim==1 ) return static_cast< const BaseBlender & >( *this )( values , p[0] );
	else
	{
		double _values[BaseBlender::N];
		unsigned int offset = 1;
		for( unsigned int d=0 ; d<Dim-1 ; d++ ) offset *= BaseBlender::N;

		for( unsigned int i=0 ; i<BaseBlender::N ; i++ ) _values[i] = _evaluate< Dim-1 >( values+i*offset , p+1 );
		return static_cast< const BaseBlender & >( *this )( _values , p[0] );
	}
}

/////////////////////////
// ConstantInterpolant //
/////////////////////////

inline ConstantInterpolant::ConstantInterpolant( void )
{
	_blendingFunctions[0].coefficient( 0 ) =  1.;
}

#if 0
inline double ConstantInterpolant::operator()( const double values[/*N*/] , double x ) const
{
	return _blendingFunctions[0](x) * values[0];
}
#endif

///////////////////////
// LinearInterpolant //
///////////////////////
inline LinearInterpolant::LinearInterpolant( void )
{
	_blendingFunctions[0].coefficient( 0 ) =  1.;
	_blendingFunctions[0].coefficient( 1 ) = -1.;

	_blendingFunctions[1].coefficient( 0 ) =  0.;
	_blendingFunctions[1].coefficient( 1 ) =  1.;
}

#if 0
inline double LinearInterpolant::operator()( const double values[/*N*/] , double x ) const
{
	return _blendingFunctions[0](x) * values[0] + _blendingFunctions[1](x) * values[1];
}
#endif

///////////////////////////
// CatmullRomInterpolant //
///////////////////////////
inline CatmullRomInterpolant::CatmullRomInterpolant( void )
{
	_blendingFunctions[0].coefficient( 0 ) =  0.0;
	_blendingFunctions[0].coefficient( 1 ) = -0.5;
	_blendingFunctions[0].coefficient( 2 ) =  1.0;
	_blendingFunctions[0].coefficient( 3 ) = -0.5;

	_blendingFunctions[1].coefficient( 0 ) =  1.0;
	_blendingFunctions[1].coefficient( 1 ) =  0.0;
	_blendingFunctions[1].coefficient( 2 ) = -2.5;
	_blendingFunctions[1].coefficient( 3 ) =  1.5;

	_blendingFunctions[2].coefficient( 0 ) =  0.0;
	_blendingFunctions[2].coefficient( 1 ) =  0.5;
	_blendingFunctions[2].coefficient( 2 ) =  2.0;
	_blendingFunctions[2].coefficient( 3 ) = -1.5;

	_blendingFunctions[3].coefficient( 0 ) =  0.0;
	_blendingFunctions[3].coefficient( 1 ) =  0.0;
	_blendingFunctions[3].coefficient( 2 ) = -0.5;
	_blendingFunctions[3].coefficient( 3 ) =  0.5;
}

#if 0
inline double CatmullRomInterpolant::operator()( const double values[/*N*/] , double x ) const
{
	return _blendingFunctions[0](x) * values[0] + _blendingFunctions[1](x) * values[1] + _blendingFunctions[2](x) * values[2] + _blendingFunctions[3](x) * values[3];
}
#endif

/////////////////////////////
// UniformCubicApproximant //
/////////////////////////////
inline UniformCubicApproximant::UniformCubicApproximant( void )
{
	_blendingFunctions[0].coefficient( 0 ) =  1./6;
	_blendingFunctions[0].coefficient( 1 ) = -1./2;
	_blendingFunctions[0].coefficient( 2 ) =  1./2;
	_blendingFunctions[0].coefficient( 3 ) = -1./6;

	_blendingFunctions[1].coefficient( 0 ) =  2./3;
	_blendingFunctions[1].coefficient( 1 ) =  0.;
	_blendingFunctions[1].coefficient( 2 ) = -1.;
	_blendingFunctions[1].coefficient( 3 ) =  1./2;

	_blendingFunctions[2].coefficient( 0 ) =  1./6;
	_blendingFunctions[2].coefficient( 1 ) =  1./2;
	_blendingFunctions[2].coefficient( 2 ) =  1./2;
	_blendingFunctions[2].coefficient( 3 ) = -1./2;

	_blendingFunctions[3].coefficient( 0 ) =  0.;
	_blendingFunctions[3].coefficient( 1 ) =  0.;
	_blendingFunctions[3].coefficient( 2 ) =  0.;
	_blendingFunctions[3].coefficient( 3 ) =  1./6;
}

#if 0
inline double UniformCubicApproximant::operator()( const double values[/*N*/] , double x ) const
{
	return _blendingFunctions[0](x) * values[0] + _blendingFunctions[1](x) * values[1] + _blendingFunctions[2](x) * values[2] + _blendingFunctions[3](x) * values[3];
}
#endif

/////////////
// BSpline //
/////////////
template< unsigned int Degree >
BSpline< Degree >::BSpline( void )
{
	for( unsigned int n=0 ; n<N ; n++ ) _blendingFunctions[n] = _BSplineComponent( n );
}

#if 0
template< unsigned int Degree >
double BSpline< Degree >::operator()( const double values[/*N*/] , double x ) const
{
	double value = 0;
	for( unsigned int n=0 ; n<N ; n++ ) value += _blendingFunctions[n](x) * values[n];
	return value;
}
#endif

template< unsigned int Degree >
Polynomial::Polynomial1D< Degree > BSpline< Degree >::_BSplineComponent( unsigned int i )
{
	Polynomial::Polynomial1D< Degree > p;
	if constexpr( Degree==0 ) p.coefficient(0) = 1.;
	else
	{
		if( i<Degree )
		{
			Polynomial::Polynomial1D< Degree > _p = Polynomial::Integral( BSpline< Degree-1 >::_BSplineComponent( i ) );
			p -= _p;
			p.coefficient(0) += _p(1);
		}
		if( i>0 )
		{
			Polynomial::Polynomial1D< Degree > _p = Polynomial::Integral( BSpline< Degree-1 >::_BSplineComponent( i-1 ) );
			p += _p;
		}
	}
	return p;
}
