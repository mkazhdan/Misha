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

#ifndef INTERVAL_MAP_INCLUDED
#define INTERVAL_MAP_INCLUDED

#include "Misha/Polynomial.h"

namespace MishaK
{
	// A C-K function that:
	// -- Evaluates to 0 for negative values,
	// -- Evaluates to 1 for values bigger than r,
	template< unsigned int K >
	struct IntervalMap
	{
		static const unsigned int Degree = 2*K+1;
		IntervalMap( double r ) : _r(r) , _P(Poly())
		{
			double scl = 1.;
			for( unsigned int d=0 ; d<=Degree ; d++ ) 
			{
				_P.coefficient(d) *= scl;
				scl /= _r;
			}
			Polynomial::Polynomial1D< Degree > __P = _P;
			for( unsigned int k=0 ; k<=K ; k++ )
			{
				std::cout << __P(0) << " : " << __P(_r) << std::endl;
				__P = __P.d(0);
			}
		}

		double operator()( double x ) const
		{
			if     ( x<0  ) return 0.;
			else if( x>_r ) return 1.;
			else            return _P(x);
		}

		static size_t Choose( unsigned int d , unsigned int k )
		{
			size_t c = 1;
			for( unsigned int _k=0 ; _k<k ; _k++ ) c *= (d-_k);
			return c;
		}

		static Polynomial::Polynomial1D< Degree > Poly( void )
		{
			Polynomial::Polynomial1D< Degree > P;

			MishaK::SquareMatrix< double , Degree+1 > M;
			Point< double , Degree+1 > x , b;
			b[0] = 0 , b[1] = 1;

			for( unsigned int d=0 ; d<=Degree ; d++ ) for( unsigned int k=0 ; k<=K ; k++ )
			{
				// The evaluation of k-th derivative of the d-th term at x = 0
				// => Choose( d , k ) x ^{d-k}
				M( d , 2*k+0 ) = (double)Choose( d , k ) * ( d==k ? 1 : 0 );
				M( d , 2*k+1 ) = (double)Choose( d , k );
			}

			x = M.inverse() * b;
			for( unsigned int d=0 ; d<=Degree ; d++ ) P.coefficient( d ) = x[d];
			return P;
		}

	protected:
		double _r;
		Polynomial::Polynomial1D< Degree > _P;
	};
}
#endif // INTERVAL_MAP_INCLUDED