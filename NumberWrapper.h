/*
Copyright (c) 2019, Michael Kazhdan
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

#ifndef NUMBER_WRAPPER_INCLUDED
#define NUMBER_WRAPPER_INCLUDED

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <limits>
#include <type_traits>
#include <Misha/UnionFind.h>

#define DEBUG_GRAPH

#define NEW_GRAPH

namespace Misha
{
	struct _EmptyNumberWrapper{};
	template< typename Number , typename Type=_EmptyNumberWrapper , size_t I=0 >
	struct NumberWrapper
	{
		NumberWrapper( Number _n=0 ) : n(_n){}
		Number n;
		NumberWrapper operator + ( NumberWrapper _n ) const { return NumberWrapper( n + _n.n ); }
		NumberWrapper operator - ( NumberWrapper _n ) const { return NumberWrapper( n - _n.n ); }
		NumberWrapper operator * ( NumberWrapper _n ) const { return NumberWrapper( n * _n.n ); }
		NumberWrapper operator / ( NumberWrapper _n ) const { return NumberWrapper( n / _n.n ); }
		NumberWrapper &operator += ( NumberWrapper _n ){ n += _n.n ; return *this; }
		NumberWrapper &operator -= ( NumberWrapper _n ){ n -= _n.n ; return *this; }
		NumberWrapper &operator *= ( NumberWrapper _n ){ n *= _n.n ; return *this; }
		NumberWrapper &operator /= ( NumberWrapper _n ){ n /= _n.n ; return *this; }
		bool operator == ( NumberWrapper _n ) const { return n==_n.n; }
		bool operator != ( NumberWrapper _n ) const { return n!=_n.n; }
		bool operator <  ( NumberWrapper _n ) const { return n<_n.n; }
		bool operator >  ( NumberWrapper _n ) const { return n>_n.n; }
		bool operator <= ( NumberWrapper _n ) const { return n<=_n.n; }
		bool operator >= ( NumberWrapper _n ) const { return n>=_n.n; }
		NumberWrapper operator ++ ( int ) { NumberWrapper _n(n) ; n++ ; return _n; }
		NumberWrapper operator -- ( int ) { NumberWrapper _n(n) ; n-- ; return _n; }
		NumberWrapper &operator ++ ( void ) { n++ ; return *this; }
		NumberWrapper &operator -- ( void ) { n-- ; return *this; }
		explicit operator Number () const { return n; }
	};

	template< typename Number , typename Type , size_t I >
	struct std::hash< NumberWrapper< Number , Type , I > >
	{
		size_t operator()( Misha::NumberWrapper< Number , Type , I > n) const { return std::hash< Number >{}( n.n ); }
	};

	template< typename Data , typename _NumberWrapper >
	struct VectorWrapper : public std::vector< Data >
	{
		VectorWrapper( void ){}
		VectorWrapper( size_t sz ) : std::vector< Data >( sz ){}
		VectorWrapper( size_t sz , Data d ) : std::vector< Data >( sz , d ){}

		std::vector< Data >::reference operator[]( _NumberWrapper n ){ return std::vector< Data >::operator[]( n.n ); }
		std::vector< Data >::const_reference operator[]( _NumberWrapper n ) const { return std::vector< Data >::operator[]( n.n ); }
	};

	template< typename Key , typename Value , typename _NumberWrapper >
	struct UnorderedMapWrapper : public std::unordered_map< Key , Value , typename _NumberWrapper::Hash >{};
}
#endif // NUMBER_WRAPPER_INCLUDED