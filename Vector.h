/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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

#ifndef __VECTOR_HPP
#define __VECTOR_HPP

#include <cassert>
#include "Array.h"

namespace MishaK
{
	template<class T>
	class Vector
	{
	public:
		// For better compatibility with std::vector
		typedef T value_type;

		Vector();
		Vector( const Vector<T>& V );
		Vector( size_t N );
		Vector( size_t N, T* pV );
		~Vector();

		//	const T& operator () (size_t i) const;
		//	T& operator () (size_t i);
#if defined WIN32 || defined WIN64 // JP: none of this works on gcc: "ambiguous overloads"
#if 1 // This is a real pain, because depending on whether
		// one compiles with 32-bit or 64-bit, size-t is unsigned int or unsigned long long
		const T& operator[]( unsigned int       i ) const { return m_pV[i]; }
		T& operator[]( unsigned int       i )       { return m_pV[i]; }
		const T& operator[](          int       i ) const { return m_pV[i]; }
		T& operator[](          int       i )       { return m_pV[i]; }
#ifdef _WIN64
		const T& operator[]( unsigned long long i ) const { return m_pV[i]; }
		T& operator[]( unsigned long long i )       { return m_pV[i]; }
		const T& operator[](          long long i ) const { return m_pV[i]; }
		T& operator[](          long long i )       { return m_pV[i]; }
#endif // _WIN64
#else
		const T& operator[]( size_t i ) const { return m_pV[i]; }
		T& operator[]( size_t i )       { return m_pV[i]; }
		const T& operator[](    int i ) const { return m_pV[i]; }
		T& operator[](    int i )       { return m_pV[i]; }
#endif
#else // Mac
		const T& operator[] (size_t i) const { return m_pV[i]; }
		T& operator[] (size_t i)       { return m_pV[i]; }
#endif // Windows


		void SetZero();

		size_t size() const;
		void resize( size_t N );

		Vector operator * (const T& A) const;
		Vector operator / (const T& A) const;
		Vector operator - (const Vector& V) const;
		Vector operator + (const Vector& V) const;

		operator Pointer( T )( );
		operator ConstPointer( T )( ) const;

		Pointer( T ) operator + ( int idx );
		Pointer( T ) operator - ( int idx );
		ConstPointer( T ) operator + ( int idx ) const;
		ConstPointer( T ) operator - ( int idx ) const;

		Vector& operator *= (const T& A);
		Vector& operator /= (const T& A);
		Vector& operator += (const Vector& V);
		Vector& operator -= (const Vector& V);

		template<class Real>
		Vector& AddScaled(const Vector& V,const Real& scale);
		template<class Real>
		Vector& SubtractScaled(const Vector& V,const Real& scale);
		template<class Real>
		static void Add(const Vector& V1,const Real& scale1,const Vector& V2,const Real& scale2,Vector& Out);
		template<class Real>
		static void Add(const Vector& V1,const Real& scale1,const Vector& V2,Vector& Out);

		Vector operator - () const;

		Vector& operator = (const Vector& V);
		T Dot( const Vector& V ) const;
		template<class Real>
		Real IPSDot( const Vector& V ) const;
		T Length() const;

		/*
		T Norm( size_t Ln ) const;
		void Normalize();
		*/
		//	T* m_pV;
		Pointer( T ) m_pV;
	protected:
		size_t m_N;

	};

	template<class T,int Dim>
	class NVector
	{
	public:
		NVector();
		NVector( const NVector& V );
		NVector( size_t N );
		NVector( size_t N, T* pV );
		~NVector();

		const T* operator () (size_t i) const;
		T* operator () (size_t i);
		const T* operator [] (size_t i) const;
		T* operator [] (size_t i);

		void SetZero();

		size_t size() const;
		void resize( size_t N );

		NVector operator * (const T& A) const;
		NVector operator / (const T& A) const;
		NVector operator - (const NVector& V) const;
		NVector operator + (const NVector& V) const;

		NVector& operator *= (const T& A);
		NVector& operator /= (const T& A);
		NVector& operator += (const NVector& V);
		NVector& operator -= (const NVector& V);

		NVector& AddScaled(const NVector& V,const T& scale);
		NVector& SubtractScaled(const NVector& V,const T& scale);
		static void Add(const NVector& V1,const T& scale1,const NVector& V2,const T& scale2,NVector& Out);
		static void Add(const NVector& V1,const T& scale1,const NVector& V2,				NVector& Out);

		NVector operator - () const;

		NVector& operator = (const NVector& V);

		T Dot( const NVector& V ) const;

		T Length() const;

		T Norm( size_t Ln ) const;
		void Normalize();

		T* m_pV;
	protected:
		size_t m_N;

	};
#include "Vector.inl"
}
#endif
