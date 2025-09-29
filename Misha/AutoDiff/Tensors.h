/*
Copyright (c) 2020, Michael Kazhdan
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

#ifndef TENSORS_INCLUDED
#define TENSORS_INCLUDED

#define NEW_TENSOR_CODE

#include <iostream>
#include "Misha/Algebra.h"
#include "Misha/ParameterPack.h"
#include "Misha/MultiDimensionalArray.h"
#include "Misha/Geometry.h"

namespace MishaK
{
	namespace AutoDiff
	{
		// A representation of a tensor, with dimensions encoded by a ParameterPack::UIntPack< ... > parameter
		template< typename Pack > struct PTensor;

		// A representation of a tensor, with dimensions encoded by the dimensions
		template< unsigned int ... Dims > using Tensor = PTensor< ParameterPack::UIntPack< Dims ... > >;


		// A zero-tensor is the same as a double value
		template<>
		struct PTensor< ParameterPack::UIntPack<> > : public InnerProductSpace< double , PTensor< ParameterPack::UIntPack<> > >
		{
			using Pack = ParameterPack::UIntPack<>;
			static const unsigned int Size = Pack::Size;
			static const unsigned int Dimension = 1;
			using Strides = ParameterPack::UIntPack<>;

			double data;

			double & operator()( void ){ return data; }
			const double & operator()( void ) const { return data; }


			PTensor( double d=0 ) : data(d) {}
			PTensor( Point< double , 1 > p ) : data(p[0]) {}

#if 1
#pragma message( "[WARNING] Should these be explicit?" )
			explicit operator       double &( void )       { return data; }
			explicit operator const double &( void ) const { return data; }
#else
			operator       double &( void )       { return data; }
			operator const double &( void ) const { return data; }
#endif

			void Add( const PTensor &t ){ data += t.data; }
			void Scale( double s ){ data *= s; }
			double InnerProduct( const PTensor &t ) const { return data * t.data; }
			template< unsigned int ... _Dims >
			PTensor< ParameterPack::UIntPack< _Dims ... > > operator * ( const PTensor< ParameterPack::UIntPack< _Dims ... > > &t ) const { return t * data; }
			template< unsigned int I , unsigned int ... _Dims >
			PTensor< ParameterPack::UIntPack< _Dims ... > > contractedOuterProduct( const PTensor< ParameterPack::UIntPack< _Dims ... > > &t ) const 
			{
				static_assert( I==0 , "[ERROR] Contraction suffix/prefix don't match" );
				return *this * t;
			}

			// Permute indices
			template< unsigned int ... PermutationValues >
			PTensor< ParameterPack::Permutation< Pack , ParameterPack::UIntPack< PermutationValues ... > > > permute( ParameterPack::UIntPack< PermutationValues ... > ) const
			{
				static_assert( sizeof ... ( PermutationValues ) == Size || ( sizeof ... ( PermutationValues ) == 0 && Size==1 ), "[ERROR] Permutation size doesn't match dimension" );
				return *this;
			}

			static PTensor Identity( void ){ return PTensor( 1. ); }

			friend std::ostream &operator << ( std::ostream &os , const PTensor &t ){ return os << t.data; }
		};

		// A general tensor
		template< unsigned int Dim , unsigned int ... Dims >
		struct PTensor< ParameterPack::UIntPack< Dim , Dims ... > > : public MultiDimensionalArray::Array< double , Dim , Dims ... > , public InnerProductSpace< double , PTensor< ParameterPack::UIntPack< Dim , Dims ... > > >
		{
			using MultiDimensionalArray::Array< double , Dim , Dims ... >::data;
			using Pack = ParameterPack::UIntPack< Dim , Dims ... >;
			static const unsigned int Size = Pack::Size;
			static const unsigned int Dimension = Dim * PTensor< ParameterPack::UIntPack< Dims ... > >::Dimension;
			using Strides = ParameterPack::Concatenation< ParameterPack::UIntPack< Tensor< Dims... >::Dimension > , typename Tensor< Dims ... >::Strides >;

			MultiDimensionalArray::     ArrayWrapper< double , Dim , Dims ... > operator()( void )       { return MultiDimensionalArray::Array< double , Dim , Dims ... >::operator()(); }
			MultiDimensionalArray::ConstArrayWrapper< double , Dim , Dims ... > operator()( void ) const { return MultiDimensionalArray::Array< double , Dim , Dims ... >::operator()(); }

			PTensor( void ){ memset( MultiDimensionalArray::Array< double , Dim , Dims ... >::data , 0 , sizeof( double ) * MultiDimensionalArray::ArraySize< Dim , Dims ... >() ); }

			static PTensor RightConcatenation( const PTensor< typename Pack::Transpose::Rest::Transpose > tensors[Pack::Last] )
			{
				using _Pack = Pack::Transpose::Rest::Transpose;
				PTensor t;
				if constexpr( std::is_same_v< _Pack , ParameterPack::UIntPack<> > ) for( unsigned int d=0 ; d<Pack::Last ; d++ ) t.data[d] = tensors[d].data;
				else for( unsigned int d=0 ; d<Pack::Last ; d++ ) for( unsigned int sz=0 ; sz<PTensor< _Pack >::Dimension ; sz++ ) t.data[sz*Pack::Last+d] = tensors[d].data[sz];
				return t;
			}

			static PTensor LeftConcatenation( const PTensor< typename Pack::Rest > tensors[Pack::First] )
			{
				using _Pack = Pack::Rest;
				static const unsigned int Size = PTensor< _Pack >::Dimension;
				PTensor t;
				if constexpr( std::is_same_v< _Pack , ParameterPack::UIntPack<> > )	for( unsigned int d=0 ; d<Pack::First ; d++ ) t.data[d] = tensors[d].data;
				else for( unsigned int d=0 ; d<Pack::First ; d++ ) for( unsigned int sz=0 ; sz<Size ; sz++ ) t.data[sz+d*Size] = tensors[d].data[sz];
				return t;
			}

			template< unsigned int _Dim , typename std::enable_if_t< Pack::Size==1 && ParameterPack::Comparison< ParameterPack::UIntPack< _Dim > , Pack >::Equal > * = nullptr >
			PTensor( Point< double , _Dim > t ){ for( unsigned int d=0 ; d<_Dim ; d++ ) operator()(d) = t[d]; }

			template< unsigned int Cols , unsigned int Rows , typename std::enable_if_t< Pack::Size==2 && ParameterPack::Comparison< ParameterPack::UIntPack< Rows , Cols > , Pack >::Equal > * = nullptr >
			PTensor( Matrix< double , Cols , Rows > value ){ for( unsigned int c=0 ; c<Cols ; c++ ) for( unsigned int r=0 ; r<Rows ; r++ ) operator()(r,c) = value(c,r); }

#pragma message( "[WARNING] Should these be explicit?" )
			template< unsigned int _Dim , typename std::enable_if_t< ParameterPack::Comparison< ParameterPack::UIntPack< _Dim > , Pack >::Equal > * = nullptr >
			explicit operator Point< double , _Dim > ( void ) const
			{
				Point< double , _Dim > p;
				for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = operator()(d);
				return p;
			}

			template< unsigned int Cols , unsigned int Rows , typename std::enable_if_t< ParameterPack::Comparison< ParameterPack::UIntPack< Rows , Cols > , Pack >::Equal > * = nullptr >
			explicit operator Matrix< double , Cols , Rows > ( void ) const
			{
				Matrix< double , Cols , Rows > m;
				for( unsigned int c=0 ; c<Cols ; c++ ) for( unsigned int r=0 ; r<Rows ; r++ ) m(c,r) = operator()(r,c);
				return m;
			}

			template< typename ... UInts >
			double &operator()( unsigned int index , UInts ... indices )
			{
				static_assert( sizeof...(indices)==Pack::Size-1 , "[ERROR] Wrong number of indices" );
				unsigned int idx[] = { index , indices ... };
				return MultiDimensionalArray::Array< double , Dim , Dims ... >::operator()( idx );
			}

			template< typename ... UInts >
			const double &operator()( unsigned int index , UInts ... indices ) const
			{
				static_assert( sizeof...(indices)==Pack::Size-1 , "[ERROR] Wrong number of indices" );
				unsigned int idx[] = { index , indices ... };
				return MultiDimensionalArray::Array< double , Dim , Dims ... >::operator()( idx );
			}

			double &operator()( const unsigned int indices[] ){ return MultiDimensionalArray::Array< double , Dim , Dims ... >::operator()( indices ); }
			const double &operator()( const unsigned int indices[] ) const { return MultiDimensionalArray::Array< double , Dim , Dims ... >::operator()( indices ); }

			// Inner-product space methods
			void Add( const PTensor &t ){ for( unsigned int i=0 ; i<Dimension ; i++ ) data[i] += t.data[i]; }
			void Scale( double s ){ for( unsigned int i=0 ; i<Dimension ; i++ ) data[i] *= s; }
			double InnerProduct( const PTensor &t ) const
			{
				double innerProduct = 0;
				for( unsigned int i=0 ; i<Dimension ; i++ ) innerProduct += data[i] * t.data[i];
				return innerProduct;
			}

			static PTensor< ParameterPack::UIntPack< Dim , Dims ... , Dim , Dims ... > > Identity( void )
			{
				static const unsigned int Size = sizeof ... ( Dims ) + 1;
				PTensor< ParameterPack::UIntPack< Dim , Dims ... , Dim , Dims ... > > id;
				unsigned int indices[ 2*Size ];
				MultiDimensionalArray::Loop< Size >::Run
				(
					ParameterPack::IsotropicUIntPack< Size >::Values , ParameterPack::UIntPack< Dim , Dims ... >::Values ,
					[&]( int d , int i ){ indices[d] = indices[ d+Size ] = i;} ,
					[&]( void ){ id( indices ) = 1.; }
				);
				return id;
			}

#ifdef NEW_TENSOR_CODE
#else // !NEW_TENSOR_CODE
			template< unsigned int ... PermutationValues >
			static auto PermutationTensor( ParameterPack::UIntPack< PermutationValues ... > )
			{
#pragma message( "[WARNING] Should avoid using PermutationTensor" )
				MK_WARN_ONCE( "Invoking PermutationTensor" );
				PTensor< ParameterPack::Concatenation< ParameterPack::Permutation< Pack , ParameterPack::UIntPack< PermutationValues ... > > , Pack > > t;
				const unsigned int permutation[] = { PermutationValues ... };
				unsigned int idx[ 2*Size ];
				MultiDimensionalArray::Loop< Size >::Run
				(
					ParameterPack::IsotropicUIntPack< Size >::Values , Pack::Values ,
					[&]( int d , int i ){ idx[ permutation[d] ] = idx[ Size+d ] = i; } ,
					[&]( void ){ t( idx ) = 1; }
				);
				return t;
			}
#endif // NEW_TENSOR_CODE

			// Permute indices
			template< unsigned int ... PermutationValues >
			PTensor< ParameterPack::Permutation< Pack , ParameterPack::UIntPack< PermutationValues ... > > > permute( ParameterPack::UIntPack< PermutationValues ... > ) const
			{
				static_assert( sizeof ... ( PermutationValues ) == Size , "[ERROR] Permutation size doesn't match dimension" );
				if constexpr( Pack::Size==1 ) return *this;
				else
				{
					using PPack = ParameterPack::UIntPack< PermutationValues ... >;
					using OutPack = ParameterPack::Permutation< Pack , PPack >;
					using OutStrides = ParameterPack::Permutation< Strides , PPack >;

					PTensor< ParameterPack::Permutation< Pack , PPack > > t;
					unsigned int inOffsets[ Size ];
					inOffsets[0] = 0;
					MultiDimensionalArray::Loop< Size-1 >::Run
					(
						ParameterPack::IsotropicUIntPack< Size >::Values , OutPack::Values ,
						[&]( int d , int i ){ inOffsets[d+1] = inOffsets[d] + OutStrides::Values[d] * i; } ,
						[&]( MultiDimensionalArray::ArrayWrapper< double , OutPack::Last > out )
						{
							const double * in = data + inOffsets[Size-1];
							unsigned int stride = OutStrides::Values[Size-1];
							for( unsigned int i=0 , ii=0 ; i<OutPack::Last ; i++ , ii+=stride ) out[i] = in[ii];
						} ,
						t
					);
					return t;
				}
			}

			// Extract slice
			template< unsigned int I >
			auto extract( const unsigned int indices[/*I*/] ) const
			{
				using Remainder = typename ParameterPack::Partition< I , Pack >::Second;
				PTensor< Remainder> t;

				if constexpr( Remainder::Size!=0 )
				{
					unsigned int _indices[ Pack::Size ];
					for( unsigned int i=0 ; i<I ; i++ ) _indices[i] = indices[i];

					MultiDimensionalArray::Loop< Remainder::Size >::Run
					(
						ParameterPack::IsotropicUIntPack< Remainder::Size >::Values , Remainder::Values ,
						[&]( int d , int i ){ _indices[d+I] = i; } ,
						[&]( double &_t ){ _t = operator()( _indices ); } ,
						t
					);
				}
				else static_cast< double & >( t ) = operator()( indices );
				return t;
			}

#ifdef NEW_TENSOR_CODE
#else // !NEW_TENSOR_CODE
			static auto TransposeTensor( void )
			{
#pragma message( "[WARNING] Should avoid using TransposeTensor" )
				MK_WARN_ONCE( "Invoking TransposeTensor" );
				PTensor< ParameterPack::Concatenation< typename Pack::Transpose , Pack > > t;
				unsigned int idx[ 2*Size ];
				MultiDimensionalArray::Loop< Size >::Run
				(
					ParameterPack::IsotropicUIntPack< Size >::Values , Pack::Transpose::Values ,
					[&]( int d , int i ){ idx[d] = idx[ 2*Size - 1 - d ] = i; } ,
					[&]( void ){ t( idx ) = 1; }
				);
				return t;
			}
#endif // NEW_TENSOR_CODE

			// Transpose operator
			PTensor< typename Pack::Transpose > transpose( void ) const
			{
				return permute( ParameterPack::SequentialPack< unsigned int , Pack::Size >::Transpose() );
			}

			// Outer product
			template< unsigned int ... _Dims >
			PTensor< ParameterPack::Concatenation< Pack , ParameterPack::UIntPack< _Dims ... > > > operator * ( const PTensor< ParameterPack::UIntPack< _Dims ... > > &t ) const 
			{
				return this->template contractedOuterProduct< 0 >( t );
			}

			PTensor< Pack > operator * ( const PTensor< ParameterPack::UIntPack<> > &t ) const { return *this * t.data; }

#ifdef NEW_TENSOR_CODE
#else // !NEW_TENSOR_CODE
		protected:
			template< unsigned int D1 , unsigned int D2 >
			static auto _ContractionTensor( void )
			{
				static_assert( D1<D2 , "[ERROR] Contraction indices are the same" );
				static_assert( D1<Pack::Size , "[ERROR] First contraction index too large" );
				static_assert( D2<Pack::Size , "[ERROR] Second contraction index too large" );
				static_assert( ParameterPack::Selection< D1 , Pack >::Value==ParameterPack::Selection< D2 , Pack >::Value , "[ERROR] Contraction dimensions differ" );
				using OutPack = typename ParameterPack::Selection< D1 , typename ParameterPack::Selection< D2 , Pack >::Complement >::Complement;

				PTensor< ParameterPack::Concatenation< OutPack , Pack > > t;

				unsigned int index[ Pack::Size+OutPack::Size ];
				if constexpr( OutPack::Size==0 )
					for( unsigned int i=0 ; i<Pack::template Get<D1>() ; i++ )
					{
						index[D1] = index[D2] = i;
						t( index ) = 1;
					}
				else
				{
					unsigned int out2in[ OutPack::Size ];
					{
						unsigned int count = 0;
						for( unsigned int i=0 ; i<Pack::Size ; i++ ) if( i!=D1 && i!=D2 ) out2in[ count++ ] = i;
					}

					MultiDimensionalArray::Loop< OutPack::Size >::Run
					(
						ParameterPack::IsotropicUIntPack< OutPack::Size >::Values , OutPack::Values ,
						[&]( int d , int i ){ index[d] = index[ out2in[d] ] = i; } ,
						[&]( void )
					{
						for( unsigned int i=0 ; i<Pack::template Get<D1>() ; i++ )
						{
							index[ OutPack::Size+D1 ] = i;
							index[ OutPack::Size+D2 ] = i;
							t( index ) = 1;
						}
					}
					);
				}

				return t;
			}

		public:
			template< unsigned int D1 , unsigned int D2 >
			static auto ContractionTensor( void )
			{
#pragma message( "[WARNING] Should avoid using ContractionTensor" )
				MK_WARN_ONCE( "Invoking ContractionTensor" );
				if constexpr( D1<D2 ) return _ContractionTensor< D1 , D2 >();
				else                  return _ContractionTensor< D2 , D1 >();
			}
#endif // NEW_TENSOR_CODE

			// Tensor contraction
			template< unsigned int I1 , unsigned int I2 >
			auto contract( void ) const
			{
				static_assert( I1!=I2 , "[ERROR] Contraction indices must differ" );
				static_assert( Pack::template Get< I1 >()==Pack::template Get< I2 >() , "[ERROR] Contraction dimensions don't match" );
				static_assert( I1<Pack::Size && I2<Pack::Size , "[ERROR] Contraction indices out of bounds" );
				if constexpr( I2<I1 ) return this->template contract< I2 , I1 >();
				using OutPack = typename ParameterPack::Selection< I1 , typename ParameterPack::Selection< I2 , Pack >::Complement >::Complement;
				PTensor< OutPack > out;
				if constexpr( Pack::Size>2 )
				{
					unsigned int indices[ OutPack::Size ];
					MultiDimensionalArray::Loop< OutPack::Size >::Run
					(
						ParameterPack::IsotropicUIntPack< OutPack::Size >::Values , OutPack::Values ,
						[&]( int d , int i ){ indices[d] = i; } ,
						[&]( double &_out )
					{
						unsigned int _indices[ Pack::Size ];
						unsigned int idx=0;
						for( unsigned int i=0 ; i<Pack::Size ; i++ ) if( i!=I1 && i!=I2 ) _indices[i] = indices[idx++];
						_out = 0;
						for( unsigned int i=0 ; i<Pack::template Get< I1 >() ; i++ )
						{
							_indices[I1] = _indices[I2] = i;
							_out += operator()( _indices );
						}
					} ,
						out
					);
				}
				else
				{
					double &_out = static_cast< double & >( out );
					unsigned int _indices[2];
					for( unsigned int i=0 ; i<Pack::template Get<I1>() ; i++ )
					{
						_indices[I1] = _indices[I2] = i;
						_out += operator()( _indices );
					}
				}
				return out;
			}

			// In1 := [ N{1} , ... , N{I} , N{I+1} , ... , N{K} ]
			// In2 :=                     [ N{I+1} , ... , N{K} , N{K+1} , ... N{M} ]
			// Out := [ N{1} , ... , N{I}             ,           N{K+1} , ... N{M} ]
			template< unsigned int I , unsigned int ... _Dims >
			PTensor< ParameterPack::Concatenation< typename ParameterPack::Partition< Size-I , Pack >::First , typename ParameterPack::Partition< I , ParameterPack::UIntPack< _Dims ... > >::Second > > contractedOuterProduct( const PTensor< ParameterPack::UIntPack< _Dims ... > > &t ) const 
			{
				static_assert( ParameterPack::Comparison< typename ParameterPack::Partition< Size-I , Pack >::Second , typename ParameterPack::Partition< I , ParameterPack::UIntPack< _Dims ... > >::First >::Equal , "[ERROR] Contraction suffix/prefix don't match" );
				using _Pack = ParameterPack::UIntPack< _Dims ... >;
				static const unsigned int _Size = _Pack::Size;
				using P1 = typename ParameterPack::Partition< Size-I ,  Pack >:: First;
				using P2 = typename ParameterPack::Partition< Size-I ,  Pack >::Second;
				using P3 = typename ParameterPack::Partition<      I , _Pack >::Second;

				using In1SliceType = MultiDimensionalArray::ConstSliceType< P1::Size , double , Dim , Dims ... >;
				using In2SliceType = MultiDimensionalArray::ConstSliceType< P2::Size , double , _Dims ... >;
				// In the case that we are collapsing completely, out is of type PTensor< ParameterPack::UIntPack<> >
				// -- Then the first and last loops are trivial and we never access the contents of out using operator[]
				using OutBaseType = typename std::conditional< ParameterPack::Concatenation< P1 , P3 >::Size!=0 , double , PTensor< ParameterPack::UIntPack<> > >::type;
				using OutSliceType = std::conditional_t< P3::Size!=0 , typename MultiDimensionalArray::SliceType< P2::Size , double , _Dims ... > , OutBaseType & >;

				const PTensor<  Pack > &in1 = *this;
				const PTensor< _Pack > &in2 = t;
				PTensor< ParameterPack::Concatenation< P1 , P3 > > out;

				// Iterate over {1,...,I} of in1 and out
				MultiDimensionalArray::Loop< P1::Size >::Run
				(
					ParameterPack::IsotropicUIntPack< P1::Size >::Values , P1::Values ,
					[]( int d , int i ){} ,
					[&]( In1SliceType _in1 , OutSliceType _out )
					{
						// Iterate over {I,...,K} of in1 and in2
						MultiDimensionalArray::Loop< P2::Size >::Run
						(
							ParameterPack::IsotropicUIntPack< P2::Size >::Values , P2::Values ,
							[]( int d , int i ){} ,
							[&]( double __in1 , In2SliceType _in2 )
							{
								// Iterate over {K+1,...,M} of in2 and out
								MultiDimensionalArray::Loop< P3::Size >::Run
								(
									ParameterPack::IsotropicUIntPack< P3::Size >::Values , P3::Values ,
									[]( int d , int i ){} ,
									[&]( double __in2 , OutBaseType &_out_ ){ _out_ += __in1 * __in2; } ,
									_in2 , _out
								);
							} ,
							_in1 , in2
						);
					} ,
					in1 , out
				);
				return out;
			}

			template< unsigned int I >
			PTensor< Pack > contractedOuterProduct( const PTensor< ParameterPack::UIntPack<> > &t ) const { return *this * t; }
		};
		Tensor<> operator + ( double s , Tensor<> t ){ return Tensor<>( s + t.data ); }
		Tensor<> operator + ( Tensor<> t , double s ){ return Tensor<>( s + t.data ); }
	}
}
#endif // TENSORS_INCLUDED
