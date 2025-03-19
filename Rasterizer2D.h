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

#ifndef RASTERIZER_2D_INCLUDED
#define RASTERIZER_2D_INCLUDED

#include "RegularGrid.h"
#include "Geometry.h"
#include "Miscellany.h"

#if 1 // NEW_CODE
#define NEW_RASTERIZER
#endif // NEW_CODE

#ifdef NEW_RASTERIZER
#else // !NEW_RASTERIZER
// [WARNING] Not sure that the un-Centered version works right.
#endif // NEW_RASTERIZER

namespace MishaK
{
	namespace Rasterizer2D
	{
		using Index = typename RegularGrid< 2 >::Index;
		template< unsigned int Dim > using Range = typename RegularGrid< Dim >::Range;
		using Triangle = Simplex< double , 2 , 2 >;

#ifdef NEW_RASTERIZER
		// CenterPoints = true:  Invokes the rasterization functor for all cells centers that are (left-bottom) inside of the triangle
		// CenterPoints = false: Invokes the rasterization functor for all cells overlapping the triangle
		template< bool CenterPoints , typename RasterizationFunctor /*=std::function< void ( Index ) > )*/ >
#else // !NEW_RASTERIZER
		// Invokes the rasterization functor for all cells centers that are (left-bottom) inside of the triangle
		template< bool Centered , typename RasterizationFunctor /*=std::function< void ( Index ) > )*/ >
#endif // NEW_RASTERIZER
		void Rasterize( Triangle triangle , RasterizationFunctor && F , Range< 2 > cellRange );

		//////////////////////////////////////////

#ifdef NEW_RASTERIZER
		template< bool CenterPoints >
#else // !NEW_RASTERIZER
		template< bool Centered >
#endif // NEW_RASTERIZER
		Range< 1 > GetCellRange( double s1 , double s2 )
		{
			Range< 1 > range;

#ifdef NEW_RASTERIZER
			if constexpr( CenterPoints )
			{
				// Solve for the smallest integer mn s.t.:
				//	mn+0.5 >= s1
				//	mn >= s1-0.5
				range.first[0] = (int)std::ceil( s1-0.5 );

				// Solve for the largest integer mx s.t.:
				//	mx+0.5 < s2
				//	mx < s2-0.5
				range.second[0] = (int)std::floor( s2-0.5 );
				if( (s2-0.5)==range.second[0] ) range.second[0]--;
			}
			else
			{
				range.first[0] = (int)std::floor( s1+0 );
				range.second[0] = (int)std::ceil( s2-1 );
			}
#else // !NEW_RASTERIZER
			if constexpr( Centered )
			{
				// Solve for the smallest integer mn s.t.:
				//	mn+0.5 >= s1
				//	mn >= s1-0.5
				range.first[0] = (int)std::ceil( s1-0.5 );

				// Solve for the largest integer mx s.t.:
				//	mx+0.5 < s2
				//	mx < s2-0.5
				range.second[0] = (int)std::floor( s2-0.5 );
				if( (s2-0.5)==range.second[0] ) range.second[0]--;
			}
			else
			{
				// Solve for the smallest integer mn s.t.:
				//	mn >= s1
				range.first[0] = (int)std::ceil( s1 );

				// Solve for the largest integer mx s.t.:
				//	mx < s2
				range.second[0] = (int)std::floor( s2 );
				if( s2==range.second[0] ) range.second[0]--;
			}
#endif // RASTERIZER

			return range;
		}

		// Rasterizes a horizontal line segment connected to a  point
#ifdef NEW_RASTERIZER
		template< bool CenterPoints , typename RasterizationFunctor /*=std::function< void ( Index ) > )*/ >
#else // !NEW_RASTERIZER
		template< bool Centered , typename RasterizationFunctor /*=std::function< void ( Index ) > )*/ >
#endif // NEW_RASTERIZER
		void _Rasterize( double y , double x0 , double x1 , Point< double , 2 > tip , RasterizationFunctor &&  F , const Range< 1 > cellRanges[2] )
		{
			double y0 = y , y1 = tip[1];
			if( y0>y1 ) std::swap( y0 , y1 );
			if( x0>x1 ) std::swap( x0 , x1 );

#ifdef NEW_RASTERIZER
			auto Intersection = [&]( double _y )
				{
					// Solve for s s.t.:
					//	y*(1-s) + tip[1]*s = _y
					//  (tip[1]-y)*s = _y - y
					//	s = (_y-y) / (tip[1]-y)
					double s = (_y-y) / (tip[1]-y);
					return std::pair< double , double >( x0*(1-s) + tip[0]*s , x1*(1-s) + tip[0]*s );
				};

			auto HorizontalCellRange = [&]( int iy )
				{
					if constexpr( CenterPoints )
					{
						std::pair< double , double > r = Intersection( iy+0.5 );
						return Range< 1 >::Intersect( cellRanges[0] , GetCellRange< CenterPoints >( r.first , r.second ) );
					}
					else
					{
						std::pair< double , double > r1 = Intersection( std::max< double >( y0 , iy+0. ) );
						std::pair< double , double > r2 = Intersection( std::min< double >( y1 , iy+1. ) );
						return Range< 1 >::Intersect( cellRanges[0] , GetCellRange< CenterPoints >( std::min< double >( r1.first , r2.first ) , std::max< double >( r1.second , r2.second ) ) );
					}
				};

			Range< 1 > iyRange = Range< 1 >::Intersect( cellRanges[1] , GetCellRange< CenterPoints >( y0 , y1 ) );

#else // !NEW_RASTERIZER
			auto Intersection = [&]( int iy )
				{
					if constexpr( Centered )
					{
						// Solve for s s.t.:
						//	y*(1-s) + tip[1]*s = iy+0.5
						//  (tip[1]-y)*s = iy + 0.5 - y
						//	s = (iy+0.5-y) / (tip[1]-y)
						double s = (iy+0.5-y) / (tip[1]-y);
						return std::pair< double , double >( x0*(1-s) + tip[0]*s , x1*(1-s) + tip[0]*s );
					}
					else
					{
						// Solve for s s.t.:
						//	y*(1-s) + tip[1]*s = iy
						//  (tip[1]-y)*s = iy - y
						//	s = (iy-y) / (tip[1]-y)
						double s = (iy-y) / (tip[1]-y);
						return std::pair< double , double >( x0*(1-s) + tip[0]*s , x1*(1-s) + tip[0]*s );
					}
				};

			auto HorizontalCellRange = [&]( int iy )
				{
					std::pair< double , double > r = Intersection( iy );
					return Range< 1 >::Intersect( cellRanges[0] , GetCellRange< Centered >( r.first , r.second ) );
				};

			Range< 1 > iyRange = Range< 1 >::Intersect( cellRanges[1] , GetCellRange< Centered >( y0 , y1 ) );
#endif // NEW_RASTERIZER

			for( int iy=iyRange.first[0] ; iy<=iyRange.second[0] ; iy++ )
			{
				Range< 1 > ixRange = HorizontalCellRange( iy );
				for( int ix=ixRange.first[0] ; ix<=ixRange.second[0] ; ix++ ) F( Index(ix,iy) );
			}
		}

#ifdef NEW_RASTERIZER
		template< bool CenterPoints , typename RasterizationFunctor /*=std::function< void ( Index ) > )*/ >
#else // !NEW_RASTERIZER
		template< bool Centered , typename RasterizationFunctor /*=std::function< void ( Index ) > )*/ >
#endif // NEW_RASTERIZER
		void Rasterize( Triangle triangle , RasterizationFunctor && F , Range< 2 > cellRange )
		{
			static_assert( std::is_convertible_v< RasterizationFunctor , std::function< void ( Index ) > > , "[ERROR] RasterizationFunctor poorly formed" );

			// Split the 2D range into two 1D ranges
			Range< 1 > cellRanges[2];
			for( unsigned int d=0 ; d<2 ; d++ ) cellRanges[d].first[0] = cellRange.first[d] , cellRanges[d].second[0] = cellRange.second[d];

			// Find the index of the middle vertex
			unsigned int i1 = -1;
			for( unsigned int d=0 ; d<=2 ; d++ )
				if( ( triangle[d][1]<=triangle[(d+1)%3][1] && triangle[d][1]>=triangle[(d+2)%3][1] ) || ( triangle[d][1]>=triangle[(d+1)%3][1] && triangle[d][1]<=triangle[(d+2)%3][1] ) )
					i1 = d;
			if( i1==-1 ) MK_ERROR_OUT( "Could not find middle vertex: " , triangle );
			unsigned int i0 = (i1+2)%3 , i2 = (i1+1)%3;

			double x1 = triangle[i1][0] , y = triangle[i1][1];

#ifdef NEW_RASTERIZER
			// All three vertices at the same height
			if     ( y==triangle[i0][1] && y==triangle[i2][1] ) return;
			// First and second vertices at the same height
			else if( y==triangle[i0][1]                       ) _Rasterize< CenterPoints >( y , x1 , triangle[i0][0] , triangle[i2] , std::forward< RasterizationFunctor >( F ) , cellRanges );
			// Second and third vertices at the same height
			else if(                       y==triangle[i2][1] ) _Rasterize< CenterPoints >( y , x1 , triangle[i2][0] , triangle[i0] , std::forward< RasterizationFunctor >( F ) , cellRanges );
			// All vertices at different heights
			else
			{
				// Solve for s s.t.:
				//	triangle[i0][1]*(1-s) + triangle[i2][1]*s = triangle[i1][1]
				//	(triangle[i2][1]-triangle[i0][1]) * s = triangle[i1][1] - triangle[i0][1]
				//	s = (triangle[i1][1] - triangle[i0][1]) / (triangle[i2][1]-triangle[i0][1])
				double s = (triangle[i1][1] - triangle[i0][1]) / (triangle[i2][1]-triangle[i0][1]);
				double x2 = triangle[i0][0]*(1-s) + triangle[i2][0]*s;
				_Rasterize< CenterPoints >( y , x1 , x2 , triangle[i0] , std::forward< RasterizationFunctor >( F ) , cellRanges );
				_Rasterize< CenterPoints >( y , x1 , x2 , triangle[i2] , std::forward< RasterizationFunctor >( F ) , cellRanges );
			}
#else // !NEW_RASTERIZER
			// All three vertices at the same height
			if     ( y==triangle[i0][1] && y==triangle[i2][1] ) return;
			// First and second vertices at the same height
			else if( y==triangle[i0][1]                       ) _Rasterize< Centered >( y , x1 , triangle[i0][0] , triangle[i2] , std::forward< RasterizationFunctor >( F ) , cellRanges );
			// Second and third vertices at the same height
			else if(                       y==triangle[i2][1] ) _Rasterize< Centered >( y , x1 , triangle[i2][0] , triangle[i0] , std::forward< RasterizationFunctor >( F ) , cellRanges );
			// All vertices at different heights
			else
			{
				// Solve for s s.t.:
				//	triangle[i0][1]*(1-s) + triangle[i2][1]*s = triangle[i1][1]
				//	(triangle[i2][1]-triangle[i0][1]) * s = triangle[i1][1] - triangle[i0][1]
				//	s = (triangle[i1][1] - triangle[i0][1]) / (triangle[i2][1]-triangle[i0][1])
				double s = (triangle[i1][1] - triangle[i0][1]) / (triangle[i2][1]-triangle[i0][1]);
				double x2 = triangle[i0][0]*(1-s) + triangle[i2][0]*s;
				_Rasterize< Centered >( y , x1 , x2 , triangle[i0] , std::forward< RasterizationFunctor >( F ) , cellRanges );
				_Rasterize< Centered >( y , x1 , x2 , triangle[i2] , std::forward< RasterizationFunctor >( F ) , cellRanges );
			}
#endif // NEW_RASTERIZER
		}
	};


	namespace _Rasterizer2D
	{
		using Index = typename RegularGrid< 2 >::Index;
		template< unsigned int Dim > using Range = typename RegularGrid< Dim >::Range;
		using Triangle = Simplex< double , 2 , 2 >;

		// Invokes the rasterization functor for all cells centers that are (left-bottom) inside of the triangle
		template< typename RasterizationFunctor /*=std::function< void ( Index ) > )*/ >
		void Rasterize( Triangle triangle , RasterizationFunctor && F , Range< 2 > cellRange );

		//////////////////////////////////////////

		Range< 1 > GetCellRange( double s1 , double s2 )
		{
			Range< 1 > range;

			range.first[0] = (int)std::floor( s1+0 );
			range.second[0] = (int)std::ceil( s2-1 );

			return range;
		}

		// Rasterizes a horizontal line segment connected to a  point
		template< typename RasterizationFunctor /*=std::function< void ( Index ) > )*/ >
		void _Rasterize( double y , double x0 , double x1 , Point< double , 2 > tip , RasterizationFunctor &&  F , const Range< 1 > cellRanges[2] )
		{
			double minY = y , maxY = tip[1];
			if( x0>x1 ) std::swap( x0 , x1 );
			if( maxY<minY ) std::swap( maxY , minY );


			auto Intersection = [&]( double _y )
				{
					// Solve for s s.t.:
					//	y*(1-s) + tip[1]*s = _y
					//  (tip[1]-y)*s = _y - y
					//	s = (_y-y) / (tip[1]-y)
					double s = (_y-y) / (tip[1]-y);
					return std::pair< double , double >( x0*(1-s) + tip[0]*s , x1*(1-s) + tip[0]*s );
				};

			auto CellRange = [&]( int iy )
				{
					std::pair< double , double > r1 = Intersection( std::max< double >( minY , iy+0. ) );
					std::pair< double , double > r2 = Intersection( std::min< double >( maxY , iy+1. ) );
					return GetCellRange( std::min< double >( r1.first , r2.first ) , std::max< double >( r1.second , r2.second ) );
				};


			// Get the range of heights containing the triangle
			Range< 1 > iyRange = Range< 1 >::Intersect( cellRanges[1] , GetCellRange( minY , maxY ) );
			for( int iy=iyRange.first[0] ; iy<=iyRange.second[0] ; iy++ )
			{
				Range< 1 > ixRange = Range< 1 >::Intersect( cellRanges[0] , CellRange( iy ) );
				for( int ix=ixRange.first[0] ; ix<=ixRange.second[0] ; ix++ ) F( Index(ix,iy) );
			}
		}

		template< bool Centered , typename RasterizationFunctor /*=std::function< void ( Index ) > )*/ >
		void Rasterize( Triangle triangle , RasterizationFunctor && F , Range< 2 > cellRange )
		{
			static_assert( std::is_convertible_v< RasterizationFunctor , std::function< void ( Index ) > > , "[ERROR] RasterizationFunctor poorly formed" );

			// Split the 2D range into two 1D ranges
			Range< 1 > cellRanges[2];
			for( unsigned int d=0 ; d<2 ; d++ ) cellRanges[d].first[0] = cellRange.first[d] , cellRanges[d].second[0] = cellRange.second[d];

			// Find the index of the middle vertex
			unsigned int i1 = -1;
			for( unsigned int d=0 ; d<=2 ; d++ )
				if( ( triangle[d][1]<=triangle[(d+1)%3][1] && triangle[d][1]>=triangle[(d+2)%3][1] ) || ( triangle[d][1]>=triangle[(d+1)%3][1] && triangle[d][1]<=triangle[(d+2)%3][1] ) )
					i1 = d;
			if( i1==-1 ) MK_ERROR_OUT( "Could not find middle vertex: " , triangle );
			unsigned int i0 = (i1+2)%3 , i2 = (i1+1)%3;

			double x1 = triangle[i1][0] , y = triangle[i1][1];

			// All three vertices at the same height
			if     ( y==triangle[i0][1] && y==triangle[i2][1] ) return;
			// First and second vertices at the same height
			else if( y==triangle[i0][1]                       ) _Rasterize( y , x1 , triangle[i0][0] , triangle[i2] , std::forward< RasterizationFunctor >( F ) , cellRanges );
			// Second and third vertices at the same height
			else if(                       y==triangle[i2][1] ) _Rasterize( y , x1 , triangle[i2][0] , triangle[i0] , std::forward< RasterizationFunctor >( F ) , cellRanges );
			// All vertices at different heights
			else
			{
				// Solve for s s.t.:
				//	triangle[i0][1]*(1-s) + triangle[i2][1]*s = triangle[i1][1]
				//	(triangle[i2][1]-triangle[i0][1]) * s = triangle[i1][1] - triangle[i0][1]
				//	s = (triangle[i1][1] - triangle[i0][1]) / (triangle[i2][1]-triangle[i0][1])
				double s = (triangle[i1][1] - triangle[i0][1]) / (triangle[i2][1]-triangle[i0][1]);
				double x2 = triangle[i0][0]*(1-s) + triangle[i2][0]*s;
				_Rasterize( y , x1 , x2 , triangle[i0] , std::forward< RasterizationFunctor >( F ) , cellRanges );
				_Rasterize( y , x1 , x2 , triangle[i2] , std::forward< RasterizationFunctor >( F ) , cellRanges );
			}
		}
	};
}
#endif // RASTERIZER_2D_INCLUDED