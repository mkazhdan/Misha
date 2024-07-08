#ifndef RASTERIZER_2D_INCLUDED
#define RASTERIZER_2D_INCLUDED

#include "RegularGrid.h"
#include "Geometry.h"
#include "Miscellany.h"

namespace Rasterizer2D
{
	using Index = typename RegularGrid< 2 >::Index;
	using Range = typename RegularGrid< 2 >::Range;
	using Triangle = Simplex< double , 2 , 2 >;

	// [NOTE] Cell centers are at the half-integers
	// Invokes the rasterization functor for all cells that are (left-bottom) inside of the triangle
	template< typename RasterizationFunctor /*=std::function< void ( Index ) > )*/ >
	void Rasterize( Triangle triangle , RasterizationFunctor F , Range cellRange );

	//////////////////////////////////////////
	template< unsigned int Dim > using _Range = typename RegularGrid< Dim >::Range;

	_Range< 1 > GetCellRange( double s1 , double s2 )
	{
		_Range< 1 > range;

		// Solve for the smallest integer mn s.t.:
		//	mn+0.5 >= s1
		//	mn >= s1-0.5
		range.first[0] = (int)std::ceil( s1-0.5 );

		// Solve for the largest integer mx s.t.:
		//	mx+0.5 < s2
		//	mx < s2-0.5
		range.second[0] = (int)std::floor( s2-0.5 );
		if( (s2-0.5)==range.second[0] ) range.second[0]--;

		return range;
	}

	template< typename RasterizationFunctor /*=std::function< void ( Index ) > )*/ >
	void _Rasterize( double y , double x0 , double x1 , Point< double , 2 > tip , RasterizationFunctor F , const _Range< 1 > cellRanges[2] )
	{
		if( x0>x1 ) std::swap( x0 , x1 );

		auto Intersection = [&]( double _y )
			{
				// Solve for s s.t.:
				//	y*(1-s) + tip[1]*s = _y
				//  (tip[1]-y)*s = _y - y
				//	s = (_y-y) / (tip[1]-y)
				double s = (_y-y) / (tip[1]-y);
				return std::pair< double , double >( x0*(1-s) + tip[0]*s , x1*(1-s) + tip[0]*s );
			};

		if( y<tip[1] )
		{
			_Range< 1 > iyRange = _Range< 1 >::Intersect( cellRanges[1] , GetCellRange( y , tip[1] ) );
			for( int iy=iyRange.first[0] ; iy<=iyRange.second[0] ; iy++ )
			{
				std::pair< double , double > xRange = Intersection( iy );
				_Range< 1 > ixRange = _Range< 1 >::Intersect( cellRanges[0] , GetCellRange( xRange.first , xRange.second ) );
				for( int ix=ixRange.first[0] ; ix<=ixRange.second[0] ; ix++ ) F( Index(ix,iy) );
			}
		}
		else if( y>tip[1] )
		{
			_Range< 1 > iyRange = _Range< 1 >::Intersect( cellRanges[1] , GetCellRange( tip[1] , y ) );
			for( int iy=iyRange.first[0] ; iy<=iyRange.second[0] ; iy++ )
			{
				std::pair< double , double > xRange = Intersection( iy );
				_Range< 1 > ixRange = _Range< 1 >::Intersect( cellRanges[0] , GetCellRange( xRange.first , xRange.second ) );
				for( int ix=ixRange.first[0] ; ix<=ixRange.second[0] ; ix++ ) F( Index(ix,iy) );
			}
		}
	}

	template< typename RasterizationFunctor /*=std::function< void ( Index ) > )*/ >
	void Rasterize( Triangle triangle , RasterizationFunctor F , Range cellRange )
	{
		_Range< 1 > cellRanges[2];
		for( unsigned int d=0 ; d<2 ; d++ ) cellRanges[d].first[0] = cellRange.first[d] , cellRanges[d].second[0] = cellRange.second[d];
		unsigned int i1 = -1;
		for( unsigned int d=0 ; d<=2 ; d++ )
			if( ( triangle[d][1]<=triangle[(d+1)%3][1] && triangle[d][1]>=triangle[(d+2)%3][1] ) || ( triangle[d][1]>=triangle[(d+1)%3][1] && triangle[d][1]<=triangle[(d+2)%3][1] ) )
				i1 = d;
		if( i1==-1 ) ERROR_OUT( "Could not find middle index: " , triangle );
		unsigned int i0 = (i1+2)%3 , i2 = (i1+1)%3;

		double x1 = triangle[i1][0] , y = triangle[i1][1];

		if     ( y==triangle[i0][1] && y==triangle[i2][1] ) return;
		else if( y==triangle[i0][1]                       ) _Rasterize( y , x1 , triangle[i0][0] , triangle[i2] , F , cellRanges );
		else if(                       y==triangle[i2][1] ) _Rasterize( y , x1 , triangle[i2][0] , triangle[i0] , F , cellRanges );
		else
		{
			// Solve for s s.t.:
			//	triangle[i0][1]*(1-s) + triangle[i2][1]*s = triangle[i1][1]
			//	(triangle[i2][1]-triangle[i0][1]) * s = triangle[i1][1] - triangle[i0][1]
			//	s = (triangle[i1][1] - triangle[i0][1]) / (triangle[i2][1]-triangle[i0][1])
			double s = (triangle[i1][1] - triangle[i0][1]) / (triangle[i2][1]-triangle[i0][1]);
			double x2 = triangle[i0][0]*(1-s) + triangle[i2][0]*s;
			_Rasterize( y , x1 , x2 , triangle[i0] , F , cellRanges );
			_Rasterize( y , x1 , x2 , triangle[i2] , F , cellRanges );
		}
	}
};
#endif // RASTERIZER_2D_INCLUDED