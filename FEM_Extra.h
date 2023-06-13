#ifndef FEM_EXTRA_INCLUDED
#define FEM_EXTRA_INCLUDED

#include <string.h>
#include "SparseMatrix.h"
#include "Geometry.h"
#include "Array.h"


namespace FEM_Extra
{
	enum
	{
		QUADRATURE_ANGULAR=1 ,
		QUADRATURE_SQUARE_LENGTH=2 ,
	};

	template< class Real > Point3D< Real > TraceWeights( const SquareMatrix< Real , 2 >& tensor , const Point2D< Real > directions[3] );
	template< class Real > void TraceWeights( const Point2D< Real >* directions , int dirCount , Real* weights );
	template< class Real > void CircularQuadratureWeights( const SquareMatrix< Real , 2 >& tensor , const Point2D< Real >* dirs , int dirCount , Real* weights , int quadratureType );

	// The matrix that takes the values at the three prescribed directions and returns the 2x2 linear transform that best matches
	template< class Real > Matrix< Real , 6 , 4 > LinearFit( const Point2D< Real > directions[3] );
	// The matrix that evaluates the best-fit matrix at the three prescribed directions
	template< class Real > SquareMatrix< Real , 6 > LinearFitEvaluation( const Point2D< Real > directions[3] );
	// The matrix that evaluates the difference between the best-fit predicted values and the actual values
	template< class Real > SquareMatrix< Real , 6 > LinearFitResidual( const Point2D< Real > directions[3] );

	template< class Real > SquareMatrix< Real , 6 > MCTraceForm( const SquareMatrix< Real , 2 >& tensor , const Point2D< Real > directions[3] , int quadratureType=0 );
	template< class Real > SquareMatrix< Real , 6 > TraceForm( const SquareMatrix< Real , 2 >& tensor , const Point2D< Real > directions[3] );


	template< class Real >
	struct RightTriangle : FEM::RightTriangle< Real >
	{
		enum
		{
			DUAL_BARYCENTRIC ,
			DUAL_CIRCUMCENTRIC ,
			DUAL_CIRCUMCENTER_PROJECTED_BARYCENTRIC ,
			DUAL_INCENTRIC ,
			DUAL_ISOGONIC ,
			DUAL_ISOGON_PROJECTED_BARYCENTRIC ,
			DUAL_COUNT
		};
		static const char* DualNames[];
		static const int DualCenters[];
	};

	// This structure represents a Riemmanian mesh, with the triangles giving the connectivity and the square (symmetric) matrices giving the metric
	template< class Real >
	struct RiemannianMesh : public FEM::RiemannianMesh< Real >
	{
		void setTriangleDerivativeDirections( int t , ConstPointer( FEM::EdgeXForm< Real > ) edges , int dualType , Point2D< Real > dirs[3] ) const;

		SparseMatrix< Real , int > gradientDualMatrix( int gradType ) const;

		SparseMatrix< Real , int > vectorFieldStiffnessMatrix( int dualType , int quadratureType=0 ) const;

		// Construct the connection
		// [NOTE] If linear fit is disabled, vectorFieldStiffnessMatrix_ gives the same matrix as vectorFieldStiffnessMatrix (though with additional zero-valued entries)
		SparseMatrix< Real , int > vectorFieldStiffnessMatrix ( ConstPointer( FEM::EdgeXForm< Real > ) edges , int dualType , int quadratureType=0 ) const;
		SparseMatrix< Real , int > vectorFieldStiffnessMatrix_( ConstPointer( FEM::EdgeXForm< Real > ) edges , int dualType , int quadratureType=0 , bool linearFit=true ) const;

		SparseMatrix< Real , int > vectorFieldStiffnessMatrix( ConstPointer( FEM::EdgeXForm< Real > ) edges , ConstPointer( Point2D< Real > ) centers ) const;
		SparseMatrix< Real , int > vectorFieldCovariantDerivativeTraceMatrix( int dualType ) const;
		SparseMatrix< Real , int > vectorFieldCovariantDerivativeTraceMatrix( ConstPointer( FEM::EdgeXForm< Real > ) edges , int dualType ) const;
		SparseMatrix< Real , int > vectorFieldCovariantDerivativeTraceMatrix2( int dualType ) const;
		SparseMatrix< Real , int > vectorFieldCovariantDerivativeTraceMatrix2( ConstPointer( FEM::EdgeXForm< Real > ) edges , int dualType ) const;
		SparseMatrix< Real , int > vectorFieldDivergenceMatrix( void ) const;
		SparseMatrix< Real , int > vectorFieldDivergenceMatrix( ConstPointer( FEM::EdgeXForm< Real > ) edges ) const;

		///////////////////////////
		// Topological Operators //
		///////////////////////////

		// Average per-triangle values into the vertices
		template< class Data > Pointer( Data ) getProlongation( ConstPointer( Data ) faceData                                                          ) const;
		template< class Data > void            setProlongation( ConstPointer( Data ) faceData , Pointer( Data ) vertexData , int systemFlag=SYSTEM_SET ) const;
	};
}
template< class Real > const char* FEM_Extra::RightTriangle< Real >::DualNames[] = { "barycentric" , "circumcentric" , "circumcenter projected barycentric" , "incentric" , "isogonic" , "isogon projected barycentric" };

template< class Real > const int FEM_Extra::RightTriangle< Real >::DualCenters[] =
{
	FEM::RightTriangle< Real >::CENTER_BARYCENTRIC ,
	FEM::RightTriangle< Real >::CENTER_CIRCUMCENTRIC ,
	FEM::RightTriangle< Real >::CENTER_BARYCENTRIC ,
	FEM::RightTriangle< Real >::CENTER_INCENTRIC ,
	FEM::RightTriangle< Real >::CENTER_ISOGONIC ,
	FEM::RightTriangle< Real >::CENTER_BARYCENTRIC
};

#include "FEM_Extra.inl"

#endif // FEM_EXTRA_INCLUDED