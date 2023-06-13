#include <math.h>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <atomic>
#include "lineqn.h"

#undef SUPPORT_LINEAR_PROGRAM

template< class Real >
inline SquareMatrix< Real , 6 > FEM_Extra::TraceForm( const SquareMatrix< Real , 2 >& tensor , const Point2D< Real > dirs[3] )
{
	SquareMatrix< Real , 6 > tForm;

	// The linear operator that takes the values of the matrix along the three directions
	// and returns the best-fit linear transform.
	Matrix< Real , 6 , 4 > lFit = LinearFit( dirs );

	// Given the vectors {v0,v1,v2} and {w0,w1,w2}, we would like to:
	// -- Compute the best fit linear operators: L_v = L(v0,v1,v2) and L_w = L(w0,w1,w2)
	// -- Compute the bilinear form: B_{v,w}(x,y) = tensor( L_v * x )[ L_w * y ] = x^t * L_v^t * tensor * L_w * y
	// -- Compute the corresponding linear operator: L_{v,w} = tensor^{-1} * B_{v,w} = tensor^{-1} * L_v^t * tensor * L_w
	// -- Take its trace: Tr( L_{v,w} )
	SquareMatrix< Real , 2 > iTensor = tensor.inverse();
	for( int i=0 ; i<6 ; i++ ) for( int j=0 ; j<6 ; j++ )
	{
		SquareMatrix< Real , 2 > L_v , L_w;
		Real *L_v_ptr = &L_v(0,0) , *L_w_ptr = &L_w(0,0);
		for( int k=0 ; k<4 ; k++ ) L_v_ptr[k] = lFit(i,k) , L_w_ptr[k] = lFit(j,k);
		SquareMatrix< Real , 2 > L_vw = iTensor * L_v.transpose() * tensor * L_w;
		tForm(i,j) = L_vw(0,0) + L_vw(1,1);
	}
	return tForm;
}
template< class Real >
inline SquareMatrix< Real , 6 > FEM_Extra::LinearFitEvaluation( const Point2D< Real > dirs[3] )
{
	// The linear operator that takes the values of the matrix along the three directions
	// and returns the best-fit linear transform.
	Matrix< Real , 6 , 4 > lFit = LinearFit( dirs );

	// Given the vectors {v0,v1,v2} , we would like to:
	// -- Compute the best fit linear operators: L_v = L(v0,v1,v2)
	// -- Compute the difference between the predicted and actual values: dv_{j} = ( L_v * directions[j] - v_{j} )
	SquareMatrix< Real , 6 > lfEvaluation;

	for( int i=0 ; i<6 ; i++ )
	{
		SquareMatrix< Real , 2 > L_v;
		Real *L_v_ptr = &L_v(0,0);
		for( int k=0 ; k<4 ; k++ ) L_v_ptr[k] = lFit(i,k);
		for( int j=0 ; j<3 ; j++ )
		{
			Point2D< Real > temp = L_v * dirs[j];
			lfEvaluation( i , j*2+0 ) = temp[0] , lfEvaluation( i , j*2+1 ) = temp[1];
		}
	}

	return lfEvaluation;
}
template< class Real > inline SquareMatrix< Real , 6 > FEM_Extra::LinearFitResidual( const Point2D< Real > dirs[3] ){ return LinearFitEvaluation( dirs ) - SquareMatrix< Real , 6 >::Identity(); }

template< class Real >
inline SquareMatrix< Real , 6 > FEM_Extra::MCTraceForm( const SquareMatrix< Real , 2 >& tensor , const Point2D< Real > dirs[3] , int quadratureType )
{
	SquareMatrix< Real , 6 > tForm;

	Point3D< Real > weights( 1 , 1 , 1 );
	CircularQuadratureWeights( tensor , dirs , 3 , &weights[0] , quadratureType );
	weights /= (Real)( M_PI );

	for( int i=0 ; i<3 ; i++ )
	{
		weights[i] /= Point2D< Real >::Dot( dirs[i] , tensor * dirs[i] );
		for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ ) tForm(2*i+j,2*i+k) = tensor(j,k) * weights[i];
	}
	return tForm;
}

#ifdef SUPPORT_LINEAR_PROGRAM
#include <Eigen/Dense>
#define IL_STD
#include <ilcplex/ilocplex.h>
template< class Real >
void FEM_Extra::TraceWeights( const Point2D< Real >* directions , int dirCount , Real* weights )
{
	{
		IloEnv env;
		IloModel model( env );
		IloNumVarArray var( env , dirCount , 0.0 , IloInfinity , ILOFLOAT );
		IloExpr ojbective( env );
		for( int i=0 ; i<dirCount ; i++ ) ojbective += 1. * var[i];
		model.add( IloMinimize( env , ojbective ) );
		IloExpr constraint1( env ) , constraint2( env ) , constraint3( env );
		for( int i=0 ; i<dirCount ; i++ ) constraint1 += ( directions[i][0] * directions[i][0] ) * var[i] , constraint2 += ( directions[i][1] * directions[i][1] ) * var[i] , constraint3 += ( directions[i][0] * directions[i][1] ) * var[i];
		model.add( constraint1==1. ) , model.add( constraint2==1. ) , model.add( constraint3==0. );

		IloCplex cplex( model );
//		cplex.setParam(IloCplex::Param::Simplex::Display, 2);
//		cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Network);
		cplex.setOut( env.getNullStream() );
		if( !cplex.solve() )
		{
			IloNumArray x( env ) ; cplex.getValues( x , var );
			fprintf( stderr , "[ERROR] Failed to solve linear program\n" );
			std::vector< Real > angles( dirCount );
			for( int i=0 ; i<dirCount ; i++ ) angles[i] = atan2( directions[i][1] , directions[i][0] );
//			std::sort( angles.begin() , angles.end() , []( Real a1 , Real a2 ){ return a1<a2; } );
			fprintf( stderr , "\tAngles: " ) ; for( int i=0 ; i<dirCount ; i++ ) fprintf( stderr , " %g" , angles[i] * 360. / 2. / M_PI ) ; fprintf( stderr , "\n" );
			fprintf( stderr , "\tWeights: " ) ; for( int i=0 ; i<dirCount ; i++ ) fprintf( stderr , " %g" , x[i] ) ; fprintf( stderr , "\n" );
			exit( 0 );
		}

		constraint1.end();
		constraint2.end();
		constraint3.end();
		ojbective.end();
		IloNumArray x( env ) ; cplex.getValues( x , var );
		for( int i=0 ; i<dirCount ; i++ ) weights[i] = x[i];
		env.end();
		return;
	}

	// Want the set of weights such that:
	//		\sum_i w[i] * dir[i] * dir[i]^t = g^{-1}
	// with w[i]>=0.
	// For uniformly distributed directions on the circle, we expect w[i] = 2 / n
	// For uniformly distributed directions, we expect w[i] = 2 / ( n * || dir[i] ||^2 )
	Eigen::MatrixXd M( dirCount+3 , dirCount+3 );
	Eigen::VectorXd b( dirCount+3 ) , x( dirCount+3 );

	SquareMatrix< Real , 2 > I = SquareMatrix< Real , 2 >::Identity();
	static const int Indices[][2] = { { 0 , 0 } , { 0 , 1 } , { 1 , 1 } };

	// Construct the equality constraints
	for( int i=0 ; i<3 ; i++ )
	{
		for( int j=0 ; j<dirCount ; j++ ) M( dirCount+i , j ) = M( j , dirCount+i ) = directions[j][ Indices[i][0] ] * directions[j][ Indices[i][1] ];
		b[ dirCount+i ] = I( Indices[i][0] ,  Indices[i][1] );
	}
	for( int i=0 ; i<dirCount ;  i++ ) b[i] = 0;

	// Add the optimization constraints
	// E = \sum_{i \neq j} ( w[i] - w[j] )^2
	//   = \sum_{i \neq j} w^2[i] + w*2[j] - 2 * w[i] * w[j]
//	for( int i=0 ; i<dirCount ; i++ ) for( int j=0 ; j<dirCount ; j++ ) M(i,j) = -1.;
	for( int i=0 ; i<dirCount ; i++ ) for( int j=0 ; j<dirCount ; j++ ) M(i,j) =  0.;
	for( int i=0 ; i<dirCount ; i++ ) M(i,i) = dirCount;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) M( i+dirCount , j+dirCount ) = 0;

	// Solve the linear sytem
	x = Eigen::FullPivLU< Eigen::MatrixXd >( M ).solve( b );

	if( dirCount<20 )
	{
		Eigen::MatrixXd M( dirCount , dirCount );
		M *= 0;
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<dirCount ; j++ ) M(i,j) = directions[j][ Indices[i][0] ] * directions[j][ Indices[i][1] ];
		Eigen::MatrixXd K = Eigen::FullPivLU< Eigen::MatrixXd >( M ).kernel();
		printf( "Kernel dim: %d x %d\n" , K.cols() , K.rows() );
std::cout << K << std::endl;
	}

	for( int i=0 ; i<dirCount ; i++ ) weights[i] = (Real)x[i];
}
#endif // SUPPORT_LINEAR_PROGRAM
template< class Real >
inline Point3D< Real > FEM_Extra::TraceWeights( const SquareMatrix< Real , 2 >& tensor , const Point2D< Real > directions[3] )
{
	// Given a linear operator L, a tensor g, and directions v[3], we want to find the weights w[3] such that:
	//		Tr( L ) = \sum w[i] * < v[i] , L( v[i] ) >_g
	// Re-writing the RHS we get:
	//		Tr( L ) = \sum w[i] * < v[i] , ( g * L ) * v[i] >
	//		        = \sum w[i] * Tr( v[i]^t * g * L * v[i] )
	//		        = Tr( L * ( \sum w[i] * v[i] * v[i]^t ) * g )
	// Thus, the weights w[3] must satisfy:
	//		\sum w[i] * v[i] * v[i]^t = g^{-1}
	static const int Indices[][2] = { { 0 , 0 } , { 0, 1 } , { 1 , 1 } };
	SquareMatrix< Real , 2 > tensor_inverse = tensor.inverse();
	SquareMatrix< Real , 2 > M[3];
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ ) M[i](j,k) = directions[i][j] * directions[i][k];
	SquareMatrix< Real , 3 > A;
	Point3D< Real > b;
	for( int i=0 ; i<3 ; i++ )
	{
		b[i] = tensor_inverse( Indices[i][0] , Indices[i][1] );
		for( int j=0 ; j<3 ; j++ ) A(i,j) = M[i]( Indices[j][0] , Indices[j][1] );
	}
	return A.inverse() * b;
}
template< class Real >
inline Matrix< Real , 6 , 4 > FEM_Extra::LinearFit( const Point2D< Real > v[3] )
{
// Given directions v[3] and values w[3], solve for the linear operator L minimizing:
//	E(L) = \sum_i ||L(v[i]) - w[i]||_g^2
//       = \sum_i ( L(v[i]) - w[i] )^t * g * ( L(v[i]) - w[i] )
//       = \sum_i v[i]^t * L^t * g * L * v[i] - 2 * v[i]^t * L^t * g*  w[i] + ...
//       = \sum_i Tr( L^t * g * L * v[i] * v[i]^t ) - 2 Tr( L^t * g * w[i] * v[i]^t ) + ...
// Setting V = \sum_i v[i] * v[i]^t and W = \sum_i w[i] * v[i]^t, this gives:
//       = Tr( L^t * g * L * V ) - 2 * Tr( L^t * g * W )
//       = [ Tr( L^t * g * L * V ) + Tr( L * V * L^t * g ) ] / 2 - 2 * Tr( L^t * g * W )
//       = [ Tr( L^t * g * L * V ) + Tr( g * L * V^t * L^t ) ] / 2 - 2 * Tr( L^t * g * W )
//       = [ Tr( L^t * g * L * V ) + Tr( L^T * g * L * V^t ] / 2 - 2 * Tr( L^t * g * W )
// Differentiating with respect to to L gives:
// g * L * ( V + V^t )/2 = g * W
// Or equivalently:
//	    L = W * V^{-1}
	auto OuterProduct = [] ( Point2D< Real > v1 , Point2D< Real > v2 )
	{
		SquareMatrix< Real , 2 > M;
		for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) M(i,j) = v1[i] * v2[j];
		return M;
	};

	// Compute the sum of the outer products of the drections and invert
	SquareMatrix< Real , 2 > V;
	for( int i=0 ; i<3 ; i++ ) V += OuterProduct( v[i] , v[i] );
	V = V.inverse();

	Matrix< Real , 6 , 4 > fitMatrix;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<2 ; j++ )
	{
		SquareMatrix< Real , 2 > W = OuterProduct( v[i] , j==0 ? Point2D< Real >(1,0) : Point2D< Real >(0,1) ) * V;
		memcpy( &fitMatrix(2*i+j,0) , &W(0,0) , sizeof(Real)*4 );
	}
	return fitMatrix;
}
template< class Real > void FEM_Extra::CircularQuadratureWeights( const SquareMatrix< Real , 2 >& tensor , const Point2D< Real >* dirs , int dirCount , Real* weights , int quadratureType )
{
	if( quadratureType & QUADRATURE_ANGULAR )
	{
		struct IndexedAngle { int idx ; Real angle; };
		Point2D< Real > x( 1 , 0 ) , y = Rotate90( tensor , x );
		x = tensor * x , y = tensor * y;
		std::vector< IndexedAngle > iAngles( 2*dirCount );
		for( int i=0 ; i<dirCount ; i++ )
		{
			weights[i] = 0;
			iAngles[2*i].idx = iAngles[2*i+1].idx = i;
			iAngles[2*i].angle = atan2( Point2D< Real >::Dot( y , dirs[i] ) , Point2D< Real >::Dot( x , dirs[i] ) ) , iAngles[2*i+1].angle = iAngles[2*i].angle + (Real)M_PI;
		}
		for( int i=0 ; i<2*dirCount ; i++ )
		{
			while( iAngles[i].angle<0       ) iAngles[i].angle += (Real)( 2. * M_PI );
			while( iAngles[i].angle>2.*M_PI ) iAngles[i].angle -= (Real)( 2. * M_PI );
		}
		std::sort( iAngles.begin() , iAngles.end() , [] ( const IndexedAngle& i1 , const IndexedAngle& i2 ){ return i1.angle<i2.angle; } );
		for( int i=0 ; i<2*dirCount ; i++ )
		{
			Real a1 , a2;
			if( i==0            ) a1 = ( iAngles[i].angle + iAngles[2*dirCount-1].angle - 2. * M_PI ) / 2;
			else                  a1 = ( iAngles[i].angle + iAngles[i-1         ].angle             ) / 2;
			if( i==2*dirCount-1 ) a2 = ( iAngles[i].angle + iAngles[0           ].angle + 2. * M_PI ) / 2;
			else                  a2 = ( iAngles[i].angle + iAngles[i+1         ].angle             ) / 2;
			weights[ iAngles[i].idx ] += a2 - a1;
		}
	}
	else for( int i=0 ; i<dirCount ; i++ ) weights[i] = (Real)( 2. * M_PI  / dirCount );
	if( quadratureType & QUADRATURE_SQUARE_LENGTH )
	{
		Real sum = 0;
		for( int i=0 ; i<dirCount ; i++ )
		{
			Real l = Point2D< Real >::Dot( dirs[i] , tensor * dirs[i] );
			weights[i] *= l , sum += l;
		}
		for( int i=0 ; i<dirCount ; i++ ) weights[i] /= sum;
	}
}

template< class Real >
inline void FEM_Extra::RiemannianMesh< Real >::setTriangleDerivativeDirections( int t , ConstPointer( FEM::EdgeXForm< Real > ) edges , int dualType , Point2D< Real > dirs[3] ) const
{
	static const double SQRT_THREE_QUARTERS = sqrt( 3./4 );
	for( int j=0 ; j<3 ; j++ )
	{
		int o = edges[t*3+j].oppositeEdge , tt = o / 3;
		if( o!=-1 ) dirs[j] = edges[o].xForm( FEM_Extra::RightTriangle< Real >::Center( g[tt] , DualCenters[dualType] ) ) - FEM_Extra::RightTriangle< Real >::Center( g[t] , DualCenters[dualType] );
		else        dirs[j] = RightTriangle< Real >::EdgeReflect( g[t] , j , FEM_Extra::RightTriangle< Real >::Center( g[t] , DualCenters[dualType] ) ) - FEM_Extra::RightTriangle< Real >::Center( g[t] , DualCenters[dualType] );
		if( dualType==RightTriangle< Real >::DUAL_CIRCUMCENTER_PROJECTED_BARYCENTRIC )
		{
			Point2D< Real > dir = Rotate90( g[t] , RightTriangle< Real >::Edges[j] );
			dirs[j] = dir * Point2D< Real >::Dot( dirs[j] , g[t] * dir ) / Point2D< Real >::Dot( dir , g[t] * dir );
		}
		else if( dualType==RightTriangle< Real >::DUAL_ISOGON_PROJECTED_BARYCENTRIC )
		{
			Point2D< Real > dir = RightTriangle< Real >::EdgeMidpoints[j] - Rotate90( g[t] ,  RightTriangle< Real >::Edges[j] ) * (Real)SQRT_THREE_QUARTERS - RightTriangle< Real >::Center( g[t] , RightTriangle< Real >::DUAL_ISOGONIC );
			dirs[j] = dir * Point2D< Real >::Dot( dirs[j] , g[t] * dir ) / Point2D< Real >::Dot( dir , g[t] * dir );
		}
	}
}

template< class Real >
SparseMatrix< Real , int > FEM_Extra::RiemannianMesh< Real >::gradientDualMatrix( int gradType ) const
{
	int _vCount = vCount();
	SparseMatrix< Real , int > grad;
	Point2D< Real > _grads[] = { Point2D< Real >( (Real)-1 , (Real)-1 ) , Point2D< Real >( (Real)1 , (Real)0 ) , Point2D< Real >( (Real)0 , (Real)1 ) };
	grad.resize( tCount*2 );
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		Real a = area(i);
		SquareMatrix< Real , 2 > gInverse = g[i].inverse();
		if( gradType==HAT_GRADIENT_AND_ROTATED_GRADIENT ) grad.SetRowSize( 2*i , 6 ) , grad.SetRowSize( 2*i+1 , 6 );
		else                                              grad.SetRowSize( 2*i , 3 ) , grad.SetRowSize( 2*i+1 , 3 );
		for( int j=0 ; j<3 ; j++ )
		{
			Point2D< Real > _grad = gInverse * _grads[j];
			Point2D< Real > _gradPerp = Rotate90( g[i] , _grad );
			_grad = g[i] * _grad * a;
			_gradPerp = g[i] * _gradPerp * a;
			int inOffset = 0 , outOffset = 0;
			if( gradType&HAT_GRADIENT )
			{
				grad[2*i+0][j+inOffset] = MatrixEntry< Real , int >( triangles[i][j] + outOffset , _grad[0] );
				grad[2*i+1][j+inOffset] = MatrixEntry< Real , int >( triangles[i][j] + outOffset , _grad[1] );
				inOffset = 3 , outOffset = _vCount;
			}
			if( gradType&HAT_ROTATED_GRADIENT )
			{
				grad[2*i+0][j+inOffset] = MatrixEntry< Real , int >( triangles[i][j] + outOffset , _gradPerp[0] );
				grad[2*i+1][j+inOffset] = MatrixEntry< Real , int >( triangles[i][j] + outOffset , _gradPerp[1] );
			}
		}
	}
	return grad.transpose( gradType==HAT_GRADIENT_AND_ROTATED_GRADIENT ? 2*_vCount : _vCount );
}

template< class Real >
SparseMatrix< Real , int > FEM_Extra::RiemannianMesh< Real >::vectorFieldStiffnessMatrix( int dualType , int quadratureType ) const
{
	std::vector< FEM::EdgeXForm< Real > > edges( tCount*3 );
	setEdgeXForms( GetPointer( edges ) );
	return vectorFieldStiffnessMatrix( edges , dualType , quadratureType );
}
template< class Real >
SparseMatrix< Real , int > FEM_Extra::RiemannianMesh< Real >::vectorFieldDivergenceMatrix( void ) const
{
	std::vector< FEM::Mesh::Edge< Real > > edges( tCount*3 );
	setEdgeXForms( GetPointer( edges ) );
	return vectorFieldDivergenceMatrix( edges );
}
template< class Real >
SparseMatrix< Real , int > FEM_Extra::RiemannianMesh< Real >::vectorFieldCovariantDerivativeTraceMatrix( int dualType ) const
{
	std::vector< FEM::Mesh::Edge< Real > > edges( tCount*3 );
	setEdgeXForms( GetPointer( edges ) );
	return vectorFieldCovariantDerivativeTraceMatrix( edges , dualType );
}
template< class Real >
SparseMatrix< Real , int > FEM_Extra::RiemannianMesh< Real >::vectorFieldCovariantDerivativeTraceMatrix2( int dualType ) const
{
	std::vector< FEM::Mesh::Edge< Real > > edges( tCount*3 );
	setEdgeXForms( GetPointer( edges ) );
	return vectorFieldCovariantDerivativeTraceMatrix2( edges , dualType );
}

template< class Real >
SparseMatrix< Real , int > FEM_Extra::RiemannianMesh< Real >::vectorFieldStiffnessMatrix( ConstPointer( FEM::EdgeXForm< Real > ) edges , ConstPointer( Point2D< Real > ) centers ) const
{
	SparseMatrix< Real , int > stiffness;
	stiffness.resize( 2*tCount );
	Pointer( Real ) edgeWeights = AllocPointer< Real >( tCount*3 );
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		Point2D< Real > dirs[3];
		for( int j=0 ; j<3 ; j++ )
		{
			int e = i*3+j , oe = edges[e].oppositeEdge;
			if( oe!=-1 ) dirs[j] = edges[oe].xForm( centers[oe/3] ) - centers[i];
		}
#if 1
		Real a = ( area(i) / 3 * 2 );
		for( int j=0 ; j<3 ; j++ ) edgeWeights[i*3+j] = a / Point2D< Real >::Dot( dirs[j] , g[i] * dirs[j] );
for( int j=0 ; j<3 ; j++ ) if( Point2D< Real >::Dot( dirs[j] , g[i] * dirs[j] )==0 ) printf( "uh oh\n" );
#else
		Real a = area(i);
		Point3D< Real > weights = TraceWeights( g[i] , dirs );
if( weights[0]<0 || weights[1]<0 || weights[2]<0 ) printf( "Weights: %g %g %g\n" , weights[0] , weights[1] , weights[2] );
		for( int j=0 ; j<3 ; j++ ) edgeWeights[3*i+j] = weights[j] * a;
#endif
	}

	for( int i=0 ; i<tCount ; i++ )
	{
		int count = 1;
		for( int j=0 ; j<3 ; j++ ) if( edges[i*3+j].oppositeEdge!=-1 ) count++;
		stiffness.SetRowSize( 2*i , 2*count ) , stiffness.SetRowSize( 2*i+1 , 2*count );
		stiffness[2*i+0][0] = MatrixEntry< Real , int >( 2*i+0 , (Real)0 ) , stiffness[2*i+1][0] = MatrixEntry< Real , int >( 2*i+0 , (Real)0 );
		stiffness[2*i+0][1] = MatrixEntry< Real , int >( 2*i+1 , (Real)0 ) , stiffness[2*i+1][1] = MatrixEntry< Real , int >( 2*i+1 , (Real)0 );

		count = 1;
		for( int j=0 ; j<3 ; j++ ) if( edges[i*3+j].oppositeEdge!=-1 )
		{
			int edge = i*3+j , oppositeEdge = edges[edge].oppositeEdge , ii = oppositeEdge / 3 , jj = oppositeEdge % 3;

			Real s = edgeWeights[edge] + edgeWeights[oppositeEdge];
			stiffness[2*i+0][0].Value += s*g[i](0,0) , stiffness[2*i+1][0].Value += s*g[i](0,1);
			stiffness[2*i+0][1].Value += s*g[i](1,0) , stiffness[2*i+1][1].Value += s*g[i](1,1);

			SquareMatrix< Real , 2 > xPort = g[i] * edges[ oppositeEdge ].xForm.linear;
			stiffness[2*i+0][2*count+0] = MatrixEntry< Real , int >( 2*ii+0 , -xPort(0,0)*s ) , stiffness[2*i+1][2*count+0] = MatrixEntry< Real , int >( 2*ii+0 , -xPort(0,1)*s );
			stiffness[2*i+0][2*count+1] = MatrixEntry< Real , int >( 2*ii+1 , -xPort(1,0)*s ) , stiffness[2*i+1][2*count+1] = MatrixEntry< Real , int >( 2*ii+1 , -xPort(1,1)*s );
			count++;
		}
	}
	FreePointer( edgeWeights );
	return stiffness;
}
template< class Real >
SparseMatrix< Real , int > FEM_Extra::RiemannianMesh< Real >::vectorFieldStiffnessMatrix( ConstPointer( FEM::EdgeXForm< Real > ) edges , int dualType , int quadratureType ) const
{
	SparseMatrix< Real , int > stiffness;
	stiffness.resize( 2*tCount );
	Pointer( Real ) edgeWeights = AllocPointer< Real >( tCount*3 );
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		Real a = area(i);
		Point2D< Real > dirs[3];
		setTriangleDerivativeDirections( i , edges , dualType , dirs );
		Point3D< Real > weights;
		CircularQuadratureWeights( g[i] , dirs , 3 , &weights[0] , quadratureType );
		weights /= (Real)( M_PI );
		for( int j=0 ; j<3 ; j++ ) edgeWeights[i*3+j] = a / Point2D< Real >::Dot( dirs[j] , g[i] * dirs[j] ) * weights[j];
	}

#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		int count = 1;
		for( int j=0 ; j<3 ; j++ ) if( edges[i*3+j].oppositeEdge!=-1 ) count++;
		stiffness.SetRowSize( 2*i , 2*count ) , stiffness.SetRowSize( 2*i+1 , 2*count );

		count = 1;
		for( int k=0 ; k<2 ; k++ ) for( int l=0 ; l<2 ; l++ ) stiffness[2*i+k][l] = MatrixEntry< Real , int >( 2*i+l , (Real)0 );
		for( int j=0 ; j<3 ; j++ ) if( edges[i*3+j].oppositeEdge!=-1 )
		{
			int ii = edges[i*3+j].oppositeEdge / 3;
			for( int k=0 ; k<2 ; k++ ) for( int l=0 ; l<2 ; l++ ) stiffness[2*i+k][2*count+l] = MatrixEntry< Real , int >( 2*ii+l , (Real)0 );
			count++;
		}

		count = 1;
		for( int j=0 ; j<3 ; j++ ) if( edges[i*3+j].oppositeEdge!=-1 )
		{
			int e = i*3+j , oe = edges[e].oppositeEdge;

			Real s = edgeWeights[e] + edgeWeights[oe];

			SquareMatrix< Real , 2 > xPort = g[i] * edges[oe].xForm.linear;
			for( int k=0 ; k<2 ; k++ ) for( int l=0 ; l<2 ; l++ ) stiffness[2*i+l][k].Value += s*g[i](k,l) , stiffness[2*i+l][2*count+k].Value -= xPort(k,l)*s;
			count++;
		}
	}
	FreePointer( edgeWeights );
	return stiffness;
}
template< class Real >
SparseMatrix< Real , int > FEM_Extra::RiemannianMesh< Real >::vectorFieldStiffnessMatrix_( ConstPointer( FEM::EdgeXForm< Real > ) edges , int dualType , int quadratureType , bool linearFit ) const
{
	struct Entry
	{
		int i , j;
		SquareMatrix< Real , 2 > v;
		Entry( int i=-1 , int j=-1 , SquareMatrix< Real , 2 > v = SquareMatrix< Real , 2 >() ){ this->i = i , this->j = j , this->v = v; }
	};
	std::vector< Entry > entries( tCount*16 );

	auto setEntries = []( const RiemannianMesh< Real >* mesh , ConstPointer( FEM::EdgeXForm< Real > ) edges , int t , int dualType , Pointer( Entry ) entries , int quadratureType , bool linearFit )
	{
		entries += 16*t;

		Point2D< Real > dirs[3];
		mesh->setTriangleDerivativeDirections( t , edges , dualType , dirs );

		int tIndices[] = { t , -1 , -1 , -1 };
		Matrix< Real , 8 , 6 > finiteDifference;
		SquareMatrix< Real , 2 > identity = SquareMatrix< Real , 2 >::Identity();
		for( int v=0 ; v<3 ; v++ )
		{
			int e = t*3+v , oe = edges[e].oppositeEdge , ot = oe / 3;
			if( oe!=-1 )
			{
				tIndices[v+1] = ot;
				for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) finiteDifference( i , 2*v+j ) = identity(i,j) , finiteDifference( 2*(v+1)+i , 2*v+j ) = -edges[oe].xForm.linear(i,j);
			}
		}
		SquareMatrix< Real , 8 > form;
		if( linearFit )
		{
			SquareMatrix< Real , 6 > tForm = TraceForm( mesh->g[t] , dirs );
			SquareMatrix< Real , 6 > res = LinearFitResidual( dirs );
			SquareMatrix< Real , 6 > d = res.transpose() * MCTraceForm( mesh->g[t] , dirs , quadratureType ) * res;
			form = finiteDifference.transpose() * ( tForm + d ) * finiteDifference * mesh->area(t);
		}
		else form = finiteDifference.transpose() * MCTraceForm( mesh->g[t] , dirs , quadratureType ) * finiteDifference * mesh->area(t);

		for( int i=0 ; i<4 ; i++ ) for( int j=0 ; j<4 ; j++ )
		{
			entries[4*i+j].i = tIndices[i] , entries[4*i+j].j = tIndices[j];
			// [NOTE] The indexing is reversed because the subsequent setting of the matrix coefficients is done by having "i" index the row, not the column
			for( int ii=0 ; ii<2 ; ii++ ) for( int jj=0 ; jj<2 ; jj++ ) entries[4*i+j].v(ii,jj) = form(2*i+jj,2*j+ii);
		}
	};

#pragma omp parallel for
	for( int t=0 ; t<tCount ; t++ ) setEntries( this , edges , t , dualType , GetPointer( entries ) , quadratureType , linearFit );

	SparseMatrix< Real , int > stiffness;
	stiffness.resize( 2*tCount );
#pragma omp parallel for
	for( int i=0 ; i<entries.size() ; i++ ) if( entries[i].i!=-1 && entries[i].j!=-1 )
	{
#pragma omp atomic
		stiffness.rowSizes[ 2*entries[i].i+0 ] += 2;
#pragma omp atomic
		stiffness.rowSizes[ 2*entries[i].i+1 ] += 2;
	}
#pragma omp parallel for
	for( int i=0 ; i<stiffness.rows ; i++ )
	{
		int temp = stiffness.rowSizes[i];
		stiffness.rowSizes[i] = 0;
		stiffness.SetRowSize( i , temp );
		stiffness.rowSizes[i] = 0;
	}

	for( int i=0 ; i<entries.size() ; i++ ) if( entries[i].i!=-1 && entries[i].j!=-1 )
	{
		int ii = entries[i].i , jj = entries[i].j;
		const SquareMatrix< Real , 2 >& temp = entries[i].v;
		stiffness[ 2*ii+0 ][ stiffness.rowSizes[2*ii+0]++ ] = MatrixEntry< Real , int >( 2*jj+0 , temp(0,0) );
		stiffness[ 2*ii+0 ][ stiffness.rowSizes[2*ii+0]++ ] = MatrixEntry< Real , int >( 2*jj+1 , temp(1,0) );
		stiffness[ 2*ii+1 ][ stiffness.rowSizes[2*ii+1]++ ] = MatrixEntry< Real , int >( 2*jj+0 , temp(0,1) );
		stiffness[ 2*ii+1 ][ stiffness.rowSizes[2*ii+1]++ ] = MatrixEntry< Real , int >( 2*jj+1 , temp(1,1) );
	}

	return stiffness;
}
template< class Real >
SparseMatrix< Real , int > FEM_Extra::RiemannianMesh< Real >::vectorFieldDivergenceMatrix( ConstPointer( FEM::EdgeXForm< Real > ) edges ) const
{
	Point2D< Real > corners[] = { Point2D< Real >( (Real)0 , (Real)0 ) , Point2D< Real >( (Real)1 , (Real)0 ) , Point2D< Real >( (Real)0 , (Real)1 ) };
	SparseMatrix< Real , int > divergence;
	divergence.resize( tCount );

	for( int i=0 ; i<tCount ; i++ )
	{
		Real a = area(i);
		int count = 0;
		for( int j=0 ; j<3 ; j++ ) if( edges[i*3+j].oppositeEdge!=-1 ) count++;
		divergence.SetRowSize( i , 2*count );

		count = 0;
		for( int j=0 ; j<3 ; j++ ) if( edges[i*3+j].oppositeEdge!=-1 )
		{
			int edge = i*3+j , oppositeEdge = edges[edge].oppositeEdge , ii = oppositeEdge / 3 , jj = oppositeEdge % 3;
			Point2D< Real > e = Rotate90( g[i] , corners[ (j+2)%3 ] - corners[ (j+1)%3 ] );
			// The contribution across edge e is:
			//		e^t * g[i] * edges[ oppositeEdge ].xForm.linear
			e = ( ( SquareMatrix< Real , 2 > )edges[oppositeEdge].xForm.linear.transpose() ) * ( g[i] * e );
			e /= a * 2;
			divergence[i][2*count+0] = MatrixEntry< Real , int >( 2*ii + 0 , e[0] );
			divergence[i][2*count+1] = MatrixEntry< Real , int >( 2*ii + 1 , e[1] );
			count++;
		}
	}
	return divergence;
}
template< class Real >
SparseMatrix< Real , int > FEM_Extra::RiemannianMesh< Real >::vectorFieldCovariantDerivativeTraceMatrix( ConstPointer( FEM::EdgeXForm< Real > ) edges , int dualType ) const
{
	SparseMatrix< Real , int > covariantDerivativeTrace;
	covariantDerivativeTrace.resize( tCount );

	Pointer( Real ) triangleAreas = AllocPointer< Real >( tCount );
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ ) triangleAreas[i] = area( i );

#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		int count = 1;
		for( int j=0 ; j<3 ; j++ ) if( edges[i*3+j].oppositeEdge!=-1 ) count++;
		covariantDerivativeTrace.SetRowSize( i , 2*count );
		for( int k=0 ; k<2 ; k++ ) covariantDerivativeTrace[i][k] = MatrixEntry< Real , int >( 2*i+k , (Real)0 );

		Point2D< Real > dirs[3];
		setTriangleDerivativeDirections( i , edges , dualType , dirs );
		Point3D< Real > traceWeights = TraceWeights( g[i] , dirs );

		count = 1;
		for( int j=0 ; j<3 ; j++ ) if( edges[i*3+j].oppositeEdge!=-1 )
		{
			int edge = i*3+j , oppositeEdge = edges[edge].oppositeEdge , ii = oppositeEdge / 3 , jj = oppositeEdge % 3;

			// Given triangles T and T', the covariant derivative across the shared edge will be:
			//		( edges[ oppositeEdge ].xForm.linear * V[T'] - V[T] ) / l
			// And the contribution to the trace will be:
			//		< ( edges[ oppositeEdge ].xForm.linear * V[T'] - V[T] ) / l , dirs[j] >_g * traceWeights[j]
			//		( < V[T'] , edges[ oppositeEdge ].xForm.linear^t * g * dirs[j] > - < V[T] , g * dirs[j] > ) * traceWeights[j] / l
			Point2D< Real > gDir = g[i] * dirs[j] * traceWeights[j];

			for( int k=0 ; k<2 ; k++ ) covariantDerivativeTrace[i][k].Value -= gDir[k];

			gDir = ( ( SquareMatrix< Real , 2 > )edges[ oppositeEdge ].xForm.linear.transpose() ) * gDir;
			for( int k=0 ; k<2 ; k++ ) covariantDerivativeTrace[i][2*count+k] = MatrixEntry< Real , int >( 2*ii+k , gDir[k] );

			count++;
		}
	}
	FreePointer( triangleAreas );
	return covariantDerivativeTrace;
}
template< class Real >
SparseMatrix< Real , int > FEM_Extra::RiemannianMesh< Real >::vectorFieldCovariantDerivativeTraceMatrix2( ConstPointer( FEM::EdgeXForm< Real > ) edges , int dualType ) const
{
	SparseMatrix< Real , int > covariantDerivativeTrace;
	covariantDerivativeTrace.resize( tCount );

#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		// Get the number of triangle neighbors
		int count = 1;
		for( int j=0 ; j<3 ; j++ ) if( edges[i*3+j].oppositeEdge!=-1 ) count++;
		covariantDerivativeTrace.SetRowSize( i , 2*count );

		// Initialize the matrix entries
		count = 1;
		for( int k=0 ; k<2 ; k++ ) covariantDerivativeTrace[i][k] = MatrixEntry< Real , int >( 2*i+k , (Real)0 );
		for( int j=0 ; j<3 ; j++ ) if( edges[i*3+j].oppositeEdge!=-1 )
		{
			int ii = edges[i*3+j].oppositeEdge / 3;
			for( int k=0 ; k<2 ; k++ ) covariantDerivativeTrace[i][2*count+k] = MatrixEntry< Real , int >( 2*ii+k , (Real)0 );
			count++;
		}

		Point2D< Real > dirs[3];
		setTriangleDerivativeDirections( i , edges , dualType , dirs );
		Matrix< Real , 6 , 4 > linearFit = LinearFit( dirs );

		count = 1;
		for( int j=0 ; j<3 ; j++ ) if( edges[i*3+j].oppositeEdge!=-1 )
		{
			// Given triangles T and T', the change in derivative across the shared edge will be:
			//		edges[ oppositeEdge ].xForm.linear * V[T'] - V[T]
			// And the contribution to the covariance matrix will be:
			//		linearFits[j] * ( edges[ oppositeEdge ].xForm.linear * V[T'] - V[T] )
			Matrix< Real , 2 , 4 > lFit , _lFit;
			for( int k=0 ; k<2 ; k++ ) for( int l=0 ; l<4 ; l++ ) lFit(k,l) = linearFit(j*2+k,l);
			_lFit = lFit * edges[ edges[i*3+j].oppositeEdge ].xForm.linear;

			for( int k=0 ; k<2 ; k++ ) covariantDerivativeTrace[i][k].Value -= lFit(k,0) + lFit(k,3) , covariantDerivativeTrace[i][2*count+k].Value += _lFit(k,0) + _lFit(k,3);

			count++;
		}
	}
	return covariantDerivativeTrace;
}
template< class Real >
template< class Data >
inline Pointer( Data ) FEM_Extra::RiemannianMesh< Real >::getProlongation( ConstPointer( Data ) faceData ) const
{
	Pointer( Data ) vertexData = NewPointer< Data >( vCount() );
	setProlongation( faceData , vertexData );
	return vertexData;
}
template< class Real >
template< class Data >
inline void FEM_Extra::RiemannianMesh< Real >::setProlongation( ConstPointer( Data ) faceData , Pointer( Data ) vertexData , int systemFlag ) const
{
	int vCount = this->vCount();
	Pointer( double ) areas = NewPointer< double >( vCount );
	if( !( systemFlag & SYSTEM_ADD ) )
#pragma omp parallel for
		for( int i=0 ; i<vCount ; i++ ) areas[i] = 0 , vertexData[i] *= (Real)0;
	else
#pragma omp parallel for
		for( int i=0 ; i<vCount ; i++ ) areas[i] = 0;

	for( int i=0 ; i<tCount ; i++ )
	{
		double a = area(i) , _a = a;
		if( systemFlag & SYSTEM_NEGATE ) _a = -_a;
		for( int j=0 ; j<3 ; j++ )
		{
			areas[ triangles[i][j] ] += a;
			vertexData[ triangles[i][j] ] += faceData[i] * (Real)_a;
		}
	}
#pragma omp parallel for
	for( int i=0 ; i<vCount ; i++ ) vertexData[i] /= (Real)areas[i];
	DeletePointer( areas );
}

#ifdef SUPPORT_LINEAR_PROGRAM
#undef SUPPORT_LINEAR_PROGRAM
#endif // SUPPORT_LINEAR_PROGRAM