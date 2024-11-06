#ifndef NORMAL_SMOOTHER
#define NORMAL_SMOOTHER

#include "Eigen/Sparse"
#include "SimplexMesh.h"

namespace NormalSmoother
{
	template< unsigned int K , unsigned int Degree=1 >
	void Smooth( const std::vector< Point< double , K+1 > > &vertices , std::vector< Point< double , K+1 > > &normals , const std::vector< SimplexIndex< K > > &simplices , unsigned int iters , double timeStep )
	{
		static const unsigned int Dim = K+1;

		using Solver = Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >;

		struct TangentFrame
		{
			Point< double , Dim > t[K];

			const Point< double , Dim > &operator[]( unsigned int idx ) const { return t[idx]; }

			void set( Point< double , Dim > n )
			{
				Point< double , Dim > frame[Dim];
				{
					frame[0] = n;

					for( unsigned int d=1 ; d<Dim-1 ; d++ )
					{
						frame[d] = RandomSpherePoint< double , Dim >();
						while( true )
						{
							for( unsigned int dd=0 ; dd<d ; dd++ ) frame[d] -= Point< double , Dim >::Dot( frame[dd] , frame[d] ) * frame[dd];
							if( frame[d].squareNorm()>1e-10 )
							{
								frame[d] /= sqrt( frame[d].squareNorm() );
								break;
							}
						}
					}
					frame[Dim-1] = Point< double , Dim >::CrossProduct( frame );
				}
				for( unsigned int k=0 ; k<K ; k++ ) t[k] = frame[k+1];
			}
		};

		auto DirectSum = []( const Eigen::SparseMatrix< double > &M )
			{
				Eigen::SparseMatrix< double > _M;
				_M.resize( M.rows() * Dim , M.cols() * Dim );
				Eigen::VectorXi sizes( M.rows() * Dim );

#pragma omp parallel for
				for( int k=0; k<M.outerSize() ; k++ )
				{
					unsigned int rowSize = 0;
					for( Eigen::SparseMatrix<double>::InnerIterator it(M,k) ; it ; ++it ) rowSize++;
					for( unsigned int d=0 ; d<Dim ; d++ ) sizes[k*Dim+d] = rowSize;
				}
				_M.reserve( sizes );

#pragma omp parallel for
				for( int k=0; k<M.outerSize() ; k++ )
					for( Eigen::SparseMatrix<double>::InnerIterator it(M,k) ; it ; ++it )
						for( unsigned int d=0 ; d<Dim ; d++ ) _M.insert( it.row()*Dim+d , it.col()*Dim+d ) = it.value();

				return _M;
			};

		auto Prolongation = []( const std::vector< TangentFrame > &frames )
			{
				Eigen::SparseMatrix< double > _P;
				_P.resize( frames.size()*Dim , frames.size()*K );
				_P.reserve( Eigen::VectorXi::Constant( frames.size()*K , Dim ) );
#pragma omp parallel for
				for( int f=0 ; f<(int)frames.size() ; f++ )
					for( unsigned int k=0 ; k<K ; k++ ) for( unsigned int d=0 ; d<Dim ; d++ ) _P.insert( f*Dim+d , f*K+k ) = frames[f][k][d];
				return _P;
			};



		Eigen::SparseMatrix< double > A , M , S;
		{
			SimplexMesh< K , Degree > sMesh = SimplexMesh< K , Degree >::template Init< Dim , unsigned int >( simplices , [&]( unsigned int v ){ return vertices[v]; } );
			sMesh.makeUnitVolume();
			Eigen::SparseMatrix< double > _M = sMesh.mass();
			Eigen::SparseMatrix< double > _S = sMesh.stiffness();
			M = DirectSum( _M );
			S = DirectSum( _S );
		}
		A = M + S * timeStep;

		std::vector< TangentFrame > tangents( vertices.size() );
		Solver *solver = nullptr;
		for( unsigned int iter=0 ; iter<iters ; iter++ )
		{
			Miscellany::Timer timer;

			// Set the tangent directions
#pragma omp parallel for
			for( int i=0 ; i<vertices.size() ; i++ ) tangents[i].set( normals[i] );

			Eigen::SparseMatrix< double > P = Prolongation( tangents );
			Eigen::SparseMatrix< double > Pt = P.transpose();


			// E( o ) = M( n + P*o - n ) + t * S( n + P*o )
			//        = o^t * P^t * ( M + t*S ) * P * o + 2 * t * o^t * P^t * S * n + ...
			// =>
			//  P^t * ( M + t*S ) * P * o = - t * P^t * S * n
			if( !solver ) solver = new Solver( Pt * A * P );
			else          solver->factorize( Pt * A * P );

			Eigen::VectorXd n( vertices.size()*Dim );
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) n[i*Dim+d] = normals[i][d];
			n += P * solver->solve( - Pt * S * n * timeStep );

			for( unsigned int i=0 ; i<vertices.size() ; i++ )
			{
				for( unsigned int d=0 ; d<Dim ; d++ ) normals[i][d] = n[i*Dim+d];
				normals[i] /= Point< double , Dim >::Length( normals[i] );
			}
		}
		delete solver;
	}
}
#endif // NORMAL_SMOOTHER