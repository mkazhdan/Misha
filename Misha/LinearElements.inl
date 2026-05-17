/*
Copyright (c) 2026, Michael Kazhdan
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


////////////////////
// LinearElements //
////////////////////
template< unsigned int K >
Polynomial::Polynomial< K , 1 , double > LinearElements< K >::Interpolant( const Point< double , K > p[K+1] , const double v[K+1] )
{
	Polynomial::Polynomial< K , 1 , double > P;

	// The evaluation matrix
	SquareMatrix< double , K+1 > E;
	Point< double , K+1 > V;
	for( unsigned int k=0 ; k<=K ; k++ )
	{
		E(0,k) = 1;
		for( unsigned int _k=0 ; _k<K ; _k++ ) E(_k+1,k) = p[k][_k];
		V[k] = v[k];
	}
	V = E.inverse() * V;
	unsigned int idx[K];
	for( unsigned int k=0 ; k<K ; k++ ) idx[k] = 0;
	P.coefficient( idx ) = V[0];
	for( unsigned int k=0 ; k<K ; k++ )
	{
		idx[k] = 1;
		P.coefficient( idx ) = V[k+1];
		idx[k] = 0;
	}

	return P;
}

template< unsigned int K >
SquareMatrix< double , K > LinearElements< K >::J( const SquareMatrix< double , K > & g ){ return _J< false >( g ); }

template< unsigned int K >
template< bool InverseMetricTensor >
SquareMatrix< double , K > LinearElements< K >::_J( const SquareMatrix< double , K > & g )
{
	static_assert( K==2 , "[ERROR] 90-degree rotation only supported for 2-dimensional geometry" );
	// Recall:
	// Given an inner-product g, the rotation by 90-degrees is given by:
	//			J = g^{-1} * A * sqrt( det( g ) )
	// with A the skew-symmetric matrix with ones on the off diagonal.

	SquareMatrix< double , K > A;
	A(0,1) = 1 , A(1,0) = -1;

	double scale = sqrt( std::max< double >( 0 , g.determinant() ) );
	if constexpr( InverseMetricTensor ) return g * A / scale;
	else                                return g.adjugate() * A / scale;
}

template< unsigned int K >
Point< double , K > LinearElements< K >::Corner( unsigned int k )
{
	if( !k ) return Point< double , K >();
	else
	{
		Point< double , K > p;
		p[k-1] = 1;
		return p;
	}
}

///////////////////////////////
// LinearElements::_Elements //
///////////////////////////////
template< unsigned int K >
template< typename EType >
LinearElements< K >::_Elements< EType >::_Elements( void ){ for( unsigned int k=0 ; k<=K ; k++ ) _elements[k] = EType::Element(k); }

template< unsigned int K >
template< typename EType >
SquareMatrix< double , K+1 > LinearElements< K >::_Elements< EType >::mass( void ) const
{
	SquareMatrix< double , K+1 > M;
	for( unsigned int k=0 ; k<=K ; k++ ) for( unsigned int _k=0 ; _k<=K ; _k++ ) M(k,_k) = ( _elements[k] * _elements[_k] ).integrateUnitRightSimplex();
	return M;
}

template< unsigned int K >
template< typename EType >
SquareMatrix< SquareMatrix< double , K > , K+1 > LinearElements< K >::_Elements< EType >::stiffness( void ) const
{
	SquareMatrix< SquareMatrix< double , K > , K+1 > S;
	for( unsigned int k=0 ; k<=K ; k++ ) for( unsigned int _k=0 ; _k<=K ; _k++ ) for( unsigned int l=0 ; l<K ; l++ ) for( unsigned int _l=0 ; _l<K ; _l++ )
		S(k,_k)(l,_l) = ( _elements[k].d(l) * _elements[_k].d(_l) ).integrateUnitRightSimplex();
	return S;
}

template< unsigned int K >
template< typename EType >
Point< Point< double , K > , K+1 , double > LinearElements< K >::_Elements< EType >::differential( void ) const
{
	Point< Point< double , K > , K+1 , double > d;
	for( unsigned int k=0 ; k<=K ; k++ ) for( unsigned int _k=0 ; _k<K ; _k++ ) d[k][_k] = _elements[k].d(_k)( Point< double , K >() );
	return d;
}

template< unsigned int K >
template< typename EType >
SquareMatrix< double , K+1 > LinearElements< K >::_Elements< EType >::Mass( void )
{
	_Elements e;
	return e.mass();
}

template< unsigned int K >
template< typename EType >
SquareMatrix< SquareMatrix< double , K > , K+1 > LinearElements< K >::_Elements< EType >::Stiffness( void )
{
	_Elements e;
	return e.stiffness();
}

template< unsigned int K >
template< typename EType >
Point< Point< double , K > , K+1 , double > LinearElements< K >::_Elements< EType >::Differential( void )
{
	_Elements e;
	return e.differential();
}


/////////////////////////
// LinearElements::Hat //
/////////////////////////
template< unsigned int K >
Polynomial::Polynomial< K , 1 , double > LinearElements< K >::Hat::Element( unsigned int k )
{
	Point< double , K > p[K+1];
	for( unsigned int k=0 ; k<K ; k++ ) p[k+1][k] = 1.;

	double v[K+1] = {0.};
	v[k] = 1;
	return Interpolant( p , v );
}

/////////////////////////////////////
// LinearElements::CrouzeixRaviart //
/////////////////////////////////////
template< unsigned int K >
Polynomial::Polynomial< K , 1 , double > LinearElements< K >::CrouzeixRaviart::Element( unsigned int k )
{
	static const typename SimplexIndex< K >::template Faces< K-1 > faces;
	Point< double , K > p[K+1];
	for( unsigned int k=0 ; k<=K ; k++ )
	{
		SimplexIndex< K-1 > si = faces[k];
		for( unsigned int _k=0 ; _k<K ; _k++ ) p[k] += Corner( si[_k] );
		p[k] /= K;
	}

	double v[K+1] = {0.};
	v[k] = 1;
	return Interpolant( p , v );
}

//////////////////////////
// LinearElements::Mesh //
//////////////////////////

template< unsigned int K >
template< enum LinearElements< K >::ScalarElementType Type >
LinearElements< K >::Mesh< Type >::Mesh( const std::vector< SimplexIndex< K > > & simplices , size_t vNum ) : SubSimplexMesh< K , SubK() >( simplices ) , _vNum(vNum)
{
}

template< unsigned int K >
template< enum LinearElements< K >::ScalarElementType Type >
size_t LinearElements< K >::Mesh< Type >::elementNum( void ) const
{
	return _SubSimplexMesh::size();
}

template< unsigned int K >
template< enum LinearElements< K >::ScalarElementType Type >
typename LinearElements< K >::Mesh< Type >::ElementSimplex LinearElements< K >::Mesh< Type >::elementSimplex( size_t i ) const
{
	return _SubSimplexMesh::operator[]( i );
}

template< unsigned int K >
template< enum LinearElements< K >::ScalarElementType Type >
typename LinearElements< K >::Mesh< Type >::ElementIndex LinearElements< K >::Mesh< Type >::operator()( size_t s , unsigned int k ) const
{
	if      constexpr( Type==ScalarElementType::HAT              ) return _SubSimplexMesh::operator()( s , k ).first; 
	else if constexpr( Type==ScalarElementType::CROUZEIX_RAVIART ) return _SubSimplexMesh::operator()( s , k );
	else static_assert( false , "[ERROR] Unrecognized Type" );
}

template< unsigned int K >
template< enum LinearElements< K >::ScalarElementType Type >
size_t LinearElements< K >::Mesh< Type >::_index( size_t s , unsigned int k ) const
{
	return _SubSimplexMesh::operator()( s , k ).first;
}


template< unsigned int K >
template< enum LinearElements< K >::ScalarElementType Type >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< Type >::mass( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );
	
	SquareMatrix< double , K+1 > _M = Elements::Mass();

	std::vector< Eigen::Triplet< double , size_t > > triplets( _SubSimplexMesh::_simplices.size() * ( K+1 ) * ( K+1 ) );

	ThreadPool::ParallelFor
		(
			0 , _SubSimplexMesh::_simplices.size() ,
			[&]( size_t s )
			{
				size_t off = s * ( K+1 ) * ( K+1 );
				double scale = sqrt( std::max< double >( 0 , metric( s ).determinant() ) );
				size_t indices[K+1];
				for( unsigned int k=0 ; k<=K ; k++ ) indices[k] = _index( s , k );
				for( unsigned int i=0 ; i<=K ; i++ ) for( unsigned int j=0 ; j<=K ; j++ ) triplets[off++] = Eigen::Triplet< double , size_t >( indices[i] , indices[j] , _M(i,j) * scale );
			}
		);

	Eigen::SparseMatrix< double > M( elementNum() , elementNum() );
	M.setFromTriplets( triplets.begin() , triplets.end() );
	return M;
}

template< unsigned int K >
template< enum LinearElements< K >::ScalarElementType Type >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< Type >::stiffness( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );

	SquareMatrix< SquareMatrix< double , K > , K+1 > _S = Elements::Stiffness();

	std::vector< Eigen::Triplet< double , size_t > > triplets( _SubSimplexMesh::_simplices.size() * ( K+1 ) * ( K+1 ) );

	ThreadPool::ParallelFor
	(
		0 , _SubSimplexMesh::_simplices.size() ,
		[&]( size_t s )
		{
			size_t off = s * ( K+1 ) * ( K+1 );
			SquareMatrix< double , K > g = metric(s);
			SquareMatrix< double , K > _g = g.adjugate() / sqrt( std::max< double >( 0 , g.determinant() ) );
			size_t indices[K+1];
			for( unsigned int k=0 ; k<=K ; k++ ) indices[k] = _index( s , k );
			for( unsigned int i=0 ; i<=K ; i++ ) for( unsigned int j=0 ; j<=K ; j++ ) triplets[off++] = Eigen::Triplet< double , size_t >( indices[i] , indices[j] , SquareMatrix< double , K >::Dot( _S(i,j) , _g ) );
		}
	);

	Eigen::SparseMatrix< double > S( elementNum() , elementNum() );
	S.setFromTriplets( triplets.begin() , triplets.end() );
	return S;
}

template< unsigned int K >
template< enum LinearElements< K >::ScalarElementType Type >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< Type >::differentialMass( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );

	std::vector< Eigen::Triplet< double , size_t > > triplets( _SubSimplexMesh::_simplices.size() * K * K );

	ThreadPool::ParallelFor
	(
		0 , _SubSimplexMesh::_simplices.size() ,
		[&]( size_t s )
		{
			size_t off = s * K * K;
			SquareMatrix< double , K > g = metric(s);
			SquareMatrix< double , K > _g = g.adjugate() / sqrt( std::max< double >( 0 , g.determinant() ) );
			for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ ) triplets[off++] = Eigen::Triplet< double , size_t >( s*K+i , s*K+j , _g(i,j) / 2. );
		}
	);

	Eigen::SparseMatrix< double > M( _SubSimplexMesh::_simplices.size()*K , _SubSimplexMesh::_simplices.size()*K );
	M.setFromTriplets( triplets.begin() , triplets.end() );
	return M;
}

template< unsigned int K >
template< enum LinearElements< K >::ScalarElementType Type >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< Type >::differential( void ) const
{
	Point< Point< double , K > , K+1 , double > _D = Elements::Differential();

	std::vector< Eigen::Triplet< double , size_t > > triplets( _SubSimplexMesh::_simplices.size() * (K+1) * K );

	ThreadPool::ParallelFor
	(
		0 , _SubSimplexMesh::_simplices.size() ,
		[&]( size_t s )
		{
			size_t off = s * (K+1) * K;
			for( unsigned int i=0 ; i<=K ; i++ )
			{
				size_t idx = _index( s , i );
				for( unsigned int j=0 ; j<K ; j++ ) triplets[off++] = Eigen::Triplet< double , size_t >( s*K+j , idx , _D[i][j] );
			}
		}
	);

	Eigen::SparseMatrix< double > D( _SubSimplexMesh::_simplices.size()*K , elementNum() );
	D.setFromTriplets( triplets.begin() , triplets.end() );
	return D;
}

template< unsigned int K >
template< enum LinearElements< K >::ScalarElementType Type >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< Type >::differentialJ( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );

	std::vector< Eigen::Triplet< double , size_t > > triplets( _SubSimplexMesh::_simplices.size() * K * K );

	ThreadPool::ParallelFor
	(
		0 , _SubSimplexMesh::_simplices.size() ,
		[&]( size_t s )
		{
			size_t off = s * K * K;
			SquareMatrix< double , K > _J = LinearElements< K >::template _J< true >( metric(s) );
			for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ ) triplets[off++] = Eigen::Triplet< double , size_t >( s*K+j , s*K+i , _J(i,j) );
		}
	);

	Eigen::SparseMatrix< double > J( _SubSimplexMesh::_simplices.size()*K , _SubSimplexMesh::_simplices.size()*K );
	J.setFromTriplets( triplets.begin() , triplets.end() );
	return J;
}

template< unsigned int K >
template< enum LinearElements< K >::ScalarElementType Type >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< Type >::elementToVertex( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );

	Eigen::SparseMatrix< double > e2v( _vNum , elementNum() );
	if constexpr( Type==ScalarElementType::HAT )
	{
		std::vector< Eigen::Triplet< double , size_t > > triplets( elementNum() );
		ThreadPool::ParallelFor( 0 , elementNum() , [&]( size_t e ){ triplets[e] = Eigen::Triplet< double , size_t >( _SubSimplexMesh::operator[]( e )[0] , e , 1. ); } );
		e2v.setFromTriplets( triplets.begin() , triplets.end() );
	}
	else if constexpr( Type==ScalarElementType::CROUZEIX_RAVIART )
	{
		std::vector< double > vWeights( _vNum , 0 );
		std::vector< Eigen::Triplet< double , size_t > > triplets( _simplices.size() * (K+1) * (K+1) );

		ThreadPool::ParallelFor
			(
				0 , _simplices.size() , 
				[&]( size_t s )
				{
					size_t idx = s * (K+1) * (K+1);
					double weight = sqrt( std::max< double >( 0 , metric( s ).determinant() ) );
					for( unsigned int k=0 ; k<=K ; k++ ) 
					{
						AddAtomic( vWeights[ _simplices[s][k] ] , weight );

						size_t e = operator()( s , k ).first;
						_SubSimplexIndex fIdx = _subSimplexIndices[e];
						// The CR functions are 1 on the face and -1 on the opposite vertex
						for( unsigned int _k=0 ; _k<K ; _k++ ) triplets[idx++] = Eigen::Triplet< double , size_t >( fIdx[_k] , e , weight );
						triplets[idx++] = Eigen::Triplet< double , size_t >( _simplices[s][k] , e , -weight );
					}
				}
			);
		for( unsigned int i=0 ; i<triplets.size() ; i++ ) const_cast< double & >( triplets[i].value() ) /= vWeights[ triplets[i].row() ];
		e2v.setFromTriplets( triplets.begin() , triplets.end() );
	}
	else static_assert( false , "[ERROR] Unrecognized Type" );
	return e2v;
}

template< unsigned int K >
template< enum LinearElements< K >::ScalarElementType Type >
template< typename SampleFunctor /* = std::function< std::pair< size_t , Point< double , K > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< Type >::evaluation( size_t sampleNum , SampleFunctor && sampleFunctor ) const
{
	static_assert( std::is_convertible_v< SampleFunctor , std::function< std::pair< size_t , Point< double , K > > ( size_t ) > > , "[ERROR] SampleFunctor is poorly formed" );

	Elements elements;

	Eigen::SparseMatrix< double > E( sampleNum , elementNum() );

	std::vector< Eigen::Triplet< double , size_t > > triplets( sampleNum * (K+1) );

	const Mesh & mesh = *this;

	ThreadPool::ParallelFor
		(
			0 , sampleNum ,
			[&]( size_t s )
			{
				size_t idx = s * (K+1);
				std::pair< size_t , Point< double , K > > sample = sampleFunctor( s );

				if      constexpr( Type==ScalarElementType::HAT              ) for( unsigned int k=0 ; k<=K ; k++ ) triplets[idx++] = Eigen::Triplet< double , size_t >( s , mesh( sample.first , k )       , elements[k]( sample.second ) );
				else if constexpr( Type==ScalarElementType::CROUZEIX_RAVIART ) for( unsigned int k=0 ; k<=K ; k++ ) triplets[idx++] = Eigen::Triplet< double , size_t >( s , mesh( sample.first , k ).first , elements[k]( sample.second ) );
				else static_assert( false , "[ERROR] Unrecognized Type" );
			}
		);

	E.setFromTriplets( triplets.begin() , triplets.end() );
	return E;
}

template< unsigned int K >
template< enum LinearElements< K >::ScalarElementType Type >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
void LinearElements< K >::Mesh< Type >::sanityCheck( MetricFunctor && metric ) const
{
	Eigen::SparseMatrix< double >  M = mass( metric );
	Eigen::SparseMatrix< double >  S = stiffness( metric );
	Eigen::SparseMatrix< double >  D = differential();
	Eigen::SparseMatrix< double > dM = differentialMass( metric );
	Eigen::SparseMatrix< double > dJ = differentialJ( metric );
	Eigen::VectorXd one( _SubSimplexMesh::size() ) , rand( _SubSimplexMesh::size() );
	Eigen::SparseMatrix< double > Id( _SubSimplexMesh::_simplices.size()*K , _SubSimplexMesh::_simplices.size()*K );
	Id.setIdentity();
	for( unsigned int i=0 ; i<one.size() ; i++ ) one[i] = 1. , rand[i] = Random< double >();

	{
		double e = sqrt( ( S - D.transpose() * dM * D ).squaredNorm() / ( S + D.transpose() * dM * D ).squaredNorm() );
		std::cout << "     Stiffness factorization error: " << e << std::endl;
	}
	{
		double e = sqrt( ( Id + dJ * dJ ).squaredNorm() / ( Id - dJ * dJ ).squaredNorm() );
		std::cout << "                   J squared error: " << e << std::endl;
	}
	{
		double e = sqrt( ( dM - dJ.transpose() * dM * dJ ).squaredNorm() / ( dM + dJ.transpose() * dM * dJ ).squaredNorm() );
		std::cout << "Differential mass invariance error: " << e << std::endl;
	}

	{
		Polynomial::Polynomial< K ,0 , double > P;
		P[0] = 1.;
		double scl = P.integrateUnitRightSimplex();
		double measure = 0;
		for( size_t s=0 ; s<_simplices.size() ; s++ ) measure += scl * sqrt( metric(s).determinant() );
		double _measure = one.dot( M * one );
		double e = sqrt( ( measure - _measure ) * ( measure - _measure ) / ( measure + _measure ) / ( measure + _measure ) );
		std::cout << "                     Measure error: " << e << std::endl;
	}

	{
		std::cout << "                   Stiffness Error: " << rand.dot( S * one ) << std::endl;
	}

	if constexpr( K==2 )
	{
		static const unsigned int Quadrature = 12;

		using Integrator = SimplexIntegrator< K , Quadrature >;

		Eigen::SparseMatrix< double > E , Et , D;
		{
			E = evaluation( _simplices.size()*Quadrature , [&]( size_t s ){ return std::pair< size_t , Point< double , K > >( s/Quadrature , Integrator::Positions[s%Quadrature] ); } );
			Et = E.transpose();
		}
		{
			Polynomial::Polynomial< K ,0 , double > P;
			P[0] = 1.;
			double scl = P.integrateUnitRightSimplex();

			Eigen::VectorXd d( _simplices.size()*Quadrature );
			for( unsigned int s=0 ; s<_simplices.size() ;  s++ )
			{
				double _scl = scl * sqrt( metric(s).determinant() );
				for( unsigned int q=0 ; q<Quadrature ; q++ ) d[s*Quadrature+q] = _scl * Integrator::Weights[q];
			}
			D = d.asDiagonal();
		}
		Eigen::SparseMatrix< double > _M = Et * D * E;
		double e = sqrt( ( M - _M ).squaredNorm() / ( M + _M ).squaredNorm() );
		std::cout << "                        Mass error: " << e << std::endl;
	}
}
