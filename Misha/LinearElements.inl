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
	Point< double , K > p[K+1];
	for( unsigned int k=0 ; k<=K ; k++ )
	{
		SimplexIndex< K-1 > si = Boundary( k );
		for( unsigned int _k=0 ; _k<K ; _k++ ) p[k] += Corner( si[_k] );
		p[k] /= K;
	}

	double v[K+1] = {0.};
	v[k] = 1;
	return Interpolant( p , v );
}

template< unsigned int K >
SimplexIndex< K-1 > LinearElements< K >::CrouzeixRaviart::Boundary( unsigned int k )
{
	SimplexIndex< K-1 > si;
	for( unsigned int _k=0 ; _k<K ; _k++ ) si[_k] = (k+1+_k)%(K+1);
	return si;
}

//////////////////////////
// LinearElements::Mesh //
//////////////////////////

template< unsigned int K >
template< bool UseHat >
LinearElements< K >::Mesh< UseHat >::Mesh( const std::vector< SimplexIndex< K > > & simplices , size_t vNum ) : _simplices( simplices ) , _vNum(vNum)
{
	if constexpr( !UseHat )
	{
		size_t idx[K];
		for( unsigned int i=0 ; i<_simplices.size() ; i++ )
		{
			for( unsigned int k=0 ; k<=K ; k++ )
			{
				SimplexIndex< K-1 > si = CrouzeixRaviart::Boundary( k );
				for( unsigned int _k=0 ; _k<K ; _k++ ) idx[_k] = simplices[i][ si[_k] ];
				_FaceIndex fi( idx );
				auto iter = _simplexFaceMap.find( fi );
				if( iter==_simplexFaceMap.end() ) _simplexFaceMap[ fi ] = _simplexFaceMap.size();
			}
		}
		_simplexFaces.resize( _simplexFaceMap.size() );
		for( auto iter : _simplexFaceMap ) _simplexFaces[ iter.second ] = iter.first;
	}
}

template< unsigned int K >
template< bool UseHat >
size_t LinearElements< K >::Mesh< UseHat >::elementNum( void ) const
{
	if constexpr( UseHat ) return _vNum;
	else                   return _simplexFaceMap.size();
}

template< unsigned int K >
template< bool UseHat >
typename LinearElements< K >::Mesh< UseHat >::ElementSimplex LinearElements< K >::Mesh< UseHat >::elementSimplex( size_t i ) const
{
	if constexpr( UseHat ) return ElementSimplex( i );
	else
	{
		SimplexIndex< K-1 > si;
		for( unsigned int k=0 ; k<K ; k++ ) si[k] = _simplexFaces[i][k];
		return si;
	}
}

template< unsigned int K >
template< bool UseHat >
std::optional< typename LinearElements< K >::Mesh< UseHat >::ElementIndex > LinearElements< K >::Mesh< UseHat >::operator()( size_t s , unsigned int k ) const
{
	if constexpr( UseHat ) return _simplices[s][k];
	else
	{
		size_t idx[K];
		SimplexIndex< K-1 > si = CrouzeixRaviart::Boundary(k);
		for( unsigned int _k=0 ; _k<K ; _k++ ) idx[_k] = _simplices[s][ si[_k] ];
		auto iter = _simplexFaceMap.find( _FaceIndex( idx ) );
		if( iter!=_simplexFaceMap.end() ) return ElementIndex( iter->second , _FaceIndex::Sign( idx ) );
		else return std::nullopt;
	}
}

template< unsigned int K >
template< bool UseHat >
size_t LinearElements< K >::Mesh< UseHat >::_index( size_t s , unsigned int k ) const
{
	if constexpr( UseHat ) return _simplices[s][k];
	else
	{
		size_t idx[K];
		SimplexIndex< K-1 > si = CrouzeixRaviart::Boundary(k);
		for( unsigned int _k=0 ; _k<K ; _k++ ) idx[_k] = _simplices[s][ si[_k] ];
		auto iter = _simplexFaceMap.find( _FaceIndex( idx ) );
		if( iter!=_simplexFaceMap.end() ) return iter->second;
		else MK_THROW( "Could not find element: { " , s , " , " , k , " } -> " , _FaceIndex(idx) );
	}
}


template< unsigned int K >
template< bool UseHat >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< UseHat >::mass( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );
	
	SquareMatrix< double , K+1 > _M = std::conditional_t< UseHat , HatElements , CrouzeixRaviartElements >::Mass();

	std::vector< Eigen::Triplet< double , size_t > > triplets( _simplices.size() * ( K+1 ) * ( K+1 ) );

	ThreadPool::ParallelFor
		(
			0 , _simplices.size() ,
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
template< bool UseHat >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< UseHat >::stiffness( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );

	SquareMatrix< SquareMatrix< double , K > , K+1 > _S = std::conditional_t< UseHat , HatElements , CrouzeixRaviartElements >::Stiffness();

	std::vector< Eigen::Triplet< double , size_t > > triplets( _simplices.size() * ( K+1 ) * ( K+1 ) );

	ThreadPool::ParallelFor
	(
		0 , _simplices.size() ,
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
template< bool UseHat >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< UseHat >::differentialMass( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );

	std::vector< Eigen::Triplet< double , size_t > > triplets( _simplices.size() * K * K );

	ThreadPool::ParallelFor
	(
		0 , _simplices.size() ,
		[&]( size_t s )
		{
			size_t off = s * K * K;
			SquareMatrix< double , K > g = metric(s);
			SquareMatrix< double , K > _g = g.adjugate() / sqrt( std::max< double >( 0 , g.determinant() ) );
			for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ ) triplets[off++] = Eigen::Triplet< double , size_t >( s*K+i , s*K+j , _g(i,j) / 2. );
		}
	);

	Eigen::SparseMatrix< double > M( _simplices.size()*K , _simplices.size()*K );
	M.setFromTriplets( triplets.begin() , triplets.end() );
	return M;
}

template< unsigned int K >
template< bool UseHat >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< UseHat >::differential( void ) const
{
	Point< Point< double , K > , K+1 , double > _D = std::conditional_t< UseHat , HatElements , CrouzeixRaviartElements >::Differential();

	std::vector< Eigen::Triplet< double , size_t > > triplets( _simplices.size() * (K+1) * K );

	ThreadPool::ParallelFor
	(
		0 , _simplices.size() ,
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

	Eigen::SparseMatrix< double > D( _simplices.size()*K , elementNum() );
	D.setFromTriplets( triplets.begin() , triplets.end() );
	return D;
}

template< unsigned int K >
template< bool UseHat >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< UseHat >::differentialJ( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );

	std::vector< Eigen::Triplet< double , size_t > > triplets( _simplices.size() * K * K );

	ThreadPool::ParallelFor
	(
		0 , _simplices.size() ,
		[&]( size_t s )
		{
			size_t off = s * K * K;
			SquareMatrix< double , K > _J = LinearElements< K >::template _J< true >( metric(s) );
			for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ ) triplets[off++] = Eigen::Triplet< double , size_t >( s*K+j , s*K+i , _J(i,j) );
		}
	);

	Eigen::SparseMatrix< double > J( _simplices.size()*K , _simplices.size()*K );
	J.setFromTriplets( triplets.begin() , triplets.end() );
	return J;
}


template< unsigned int K >
template< bool UseHat >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
void LinearElements< K >::Mesh< UseHat >::sanityCheck( MetricFunctor && metric ) const
{
	Eigen::SparseMatrix< double >  S = stiffness( metric );
	Eigen::SparseMatrix< double >  D = differential();
	Eigen::SparseMatrix< double > dM = differentialMass( metric );
	Eigen::SparseMatrix< double > dJ = differentialJ( metric );
	Eigen::SparseMatrix< double > Id( _simplices.size()*K , _simplices.size()*K );
	Id.setIdentity();

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
}