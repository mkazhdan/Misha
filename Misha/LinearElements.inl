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
SquareMatrix< SquareMatrix< double , K > , K+1 > LinearElements< K >::_Elements< EType >::stiffness( void )
{
	SquareMatrix< SquareMatrix< double , K > , K+1 > S;
	for( unsigned int k=0 ; k<=K ; k++ ) for( unsigned int _k=0 ; _k<=K ; _k++ ) for( unsigned int l=0 ; l<K ; l++ ) for( unsigned int _l=0 ; _l<K ; _l++ )
		S(k,_k)(l,_l) = ( _elements[k].d(l) * _elements[_k].d(_l) ).integrateUnitRightSimplex();
	return S;
}

template< unsigned int K >
template< typename EType >
Point< Point< double , K > , K+1 , double > LinearElements< K >::_Elements< EType >::differential( void )
{
	Point< Point< double , K > , K+1 , double > d;
	for( unsigned int k=0 ; k<=K ; k++ ) for( unsigned int _k=0 ; _k<K ; _k++ ) d[k][_k] = _elements[k].d(_k)( Point< double , K >() );
	return d;
}



//////////////////////////
// LinearElements::_Hat //
//////////////////////////
template< unsigned int K >
Polynomial::Polynomial< K , 1 , double > LinearElements< K >::_Hat::Element( unsigned int k )
{
	Point< double , K > p[K+1];
	for( unsigned int k=0 ; k<K ; k++ ) p[k+1][k] = 1.;

	double v[K+1] = {0.};
	v[k] = 1;
	return Interpolant( p , v );
}

//////////////////////////////////////
// LinearElements::_CrouzeixRaviart //
//////////////////////////////////////
template< unsigned int K >
Polynomial::Polynomial< K , 1 , double > LinearElements< K >::_CrouzeixRaviart::Element( unsigned int k )
{
	auto Corner = [&]( unsigned int k )
		{
			if( !k ) return Point< double , K >();
			else
			{
				Point< double , K > p;
				p[k-1] = 1;
				return p;
			}
		};

	Point< double , K > p[K+1];
	for( unsigned int k=0 ; k<=K ; k++ )
	{
		for( unsigned int _k=0 ; _k<=K ; _k++ ) if( k!=_k ) p[k] += Corner( _k );
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
template< bool Hat >
LinearElements< K >::Mesh< Hat >::Mesh( const std::vector< SimplexIndex< K > > & simplices , size_t vNum ) : _simplices( simplices ) , _vNum(vNum)
{
	if constexpr( !Hat )
	{
		size_t idx[K];
		for( unsigned int i=0 ; i<_simplices.size() ; i++ )
		{
			for( unsigned int k=0 ; k<=K ; k++ )
			{
				for( unsigned int j=0 , _k=0 ; _k<=K ; _k++ ) if( k!=_k ) idx[j++] = simplices[i][_k];
				FaceIndex fi( idx );
				auto iter = _simplexFaceMap.find( fi );
				if( iter==_simplexFaceMap.end() ) _simplexFaceMap[ fi ] = _simplexFaceMap.size();
			}
		}
	}
}

template< unsigned int K >
template< bool Hat >
std::optional< typename LinearElements< K >::Mesh< Hat >::ElementIndex > LinearElements< K >::Mesh< Hat >::operator()( size_t s , unsigned int k ) const
{
	if constexpr( Hat ) return _simplices[s][k];
	else
	{
		size_t idx[K];
		for( unsigned int _k=0 , j=0 ; _k<=K ; _k++ ) if( _k!=k ) idx[j++] = _simplices[s][_k];
		auto iter = _simplexFaceMap.find( FaceIndex( idx ) );
		if( iter!=_simplexFaceMap.end() ) return ElementIndex( iter->second , FaceIndex::Sign( idx ) );
		else return std::nullopt;
	}
}

template< unsigned int K >
template< bool Hat >
size_t LinearElements< K >::Mesh< Hat >::_index( size_t s , unsigned int k ) const
{
	if constexpr( Hat ) return _simplices[s][k];
	else
	{
		size_t idx[K];
		for( unsigned int _k=0 , j=0 ; _k<=K ; _k++ ) if( _k!=k ) idx[j++] = _simplices[s][_k];
		auto iter = _simplexFaceMap.find( FaceIndex( idx ) );
		if( iter!=_simplexFaceMap.end() ) return iter->second;
		else MK_THROW( "Could not find element: { " , s , " , " , k , " } -> " , FaceIndex(idx) );
	}
}


template< unsigned int K >
template< bool Hat >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< Hat >::mass( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );
	
	std::conditional_t< Hat , HatElements , CrouzeixRaviartElements > elements;
	SquareMatrix< double , K+1 > _M = elements.mass();

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

	Eigen::SparseMatrix< double > M;
	if constexpr( Hat ) M.resize( _vNum , _vNum );
	else                M.resize( _simplexFaceMap.size() , _simplexFaceMap.size() );
	M.setFromTriplets( triplets.begin() , triplets.end() );
	return M;
}

template< unsigned int K >
template< bool Hat >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< Hat >::stiffness( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );

	std::conditional_t< Hat , HatElements , CrouzeixRaviartElements > elements;
	SquareMatrix< SquareMatrix< double , K > , K+1 > _S = elements.stiffness();

	std::vector< Eigen::Triplet< double , size_t > > triplets( _simplices.size() * ( K+1 ) * ( K+1 ) );

	ThreadPool::ParallelFor
	(
		0 , _simplices.size() ,
		[&]( size_t s )
		{
			size_t off = s * ( K+1 ) * ( K+1 );
			SquareMatrix< double , K > g = metric(s);
#ifdef NEW_CODE
			SquareMatrix< double , K > _g = g.adjugate() / sqrt( std::max< double >( 0 , g.determinant() ) );
#else // !NEW_CODE
			SquareMatrix< double , K > gInv = g.inverse();
			double scale = sqrt( std::max< double >( 0 , g.determinant() ) );
#endif // NEW_CODE
			size_t indices[K+1];
			for( unsigned int k=0 ; k<=K ; k++ ) indices[k] = _index( s , k );
#ifdef NEW_CODE
			for( unsigned int i=0 ; i<=K ; i++ ) for( unsigned int j=0 ; j<=K ; j++ ) triplets[off++] = Eigen::Triplet< double , size_t >( indices[i] , indices[j] , SquareMatrix< double , K+1 >::Dot( _S(i,j) , _g ) );
#else // !NEW_CODE
			for( unsigned int i=0 ; i<=K ; i++ ) for( unsigned int j=0 ; j<=K ; j++ ) triplets[off++] = Eigen::Triplet< double , size_t >( indices[i] , indices[j] , SquareMatrix< double , K+1 >::Dot( _S(i,j) , gInv ) * scale );
#endif // NEW_CODE
		}
	);

	Eigen::SparseMatrix< double > S;
	if constexpr( Hat ) S.resize( _vNum , _vNum );
	else                S.resize( _simplexFaceMap.size() , _simplexFaceMap.size() );
	S.setFromTriplets( triplets.begin() , triplets.end() );
	return S;
}

template< unsigned int K >
template< bool Hat >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< Hat >::differentialMass( MetricFunctor && metric ) const
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
#ifdef NEW_CODE
			SquareMatrix< double , K > _g = g.adjugate() / sqrt( std::max< double >( 0 , g.determinant() ) );
#if 0
{
	// Char poly:
	//		X(g) = det( g - \lambda * Id )
	//		     = ( g(0,0) - \lambda ) * ( g(1,1) - \lambda ) - g(0,1) * g(1,0)
	//		     = \lambda^2 - \lambda * tr(g) + det(g)
	double a = g.determinant() , b = - g.trace() , c = 1.;
	double disc = b*b - 4. * a * c;
	if( disc<0 ) std::cout << "Negative discriminant: " << disc << std::endl;
	else
	{
		double x0 = ( - b - sqrt( disc ) ) / ( 2. * a );
		if( x0<0 ) std::cout << "Negative eigenvalue: " << x0 << std::endl;
	}
}
#endif
			for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ ) triplets[off++] = Eigen::Triplet< double , size_t >( s*K+i , s*K+j , _g(i,j) / 2. );
#else // !NEW_CODE
			SquareMatrix< double , K > gInv = g.inverse();
			double scale = sqrt( std::max< double >( 0 , g.determinant() ) );
			for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ ) triplets[off++] = Eigen::Triplet< double , size_t >( s*K+i , s*K+j , gInv(i,j) * scale / 2. );
#endif // NEW_CODE
		}
	);

	Eigen::SparseMatrix< double > M( _simplices.size()*K , _simplices.size()*K );
	M.setFromTriplets( triplets.begin() , triplets.end() );
	return M;
}

template< unsigned int K >
template< bool Hat >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< Hat >::differential( void ) const
{
	std::conditional_t< Hat , HatElements , CrouzeixRaviartElements > elements;
	Point< Point< double , K > , K+1 , double > _D = elements.differential();

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

	Eigen::SparseMatrix< double > D;
	if constexpr( Hat ) D.resize( _simplices.size()*K , _vNum );
	else                D.resize( _simplices.size()*K , _simplexFaceMap.size() );
	D.setFromTriplets( triplets.begin() , triplets.end() );
	return D;
}

template< unsigned int K >
template< bool Hat >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< Hat >::differentialJ( MetricFunctor && metric ) const
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
template< bool Hat >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
void LinearElements< K >::Mesh< Hat >::sanityCheck( MetricFunctor && metric ) const
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