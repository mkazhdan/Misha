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
#ifdef NEW_ELEMENT_MESH
	static const typename SimplexIndex< K >::template Faces< K-1 > faces;
#endif // NEW_ELEMENT_MESH
	Point< double , K > p[K+1];
	for( unsigned int k=0 ; k<=K ; k++ )
	{
#ifdef NEW_ELEMENT_MESH
		SimplexIndex< K-1 > si = faces[k];
#else // !NEW_ELEMENT_MESH
		SimplexIndex< K-1 > si = Boundary( k );
#endif // NEW_ELEMENT_MESH
		for( unsigned int _k=0 ; _k<K ; _k++ ) p[k] += Corner( si[_k] );
		p[k] /= K;
	}

	double v[K+1] = {0.};
	v[k] = 1;
	return Interpolant( p , v );
}

#ifdef NEW_ELEMENT_MESH
#else // !NEW_ELEMENT_MESH
template< unsigned int K >
SimplexIndex< K-1 > LinearElements< K >::CrouzeixRaviart::Boundary( unsigned int k )
{
	SimplexIndex< K-1 > si;
	for( unsigned int _k=0 ; _k<K ; _k++ ) si[_k] = (k+1+_k)%(K+1);
	return si;
}
#endif // NEW_ELEMENT_MESH

#ifdef NEW_ELEMENT_MESH
/////////////////////////////////
// LinearElements::ElementMesh //
/////////////////////////////////
template< unsigned int K >
template< unsigned int SubK >
const typename SimplexIndex< K >::template Faces< SubK > LinearElements< K >::ElementMesh< SubK >::_SimplexFaces;

template< unsigned int K >
template< unsigned int SubK >
LinearElements< K >::ElementMesh< SubK >::ElementMesh( const std::vector< SimplexIndex< K > > & simplices ) : _simplices( simplices )
{
	for( unsigned int i=0 ; i<_simplices.size() ; i++ )
	{
		size_t idx[SubK+1];
		for( unsigned int ii=0 ; ii<SimplexIndex< K >::template FaceNum< SubK >() ; ii++ )
		{
			for( unsigned int k=0 ; k<=SubK ; k++ ) idx[k] = simplices[i][ _SimplexFaces[ii][k] ];
			_SubSimplexIndex fi( idx );
			auto iter = _subSimplexMap.find( fi );
			if( iter==_subSimplexMap.end() ) _subSimplexMap[ fi ] = _subSimplexMap.size();
		}
	}
	_subSimplexIndices.resize( _subSimplexMap.size() );
	for( auto iter : _subSimplexMap ) _subSimplexIndices[ iter.second ] = iter.first;
}

template< unsigned int K >
template< unsigned int SubK >
size_t LinearElements< K >::ElementMesh< SubK >::size( void ) const { return _subSimplexIndices.size(); }

template< unsigned int K >
template< unsigned int SubK >
#if 1 // NEW_CODE
SimplexIndex< SubK > LinearElements< K >::ElementMesh< SubK >::operator[]( size_t i ) const
#else // !NEW_CODE
const SimplexIndex< SubK > & LinearElements< K >::ElementMesh< SubK >::operator[]( size_t i ) const
#endif // NEW_CODE
{
	SimplexIndex< SubK > si;
	for( unsigned int k=0 ; k<=SubK ; k++ ) si[k] = static_cast< unsigned int >( _subSimplexIndices[i][k] );
	return si;
}

template< unsigned int K >
template< unsigned int SubK >
typename LinearElements< K >::ElementMesh< SubK >::SignedIndex LinearElements< K >::ElementMesh< SubK >::operator()( size_t s , unsigned int n ) const
{
	size_t idx[SubK+1];
	for( unsigned int k=0 ; k<=SubK ; k++ ) idx[k] = _simplices[s][ _SimplexFaces[n][k] ];
	auto iter = _subSimplexMap.find( _SubSimplexIndex( idx ) );
	if( iter==_subSimplexMap.end() ) MK_THROW( "Failed to find sub-simplex: " , _SubSimplexIndex(idx) );
	return SignedIndex( iter->second , _SubSimplexIndex::Sign( idx ) );
}
#endif // NEW_ELEMENT_MESH

//////////////////////////
// LinearElements::Mesh //
//////////////////////////

template< unsigned int K >
template< bool UseHat >
#ifdef NEW_ELEMENT_MESH
LinearElements< K >::Mesh< UseHat >::Mesh( const std::vector< SimplexIndex< K > > & simplices , size_t vNum ) : LinearElements< K >::ElementMesh< SubK() >( simplices ) , _vNum(vNum)
#else // !NEW_ELEMENT_MESH
LinearElements< K >::Mesh< UseHat >::Mesh( const std::vector< SimplexIndex< K > > & simplices , size_t vNum ) : _simplices( simplices ) , _vNum(vNum)
#endif // NEW_ELEMENT_MESH
{
#ifdef NEW_ELEMENT_MESH
#else // !NEW_ELEMENT_MESH
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
#endif // NEW_ELEMENT_MESH
}

template< unsigned int K >
template< bool UseHat >
size_t LinearElements< K >::Mesh< UseHat >::elementNum( void ) const
{
#ifdef NEW_ELEMENT_MESH
	return _ElementMesh::size();
#else // !NEW_ELEMENT_MESH
	if constexpr( UseHat ) return _vNum;
	else                   return _simplexFaceMap.size();
#endif // NEW_ELEMENT_MESH
}

template< unsigned int K >
template< bool UseHat >
typename LinearElements< K >::Mesh< UseHat >::ElementSimplex LinearElements< K >::Mesh< UseHat >::elementSimplex( size_t i ) const
{
#ifdef NEW_ELEMENT_MESH
	return _ElementMesh::operator[]( i );
#else // !NEW_ELEMENT_MESH
	if constexpr( UseHat ) return ElementSimplex( i );
	else
	{
//std::cout << i << " / " << _simplexFaces.size() << std::endl;
		SimplexIndex< K-1 > si;
		for( unsigned int k=0 ; k<K ; k++ ) si[k] = static_cast< unsigned int >( _simplexFaces[i][k] );
//std::cout << "Returning: " << _simplexFaces[i] << " -> " << si << std::endl;
		return si;
	}
#endif // NEW_ELEMENT_MESH
}

template< unsigned int K >
template< bool UseHat >
typename LinearElements< K >::Mesh< UseHat >::ElementIndex LinearElements< K >::Mesh< UseHat >::operator()( size_t s , unsigned int k ) const
{
#ifdef NEW_ELEMENT_MESH
	if constexpr( UseHat ) return _ElementMesh::operator()( s , k ).first; 
	else                   return _ElementMesh::operator()( s , k );
#else // !NEW_ELEMENT_MESH
	if constexpr( UseHat ) return _simplices[s][k];
	else
	{
		size_t idx[K];
		SimplexIndex< K-1 > si = CrouzeixRaviart::Boundary(k);
		for( unsigned int _k=0 ; _k<K ; _k++ ) idx[_k] = _simplices[s][ si[_k] ];
		auto iter = _simplexFaceMap.find( _FaceIndex( idx ) );
		if( iter==_simplexFaceMap.end() ) MK_THROW( "Failed to find face: " , _FaceIndex(idx) );
		return ElementIndex( iter->second , _FaceIndex::Sign( idx ) );
	}
#endif // NEW_ELEMENT_MESH
}

template< unsigned int K >
template< bool UseHat >
size_t LinearElements< K >::Mesh< UseHat >::_index( size_t s , unsigned int k ) const
{
#ifdef NEW_ELEMENT_MESH
	return _ElementMesh::operator()( s , k ).first;
#else // !NEW_ELEMENT_MESH
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
#endif // NEW_ELEMENT_MESH
}


template< unsigned int K >
template< bool UseHat >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< UseHat >::mass( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );
	
	SquareMatrix< double , K+1 > _M = std::conditional_t< UseHat , HatElements , CrouzeixRaviartElements >::Mass();

#ifdef NEW_ELEMENT_MESH
	std::vector< Eigen::Triplet< double , size_t > > triplets( _ElementMesh::_simplices.size() * ( K+1 ) * ( K+1 ) );
#else // !NEW_ELEMENT_MESH
	std::vector< Eigen::Triplet< double , size_t > > triplets( _simplices.size() * ( K+1 ) * ( K+1 ) );
#endif // NEW_ELEMENT_MESH

	ThreadPool::ParallelFor
		(
#ifdef NEW_ELEMENT_MESH
			0 , _ElementMesh::_simplices.size() ,
#else // !NEW_ELEMENT_MESH
			0 , _simplices.size() ,
#endif // NEW_ELEMENT_MESH
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

#ifdef NEW_ELEMENT_MESH
	std::vector< Eigen::Triplet< double , size_t > > triplets( _ElementMesh::_simplices.size() * ( K+1 ) * ( K+1 ) );
#else // !NEW_ELEMENT_MESH
	std::vector< Eigen::Triplet< double , size_t > > triplets( _simplices.size() * ( K+1 ) * ( K+1 ) );
#endif // NEW_ELEMENT_MESH

	ThreadPool::ParallelFor
	(
#ifdef NEW_ELEMENT_MESH
		0 , _ElementMesh::_simplices.size() ,
#else // !NEW_ELEMENT_MESH
		0 , _simplices.size() ,
#endif // NEW_ELEMENT_MESH
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

#ifdef NEW_ELEMENT_MESH
	std::vector< Eigen::Triplet< double , size_t > > triplets( _ElementMesh::_simplices.size() * K * K );
#else // !NEW_ELEMENT_MESH
	std::vector< Eigen::Triplet< double , size_t > > triplets( _simplices.size() * K * K );
#endif // NEW_ELEMENT_MESH

	ThreadPool::ParallelFor
	(
#ifdef NEW_ELEMENT_MESH
		0 , _ElementMesh::_simplices.size() ,
#else // !NEW_ELEMENT_MESH
		0 , _simplices.size() ,
#endif // NEW_ELEMENT_MESH
		[&]( size_t s )
		{
			size_t off = s * K * K;
			SquareMatrix< double , K > g = metric(s);
			SquareMatrix< double , K > _g = g.adjugate() / sqrt( std::max< double >( 0 , g.determinant() ) );
			for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ ) triplets[off++] = Eigen::Triplet< double , size_t >( s*K+i , s*K+j , _g(i,j) / 2. );
		}
	);

#ifdef NEW_ELEMENT_MESH
	Eigen::SparseMatrix< double > M( _ElementMesh::_simplices.size()*K , _ElementMesh::_simplices.size()*K );
#else // !NEW_ELEMENT_MESH
	Eigen::SparseMatrix< double > M( _simplices.size()*K , _simplices.size()*K );
#endif // NEW_ELEMENT_MESH
	M.setFromTriplets( triplets.begin() , triplets.end() );
	return M;
}

template< unsigned int K >
template< bool UseHat >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< UseHat >::differential( void ) const
{
	Point< Point< double , K > , K+1 , double > _D = std::conditional_t< UseHat , HatElements , CrouzeixRaviartElements >::Differential();

#ifdef NEW_ELEMENT_MESH
	std::vector< Eigen::Triplet< double , size_t > > triplets( _ElementMesh::_simplices.size() * (K+1) * K );
#else // !NEW_ELEMENT_MESH
	std::vector< Eigen::Triplet< double , size_t > > triplets( _simplices.size() * (K+1) * K );
#endif // NEW_ELEMENT_MESH

	ThreadPool::ParallelFor
	(
#ifdef NEW_ELEMENT_MESH
		0 , _ElementMesh::_simplices.size() ,
#else // !NEW_ELEMENT_MESH
		0 , _simplices.size() ,
#endif // NEW_ELEMENT_MESH
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

#ifdef NEW_ELEMENT_MESH
	Eigen::SparseMatrix< double > D( _ElementMesh::_simplices.size()*K , elementNum() );
#else // !NEW_ELEMENT_MESH
	Eigen::SparseMatrix< double > D( _simplices.size()*K , elementNum() );
#endif // NEW_ELEMENT_MESH
	D.setFromTriplets( triplets.begin() , triplets.end() );
	return D;
}

template< unsigned int K >
template< bool UseHat >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< UseHat >::differentialJ( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );

#ifdef NEW_ELEMENT_MESH
	std::vector< Eigen::Triplet< double , size_t > > triplets( _ElementMesh::_simplices.size() * K * K );
#else // !NEW_ELEMENT_MESH
	std::vector< Eigen::Triplet< double , size_t > > triplets( _simplices.size() * K * K );
#endif // NEW_ELEMENT_MESH

	ThreadPool::ParallelFor
	(
#ifdef NEW_ELEMENT_MESH
		0 , _ElementMesh::_simplices.size() ,
#else // !NEW_ELEMENT_MESH
		0 , _simplices.size() ,
#endif // NEW_ELEMENT_MESH
		[&]( size_t s )
		{
			size_t off = s * K * K;
			SquareMatrix< double , K > _J = LinearElements< K >::template _J< true >( metric(s) );
			for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ ) triplets[off++] = Eigen::Triplet< double , size_t >( s*K+j , s*K+i , _J(i,j) );
		}
	);

#ifdef NEW_ELEMENT_MESH
	Eigen::SparseMatrix< double > J( _ElementMesh::_simplices.size()*K , _ElementMesh::_simplices.size()*K );
#else // !NEW_ELEMENT_MESH
	Eigen::SparseMatrix< double > J( _simplices.size()*K , _simplices.size()*K );
#endif // NEW_ELEMENT_MESH
	J.setFromTriplets( triplets.begin() , triplets.end() );
	return J;
}

template< unsigned int K >
template< bool UseHat >
template< typename MetricFunctor /* = std::function< SquareMatrix< double , K+1 > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< UseHat >::elementToVertex( MetricFunctor && metric ) const
{
	static_assert( std::is_convertible_v< MetricFunctor , std::function< SquareMatrix< double , K+1 > ( size_t ) > > , "[ERROR] MetricFunctor is poorly formed" );

	Eigen::SparseMatrix< double > e2v( _vNum , elementNum() );
	if constexpr( UseHat ) e2v.setIdentity();
	else
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
#ifdef NEW_ELEMENT_MESH
						_FaceIndex fIdx = _subSimplexIndices[e];
#else // !NEW_ELEMENT_MESH
						_FaceIndex fIdx = _simplexFaces[e];
#endif // NEW_ELEMENT_MESH
						// The CR functions are 1 on the face and -1 on the opposite vertex
						for( unsigned int _k=0 ; _k<K ; _k++ ) triplets[idx++] = Eigen::Triplet< double , size_t >( fIdx[_k] , e , weight );
						triplets[idx++] = Eigen::Triplet< double , size_t >( _simplices[s][k] , e , -weight );
					}
				}
			);
		for( unsigned int i=0 ; i<triplets.size() ; i++ ) const_cast< double & >( triplets[i].value() ) /= vWeights[ triplets[i].row() ];
		e2v.setFromTriplets( triplets.begin() , triplets.end() );
	}
	return e2v;
}

template< unsigned int K >
template< bool UseHat >
template< typename SampleFunctor /* = std::function< std::pair< size_t , Point< double , Dim > ( size_t ) > */ >
Eigen::SparseMatrix< double > LinearElements< K >::Mesh< UseHat >::evaluation( size_t sampleNum , SampleFunctor && sampleFunctor ) const
{
	static_assert( std::is_convertible_v< SampleFunctor , std::function< std::pair< size_t , Point< double , K > > ( size_t ) > > , "[ERROR] SampleFunctor is poorly formed" );

	using Elements = std::conditional_t< UseHat , HatElements , CrouzeixRaviartElements >;
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

				if constexpr( UseHat ) for( unsigned int k=0 ; k<=K ; k++ ) triplets[idx++] = Eigen::Triplet< double , size_t >( s , mesh( sample.first , k ) , elements[k]( sample.second ) );
				else                   for( unsigned int k=0 ; k<=K ; k++ ) triplets[idx++] = Eigen::Triplet< double , size_t >( s , mesh( sample.first , k ).first , elements[k]( sample.second ) );
			}
		);

	E.setFromTriplets( triplets.begin() , triplets.end() );
	return E;
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
#ifdef NEW_ELEMENT_MESH
	Eigen::SparseMatrix< double > Id( _ElementMesh::_simplices.size()*K , _ElementMesh::_simplices.size()*K );
#else // !NEW_ELEMENT_MESH
	Eigen::SparseMatrix< double > Id( _simplices.size()*K , _simplices.size()*K );
#endif // NEW_ELEMENT_MESH
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

