/*
Copyright (c) 2025, Michael Kazhdan and Matthew Bolitho
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

#if 1 // NEW_CODE
///////////////////
// JacobiRelaxer //
///////////////////

JacobiRelaxer::JacobiRelaxer( const Eigen::SparseMatrix< double > &M ) : _M(M)
{
	_D.resize( _M.rows() );
	for( unsigned int i=0 ; i<_M.outerSize() ; i++ ) for( Eigen::InnerIterator it(_M,i) ; it ; ++it )
		if( it.row()==it.col() ) _D[ it.row() ] = it.value();
}

void JacobiRelaxer::operator()( const Eigen::VectorXd &b , Eigen::VectorXd &x ) const
{
	// \sum_j M_ij * x_j = b_i   ->   x_i <- ( \sum_j M_ij * x_j - b_i - D_i*x_i ) / D_i
	Eigen::VectorXd _x = x;
	ThreadPool::ParallelFor
		(
			0 , _M.outerSize() ,
			[&]( size_t i )
			{
				double value = 0;
				for( Eigen::InnerIterator it(_M,i) ; it ; ++it ) value += it.value() * x[ it.row() ];
				_x[i] += ( b[i]  - value ) / _D[i];
			}
		);
	x = _x;
}

void JacobiRelaxer::operator()( const Eigen::MatrixXd &B , Eigen::MatrixXd &X ) const
{
	// \sum_j M_ij * x_j = b_i   ->   x_i <- ( \sum_j M_ij * x_j - b_i - D_i*x_i ) / D_i
	unsigned int cols = static_cast< unsigned int >( B.cols() );
	Eigen::MatrixXd _X = X;
	ThreadPool::ParallelFor
		(
			0 , _M.outerSize() ,
			[&]( size_t i )
			{
				Eigen::VectorXd values = Eigen::VectorXd::Zero( cols );
				for( Eigen::InnerIterator it(_M,i) ; it ; ++it ) for( unsigned int c=0 ; c<cols ; c++ ) values[c] += it.value() * X( it.row() , c );
				for( unsigned int c=0 ; c<cols ; c++ ) _X(i,c) += ( B(i,c) - values[c] ) / _D[i];
			}
		);
	X = _X;
}

#endif // NEW_CODE

///////////////
// GSRelaxer //
///////////////

GSRelaxer::GSRelaxer( const Eigen::SparseMatrix< double > &M ) : _M(M)
{
	_D.resize( _M.rows() );
	for( unsigned int i=0 ; i<_M.outerSize() ; i++ ) for( Eigen::InnerIterator it(_M,i) ; it ; ++it )
		if( it.row()==it.col() ) _D[ it.row() ] = it.value();
}

void GSRelaxer::operator()( const Eigen::VectorXd &b , Eigen::VectorXd &x ) const
{
	// \sum_j M_ij * x_j = b_i   ->   x_i <- ( \sum_j M_ij * x_j - b_i - D_i*x_i ) / D_i
	for( unsigned int i=0 ; i<_M.outerSize() ; i++ )
	{
		double value = 0;
		for( Eigen::InnerIterator it(_M,i) ; it ; ++it ) value += it.value() * x[ it.row() ];
		x[i] += ( b[i]  - value ) / _D[i];
	}
}

void GSRelaxer::operator()( const Eigen::VectorXd &b , Eigen::VectorXd &x , const std::vector< std::vector< unsigned int > > &mcIndices ) const
{
	// \sum_j M_ij * x_j = b_i   ->   x_i <- ( \sum_j M_ij * x_j - b_i - D_i*x_i ) / D_i
	for( unsigned int i=0 ; i<mcIndices.size() ; i++ )
		ThreadPool::ParallelFor
		(
			0 , mcIndices[i].size() ,
			[&]( size_t j )
			{
				unsigned int idx = mcIndices[i][j];
				double value = 0;
				for( Eigen::InnerIterator it(_M,idx) ; it ; ++it ) value += it.value() * x[ it.row() ];
				x[idx] += ( b[idx]  - value ) / _D[idx];
			}
		);
}

#if 1 // NEW_CODE
void GSRelaxer::operator()( const Eigen::MatrixXd &B , Eigen::MatrixXd &X ) const
{
	// \sum_j M_ij * x_j = b_i   ->   x_i <- ( \sum_j M_ij * x_j - b_i - D_i*x_i ) / D_i
	unsigned int cols = static_cast< unsigned int >( B.cols() );
	for( unsigned int i=0 ; i<_M.outerSize() ; i++ )
	{
		Eigen::VectorXd values = Eigen::VectorXd::Zero( cols );
		for( Eigen::InnerIterator it(_M,i) ; it ; ++it ) for( unsigned int c=0 ; c<cols ; c++ ) values[c] += it.value() * X( it.row() , c );
		for( unsigned int c=0 ; c<cols ; c++ ) X(i,c) += ( B(i,c) - values[c] ) / _D[i];
	}
}

void GSRelaxer::operator()( const Eigen::MatrixXd &B , Eigen::MatrixXd &X , const std::vector< std::vector< unsigned int > > &mcIndices ) const
{
	// \sum_j M_ij * x_j = b_i   ->   x_i <- ( \sum_j M_ij * x_j - b_i - D_i*x_i ) / D_i
	unsigned int cols = static_cast< unsigned int >( B.cols() );
	for( unsigned int i=0 ; i<mcIndices.size() ; i++ )
		ThreadPool::ParallelFor
		(
			0 , mcIndices[i].size() ,
			[&]( size_t j )
			{
				unsigned int idx = mcIndices[i][j];
				Eigen::VectorXd values = Eigen::VectorXd::Zero( cols );
				for( Eigen::InnerIterator it(_M,idx) ; it ; ++it ) for( unsigned int c=0 ; c<cols ; c++ ) values[c] += it.value() * X( it.row() , c );
				for( unsigned int c=0 ; c<cols ; c++ ) X(idx,c) += ( B(idx,c)  - values[c] ) / _D[idx];
			}
		);
}
#endif // NEW_CODE

//////////////
// MGSolver //
//////////////

#if 1 // NEW_CODE
template< typename Solver , bool Jacobi >
template< typename ProlongationFunctor /* = std::function< Eigen::SparseMatrix< double > ( unsigned int ) > */ >
MGSolver< Solver , Jacobi >::MGSolver( const Eigen::SparseMatrix< double > &M , ProlongationFunctor Ps , size_t pNum ) : _M(M)
#else // !NEW_CODE
template< typename Solver >
template< typename ProlongationFunctor /* = std::function< Eigen::SparseMatrix< double > ( unsigned int ) > */ >
MGSolver< Solver >::MGSolver( const Eigen::SparseMatrix< double > &M , ProlongationFunctor Ps , size_t pNum ) : _M(M)
#endif // NEW_CODE
{
	if( pNum )
	{
		_Ms.resize( pNum );
		_Ps.resize( pNum );
		_Rs.resize( pNum );
#if 1 // NEW_CODE
#else // !NEW_CODE
		_Xs.resize( pNum );
		_Bs.resize( pNum );
#endif // NEW_CODE
		for( unsigned int i=0 ; i<pNum ; i++ )
		{
			_Ps[i] = Ps( i );
			_Rs[i] = _Ps[i].transpose();
		}

		_Ms.back() = _Rs[pNum-1] * _M * _Ps[pNum-1];
		for( int i=(int)pNum-2 ; i>=0 ; i-- ) _Ms[i] = _Rs[i] * _Ms[i+1] * _Ps[i];
#if 1 // NEW_CODE
		for( int i=1 ; i<_Ms.size() ; i++ ) _relaxers.emplace_back( _Ms[i] );
		_relaxers.emplace_back( _M );
#else // !NEW_CODE
		for( int i=1 ; i<_Ms.size() ; i++ ) _gsRelaxers.emplace_back( _Ms[i] );
		_gsRelaxers.emplace_back( _M );
		for( unsigned int i=0 ; i<_Ms.size() ; i++ ) _Xs[i].resize( _Ms[i].cols() );
#endif // NEW_CODE

		_solver = new Solver( _Ms[0] );
	}
	else _solver = new Solver( _M );
}

#if 1 // NEW_CODE
template< typename Solver , bool Jacobi >
MGSolver< Solver , Jacobi >::~MGSolver( void ){ delete _solver; }
#else // !NEW_CODE
template< typename Solver >
MGSolver< Solver >::~MGSolver( void ){ delete _solver; }
#endif // NEW_CODE

#if 1 // NEW_CODE
template< typename Solver , bool Jacobi >
void MGSolver< Solver , Jacobi >::vCycle( const Eigen::VectorXd & b , Eigen::VectorXd & x , unsigned int iters ){ return _vCycle( b , x , iters ); }

template< typename Solver , bool Jacobi >
void MGSolver< Solver , Jacobi >::vCycle( const Eigen::MatrixXd & B , Eigen::MatrixXd & X , unsigned int iters ){ return _vCycle( B , X , iters ); }

template< typename Solver , bool Jacobi >
template< typename MCIndicesFunctor /* = std::function< const std::vector< std::vector< unsigned int > > & ( unsigned int ) > */ >
void MGSolver< Solver , Jacobi >::vCycle( const Eigen::VectorXd &b , Eigen::VectorXd &x , unsigned int iters , MCIndicesFunctor MCI ){ return _vCycle( b , x , iters , MCI ); }

template< typename Solver , bool Jacobi >
template< typename MCIndicesFunctor /* = std::function< const std::vector< std::vector< unsigned int > > & ( unsigned int ) > */ >
void MGSolver< Solver , Jacobi >::vCycle( const Eigen::MatrixXd & B , Eigen::MatrixXd & X , unsigned int iters , MCIndicesFunctor MCI ){ return _vCycle( B , X , iters , MCI ); }
#else // !NEW_CODE
template< typename Solver >
void MGSolver< Solver >::vCycle( const Eigen::VectorXd &b , Eigen::VectorXd &x , unsigned int iters )
{
	// Fine-to-coarse
	for( unsigned int l=(unsigned int)_Ps.size() ; l>0 ; l-- )
	{
		if( l==_Ps.size() )
		{
			for( unsigned int i=0 ; i<iters ; i++ ) _gsRelaxers[l-1]( b , x );
			_Bs[l-1] = _Rs[l-1] * ( b - _M * x );
			for( unsigned int i=0 ; i<_Xs[l-1].size() ; i++ ) _Xs[l-1][i] = 0;
		}
		else
		{
			for( unsigned int i=0 ; i<iters ; i++ ) _gsRelaxers[l-1]( _Bs[l] , _Xs[l] );
			_Bs[l-1] = _Rs[l-1] * ( _Bs[l] - _Ms[l] * _Xs[l] );
			for( unsigned int i=0 ; i<_Xs[l-1].size() ; i++ ) _Xs[l-1][i] = 0;
		}
	}
	// Low-res solve
	if( _Ps.size() ) _Xs[0] = _solver->solve( _Bs[0] );
	else             x = _solver->solve( b );

	// Coarse-to-fine
	for( unsigned int l=1 ; l<=_Ps.size() ; l++ )
	{
		if( l==_Ps.size() )
		{
			x += _Ps[l-1] * _Xs[l-1];
			for( unsigned int i=0 ; i<iters ; i++ ) _gsRelaxers[l-1]( b , x );
		}
		else
		{
			_Xs[l] += _Ps[l-1] * _Xs[l-1];
			for( unsigned int i=0 ; i<iters ; i++ ) _gsRelaxers[l-1]( _Bs[l] , _Xs[l] );
		}
	}
}

template< typename Solver >
template< typename MCIndicesFunctor /* = std::function< const std::vector< std::vector< unsigned int > > & ( unsigned int ) > */ >
void MGSolver< Solver >::vCycle( const Eigen::VectorXd &b , Eigen::VectorXd &x , unsigned int iters , MCIndicesFunctor MCI )
{
	// Fine-to-coarse
	for( unsigned int l=(unsigned int)_Ps.size() ; l>0 ; l-- )
	{
		if( l==_Ps.size() )
		{
			for( unsigned int i=0 ; i<iters ; i++ ) _gsRelaxers[l-1]( b , x , MCI(l) );
			_Bs[l-1] = _Rs[l-1] * ( b - _M * x );
			for( unsigned int i=0 ; i<_Xs[l-1].size() ; i++ ) _Xs[l-1][i] = 0;
		}
		else
		{
			for( unsigned int i=0 ; i<iters ; i++ ) _gsRelaxers[l-1]( _Bs[l] , _Xs[l] , MCI(l) );
			_Bs[l-1] = _Rs[l-1] * ( _Bs[l] - _Ms[l] * _Xs[l] );
			for( unsigned int i=0 ; i<_Xs[l-1].size() ; i++ ) _Xs[l-1][i] = 0;
		}
	}
	// Low-res solve
	if( _Ps.size() ) _Xs[0] = _solver->solve( _Bs[0] );
	else             x = _solver->solve( b );

	// Coarse-to-fine
	for( unsigned int l=1 ; l<=_Ps.size() ; l++ )
	{
		if( l==_Ps.size() )
		{
			x += _Ps[l-1] * _Xs[l-1];
			for( unsigned int i=0 ; i<iters ; i++ ) _gsRelaxers[l-1]( b , x , MCI(l) );
		}
		else
		{
			_Xs[l] += _Ps[l-1] * _Xs[l-1];
			for( unsigned int i=0 ; i<iters ; i++ ) _gsRelaxers[l-1]( _Bs[l] , _Xs[l] , MCI(l) );
		}
	}
}
#endif // NEW_CODE

#if 1 // NEW_CODE
template< typename Solver , bool Jacobi >
template< typename T >
void MGSolver< Solver , Jacobi >::_vCycle( const T & b , T & x , unsigned int iters )
{
	std::vector< T > _Xs( _Ms.size() ) , _Bs( _Ms.size() );
	if      constexpr( std::is_same_v< T , Eigen::MatrixXd > ) for( unsigned int i=0 ; i<_Ms.size() ; i++ ) _Xs[i].resize( _Ms[i].cols() , b.cols() );
	else if constexpr( std::is_same_v< T , Eigen::VectorXd > ) for( unsigned int i=0 ; i<_Ms.size() ; i++ ) _Xs[i].resize( _Ms[i].cols() );
	else static_assert( false , "[ERROR] Expected either Eigen vector or matrix type" );

	// Fine-to-coarse
	for( unsigned int l=(unsigned int)_Ps.size() ; l>0 ; l-- )
	{
		if( l==_Ps.size() )
		{
			for( unsigned int i=0 ; i<iters ; i++ ) _relaxers[l-1]( b , x );
			_Bs[l-1] = _Rs[l-1] * ( b - _M * x );
			_Xs[l-1].setZero();
		}
		else
		{
			for( unsigned int i=0 ; i<iters ; i++ ) _relaxers[l-1]( _Bs[l] , _Xs[l] );
			_Bs[l-1] = _Rs[l-1] * ( _Bs[l] - _Ms[l] * _Xs[l] );
			_Xs[l-1].setZero();
		}
	}

	// Low-res solve
	if( _Ps.size() ) _Xs[0] = _solver->solve( _Bs[0] );
	else             x = _solver->solve( b );

	// Coarse-to-fine
	for( unsigned int l=1 ; l<=_Ps.size() ; l++ )
	{
		if( l==_Ps.size() )
		{
			x += _Ps[l-1] * _Xs[l-1];
			for( unsigned int i=0 ; i<iters ; i++ ) _relaxers[l-1]( b , x );
		}
		else
		{
			_Xs[l] += _Ps[l-1] * _Xs[l-1];
			for( unsigned int i=0 ; i<iters ; i++ ) _relaxers[l-1]( _Bs[l] , _Xs[l] );
		}
	}
}

template< typename Solver , bool Jacobi >
template< typename T , typename MCIndicesFunctor /* = std::function< const std::vector< std::vector< unsigned int > > & ( unsigned int ) > */ >
void MGSolver< Solver , Jacobi >::_vCycle( const T & b , T & x , unsigned int iters , MCIndicesFunctor MCI )
{
	std::vector< T > _Xs( _Ms.size() ) , _Bs( _Ms.size() );

	if      constexpr( std::is_same_v< T , Eigen::MatrixXd > ) for( unsigned int i=0 ; i<_Ms.size() ; i++ ) _Xs[i].resize( _Ms[i].cols() , b.cols() );
	else if constexpr( std::is_same_v< T , Eigen::VectorXd > ) for( unsigned int i=0 ; i<_Ms.size() ; i++ ) _Xs[i].resize( _Ms[i].cols() );
	else static_assert( false , "[ERROR] Expected either Eigen vector or matrix type" );

	// Fine-to-coarse
	for( unsigned int l=(unsigned int)_Ps.size() ; l>0 ; l-- )
	{
		if( l==_Ps.size() )
		{
			if constexpr( Jacobi ) for( unsigned int i=0 ; i<iters ; i++ ) _relaxers[l-1]( b , x );
			else                   for( unsigned int i=0 ; i<iters ; i++ ) _relaxers[l-1]( b , x , MCI(l) );
			_Bs[l-1] = _Rs[l-1] * ( b - _M * x );
			_Xs[l-1].setZero();
		}
		else
		{
			if constexpr( Jacobi ) for( unsigned int i=0 ; i<iters ; i++ ) _relaxers[l-1]( _Bs[l] , _Xs[l] );
			else                   for( unsigned int i=0 ; i<iters ; i++ ) _relaxers[l-1]( _Bs[l] , _Xs[l] , MCI(l) );
			_Bs[l-1] = _Rs[l-1] * ( _Bs[l] - _Ms[l] * _Xs[l] );
			_Xs[l-1].setZero();
		}
	}
	// Low-res solve
	if( _Ps.size() ) _Xs[0] = _solver->solve( _Bs[0] );
	else             x = _solver->solve( b );

	// Coarse-to-fine
	for( unsigned int l=1 ; l<=_Ps.size() ; l++ )
	{
		if( l==_Ps.size() )
		{
			x += _Ps[l-1] * _Xs[l-1];
			if constexpr( Jacobi ) for( unsigned int i=0 ; i<iters ; i++ ) _relaxers[l-1]( b , x );
			else                   for( unsigned int i=0 ; i<iters ; i++ ) _relaxers[l-1]( b , x , MCI(l) );
		}
		else
		{
			_Xs[l] += _Ps[l-1] * _Xs[l-1];
			if constexpr( Jacobi ) for( unsigned int i=0 ; i<iters ; i++ ) _relaxers[l-1]( _Bs[l] , _Xs[l] );
			else                   for( unsigned int i=0 ; i<iters ; i++ ) _relaxers[l-1]( _Bs[l] , _Xs[l] , MCI(l) );
		}
	}
}
#endif // NEW_CODE
