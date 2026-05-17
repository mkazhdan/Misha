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
// SubSimplexMesh //
////////////////////
template< unsigned int K , unsigned int SubK >
const typename SimplexIndex< K >::template Faces< SubK > SubSimplexMesh< K , SubK >::_SubSimplices;

template< unsigned int K , unsigned int SubK >
SubSimplexMesh< K , SubK >::SubSimplexMesh( const std::vector< SimplexIndex< K > > & simplices ) : _simplices( simplices )
{
	for( unsigned int i=0 ; i<_simplices.size() ; i++ )
	{
		size_t idx[SubK+1];
		for( unsigned int ii=0 ; ii<SimplexIndex< K >::template FaceNum< SubK >() ; ii++ )
		{
			for( unsigned int k=0 ; k<=SubK ; k++ ) idx[k] = simplices[i][ _SubSimplices[ii][k] ];
			_SubSimplexIndex fi( idx );
			auto iter = _subSimplexMap.find( fi );
			if( iter==_subSimplexMap.end() ) _subSimplexMap[ fi ] = _subSimplexMap.size();
		}
	}
	_subSimplexIndices.resize( _subSimplexMap.size() );
	for( auto iter : _subSimplexMap ) _subSimplexIndices[ iter.second ] = iter.first;
}

template< unsigned int K , unsigned int SubK >
size_t SubSimplexMesh< K , SubK >::size( void ) const { return _subSimplexIndices.size(); }

template< unsigned int K , unsigned int SubK >
SimplexIndex< SubK > SubSimplexMesh< K , SubK >::operator[]( size_t i ) const
{
	SimplexIndex< SubK > si;
	for( unsigned int k=0 ; k<=SubK ; k++ ) si[k] = static_cast< unsigned int >( _subSimplexIndices[i][k] );
	return si;
}

template< unsigned int K , unsigned int SubK >
typename SubSimplexMesh< K , SubK >::SignedIndex SubSimplexMesh< K , SubK >::operator()( size_t s , unsigned int n ) const
{
	size_t idx[SubK+1];
	for( unsigned int k=0 ; k<=SubK ; k++ ) idx[k] = _simplices[s][ _SubSimplices[n][k] ];
	auto iter = _subSimplexMap.find( _SubSimplexIndex( idx ) );
	if( iter==_subSimplexMap.end() ) MK_THROW( "Failed to find sub-simplex: " , _SubSimplexIndex(idx) );
	return SignedIndex( iter->second , _SubSimplexIndex::Sign( idx ) );
}

