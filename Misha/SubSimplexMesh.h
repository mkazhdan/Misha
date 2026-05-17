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

#ifndef SUB_SIMPLEX_MESH_INCLUDED
#define SUB_SIMPLEX_MESH_INCLUDED


#include <vector>
#include <Misha/Geometry.h>

namespace MishaK
{
	template< unsigned int K , unsigned int SubK >
	struct SubSimplexMesh
	{
		static_assert( SubK<=K, "[ERROR] Sub-dimension exceeds dimension" );

		using SignedIndex = std::pair< size_t , bool >;

		// Constructor:
		// [NOTE] The object accesses the reference to the simplices throughout its scope
		SubSimplexMesh( const std::vector< SimplexIndex< K > > & simplices );

		// Indexing functionality
		size_t size( void ) const;
		SimplexIndex< SubK > operator[]( size_t i ) const;

		// Get the index associated with the specified simplex of the specified simplex
		SignedIndex operator()( size_t s , unsigned int n )  const;

	protected:
		using _SubSimplexIndex = MultiIndex< SubK+1 , size_t , true >;

		static const typename SimplexIndex< K >::template Faces< SubK > _SubSimplices;

		const std::vector< SimplexIndex< K > > & _simplices;
		std::map< _SubSimplexIndex , size_t > _subSimplexMap;
		std::vector< _SubSimplexIndex > _subSimplexIndices;
	};
#include "SubSimplexMesh.inl"
}
#endif // SUB_SIMPLEX_MESH_INCLUDED