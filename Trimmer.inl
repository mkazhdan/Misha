/*
Copyright (c) 2019, Michael Kazhdan
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

#include <unordered_map>

namespace Misha
{
	/////////////
	// Trimmer //
	/////////////
	template< typename Real , typename Vertex , unsigned int K , typename IndexType >
	void Trimmer::Trim( std::vector< Vertex > &vertices , std::vector< SimplexIndex< K , IndexType > > &simplices , const std::vector< Real > &values , Real isoValue , bool keepGreater )
	{
		if( vertices.size()!=values.size() ) ERROR_OUT( "Vertex and value count don't match: " , vertices.size() , "!=" , values.size() );

		EdgeTable< IndexType > edgeTable;
		std::vector< Vertex > _vertices;
		std::vector< SimplexIndex< K , IndexType > > _simplices;

		for( size_t i=0 ; i<simplices.size() ; i++ )
		{
			Real v[K+1];
			for( unsigned int k=0 ; k<=K ; k++ ) v[k] = values[ simplices[i][k] ] - isoValue;
			std::vector< SimplexIndex< K , IndexType > > back , front;
			simplices[i].split( v , vertices , edgeTable , back , front );
			if( keepGreater ) for( int i=0 ; i<front.size() ; i++ ) _simplices.push_back( front[i] );
			else              for( int i=0 ; i<back.size() ; i++ ) _simplices.push_back( back[i] );
		}

		std::vector< size_t > vertexMap( vertices.size() , -1 );
		for( size_t i=0 ; i<_simplices.size() ; i++ ) for( int k=0 ; k<=K ; k++ ) vertexMap[ _simplices[i][k] ] = 0;

		size_t count = 0;
		for( size_t i=0 ; i<vertices.size() ; i++ ) if( vertexMap[i]!=-1 ) vertexMap[i] = count++;

		_vertices.resize( count );
		for( size_t i=0 ; i<vertices.size() ; i++ ) if( vertexMap[i]!=-1 ) _vertices[ vertexMap[i] ] = vertices[i];

		for( size_t i=0 ; i<_simplices.size() ; i++ ) for( int k=0 ; k<=K ; k++ ) _simplices[i][k] = (IndexType)vertexMap[ _simplices[i][k] ];

		vertices = _vertices;
		simplices = _simplices;
	}
}