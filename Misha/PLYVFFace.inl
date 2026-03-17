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

template< typename Index , typename Real , unsigned int Dim >
PlyVFFace< Index , Real , Dim >::PlyVFFace( void ) : _vertices(nullptr) , _sz(0) {}

template< typename Index , typename Real , unsigned int Dim >
PlyVFFace< Index , Real , Dim >::~PlyVFFace( void ){ resize(0); }

template< typename Index , typename Real , unsigned int Dim >
PlyVFFace< Index , Real , Dim >::PlyVFFace( const PlyVFFace & face )
{
	_vertices = nullptr;
	(*this) = face;
}

template< typename Index , typename Real , unsigned int Dim >
PlyVFFace< Index , Real , Dim > & PlyVFFace< Index , Real , Dim >::operator = ( const PlyVFFace& face )
{
	if( _vertices ) free( _vertices ) , _vertices = nullptr;
	_sz = face._sz;
	if( _sz ) _vertices = (Index*)malloc( sizeof(Index)*_sz );
	else      _vertices = nullptr;
	memcpy( _vertices , face._vertices , sizeof(Index)*_sz );
	v = face.v;
	return *this;
}

template< typename Index , typename Real , unsigned int Dim >
void PlyVFFace< Index , Real , Dim >::resize( unsigned int count )
{
	if( _vertices ) free( _vertices ) , _vertices = nullptr;
	_sz = 0;
	if( count ) _vertices = (Index*)malloc( sizeof(Index)*count ) , _sz = count;
}

template< typename Index , typename Real , unsigned int Dim >
unsigned int PlyVFFace< Index , Real , Dim >::size( void ) const { return _sz; }

template< typename Index , typename Real , unsigned int Dim >
Index& PlyVFFace< Index , Real , Dim >::operator[] ( unsigned int idx ){ return _vertices[idx]; }

template< typename Index , typename Real , unsigned int Dim >
const Index & PlyVFFace< Index , Real , Dim >::operator[] ( unsigned int idx ) const { return _vertices[idx]; }

template< typename Index , typename Real , unsigned int Dim >
std::vector< GregTurk::PlyProperty > PlyVFFace< Index , Real , Dim >::Properties( void )
{
	std::vector< GregTurk::PlyProperty > properties( Dim+1 );
	properties[0] = GregTurk::PlyProperty( "vertex_indices" , PLY::Type< Index >() ,  PLY::Type< Index >() , (int)offsetof( PlyVFFace , _vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyVFFace , _sz ) );
	for( unsigned int d=0 ; d<Dim ; d++ ) properties[d+1] = GregTurk::PlyProperty( "vf_" + std::to_string(d) , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyVFFace , v.coords ) + sizeof(Real)*d , 0 , 0 , 0 , 0 );
	return properties;
};

template< typename Index , typename Real , unsigned int Dim >
std::vector< GregTurk::PlyProperty > PlyVFFace< Index , Real , Dim >::ReadProperties( void )
{
	std::vector< GregTurk::PlyProperty > properties( Dim+1 );
	properties[0] = GregTurk::PlyProperty( "vertex_indices" , PLY::Type< Index >() ,  PLY::Type< Index >() , (int)offsetof( PlyVFFace , _vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyVFFace , _sz ) );
	for( unsigned int d=0 ; d<Dim ; d++ ) properties[d+1] = GregTurk::PlyProperty( "vf_" + std::to_string(d) , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyVFFace , v.coords ) + sizeof(Real)*d , 0 , 0 , 0 , 0 );
	if( Dim>0 ) properties.push_back( GregTurk::PlyProperty( "vx" , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyVFFace , v.coords ) + sizeof(Real)*0 , 0 , 0 , 0 , 0 ) );
	if( Dim>1 ) properties.push_back( GregTurk::PlyProperty( "vy" , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyVFFace , v.coords ) + sizeof(Real)*1 , 0 , 0 , 0 , 0 ) );
	if( Dim>2 ) properties.push_back( GregTurk::PlyProperty( "vz" , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyVFFace , v.coords ) + sizeof(Real)*2 , 0 , 0 , 0 , 0 ) );
	if( Dim>3 ) properties.push_back( GregTurk::PlyProperty( "vw" , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyVFFace , v.coords ) + sizeof(Real)*3 , 0 , 0 , 0 , 0 ) );
	return properties;
}
