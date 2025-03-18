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

#ifndef RASTERIZER_INCLUDED
#define RASTERIZER_INCLUDED

#include <mutex>
#include "Array.h"
#include "Geometry.h"
#include "RegularGrid.h"
#include "MultiThreading.h"

namespace MishaK
{
	template< typename Real , unsigned int Dim >
	struct Rasterizer
	{
		template< typename IndexType , unsigned int K > using SimplexRasterizationGrid = RegularGrid< std::vector< std::pair< IndexType , Simplex< Real , Dim , K > > > , Dim >;

		// This templated function rasterizes a simplicial complex into the cube [0,1]^3.
		// Template parameters:
		//		IndexType: specifies the storage for vertex/simplex indices
		// Input:
		//		vertices: the vertices of the mesh
		//		simplices: the connectivity of the mesh
		//		depth: the depth of the voxel grid, generating a grid of size (2^depth) x (2^depth) x (2^depth)
		//		lockDepth: the depth of the voxel grid storing the locks
		//		bBoxScale: The ratio of unit-cube-size to bounding-box-size
		// Output:
		//		A RegularGrid object where each cell stores the list of pairs containing the index into the original simplex list and the (clipped) simplex
		//		gridToModel: the transformation taking the grid coordinates back to the model coordinates (composed of an isotropic scale and a translation)
		template< typename IndexType , unsigned int K >
		static SimplexRasterizationGrid< IndexType , K > Rasterize( const SimplicialComplex< Real , Dim , K > &simplicialComplex , unsigned int depth , unsigned int lockDepth , XForm< Real , Dim+1 > &unitCubeToModel , Real bBoxScale );
		template< typename IndexType , unsigned int K >
		static SimplexRasterizationGrid< IndexType , K > Rasterize( const SimplicialComplex< Real , Dim , K > &simplicialComplex , unsigned int depth , unsigned int lockDepth ,                                          Real bBoxScale );

	protected:
		struct _RegularGridIndex
		{
			unsigned int depth , index[Dim];
			_RegularGridIndex( void );
			_RegularGridIndex( unsigned int d , Point< Real , Dim > p );
			template< unsigned int K > _RegularGridIndex( unsigned int maxDepth , Simplex< Real , Dim , K > simplex );

			bool operator != ( const _RegularGridIndex &idx ) const;
			bool operator == ( const _RegularGridIndex &idx ) const { return !( (*this)!=idx ); }

			_RegularGridIndex child( unsigned int c ) const;
		};

		struct _RegularGridLocks
		{
			_RegularGridLocks( unsigned int lockDepth , unsigned int maxDepth )
			{
				if( lockDepth>maxDepth )
				{
					WARN( "Lock depth exceeds max depth: " , lockDepth , " > " ,  maxDepth );
					lockDepth = maxDepth;
				}
				_bitShift = maxDepth - lockDepth;
				unsigned int _res = 1<<lockDepth;
				unsigned int res[Dim];
				for( int d=0 ; d<Dim ; d++ ) res[d] = _res;
				_locks.resize( res );

			}

			std::mutex &operator() ( const unsigned int idx[Dim] )
			{
				unsigned int _idx[Dim];
				for( int d=0 ; d<Dim ; d++ ) _idx[d] = idx[d] >> _bitShift;
				return _locks( _idx );
			}
			std::mutex &operator() ( unsigned int idx[Dim] )
			{
				unsigned int _idx[Dim];
				for( int d=0 ; d<Dim ; d++ ) _idx[d] = idx[d] >> _bitShift;
				return _locks( _idx );
			}
			template< typename ... UnsignedInts >
			std::mutex &operator()( UnsignedInts ... idx )
			{
				unsigned int _idx[] = { idx ... };
				for( int d=0 ; d<Dim ; d++ ) _idx[d] = _idx[d] >> _bitShift;
				return _locks( _idx );
			}

		protected:
			RegularGrid< std::mutex , Dim > _locks;
			size_t _bitShift;
		};

		template< typename IndexType , unsigned int K >
		static size_t _Rasterize( _RegularGridLocks &locks , SimplexRasterizationGrid< IndexType , K > &raster , IndexType simplexIndex , Simplex< Real , Dim , K > simplex , unsigned int depth , _RegularGridIndex idx );

		template< unsigned int K >
		static XForm< Real , Dim+1 > _ModelToUnitCube( const SimplicialComplex< Real , Dim , K > &simplicialComplex , Real bBoxScale );
	};
#include "Rasterizer.inl"
}

#endif // RASTERIZER_INCLUDED