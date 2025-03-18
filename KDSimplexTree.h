#ifndef KD_SIMPLEX_TREE_INCLUDE
#define KD_SIMPLEX_TREE_INCLUDE

/*
* Code borrows _heavily_ from:
Szymon Rusinkiewicz
Princeton University

KDtree.h
A K-D tree for points, with limited capabilities (find nearest point to a given point). 
*/

#include <cmath>
#include "mempool.h"
#include "Geometry.h"
#include "MultiThreading.h"

namespace MishaK
{
	// Dim is the embedding dimension
	// K is the dimension of the simplex
	template< unsigned int Dim , unsigned int K , unsigned int MaxSimplicesPerNode=7 , typename Index=unsigned int >
	struct KDSimplexTree
	{
		static Point< double , Dim > ClosestPointOnSimplex( Point< double , Dim > p , Simplex< double , Dim , K > s )
		{
			if constexpr( K==0 ) return s[0];
			else
			{
				Point< double , Dim > q;
#if 1
				Point< double , K+1 > b = s.barycentricCoordinates( p );
#else
				Point< double , K+1 > b = Barycenter( p , s );
#endif
				for( unsigned int k=0 ; k<=K ; k++ )
					if( b[k]<0. )
					{
						Simplex< double , Dim , K-1 > _s;
						for( unsigned int i=0 ; i<K ; i++ ) _s[i] = s[(k+1+i)%K];
						return KDSimplexTree< Dim , K-1 , MaxSimplicesPerNode , Index >::ClosestPointOnSimplex( p , _s );
					}
					else q += s[k] * b[k];
				return q;
			}
		}

		struct BoundingBox
		{
			BoundingBox( void )
			{
				for( unsigned int d=0 ; d<Dim ; d++ ) bBox[0][d] = std::numeric_limits< double >::infinity() , bBox[1][d] = -std::numeric_limits< double >::infinity();
			}
			BoundingBox( const BoundingBox *bBoxes , unsigned int N ) : BoundingBox()
			{
				for( unsigned int n=0 ; n<N ; n++ ) (*this) += bBoxes[n];
			}
			template< unsigned int K >
			BoundingBox( Simplex< double , Dim , K > s )
			{
				bBox[0] = bBox[1] = s[0];
				for( unsigned int k=1 ; k<=K ; k++ ) (*this) += s[k];
			}
			Point< double , Dim > bBox[2];
			Point< double , Dim > center( void ) const { return ( bBox[0] + bBox[1] )/2; }
			Point< double , Dim > diagonal( void ) const { return bBox[1] - bBox[0]; }
			Point< double , Dim > &operator[] ( unsigned int idx ){ return bBox[idx]; }
			const Point< double , Dim > &operator[] ( unsigned int idx ) const { return bBox[idx]; }

			bool contains( Point< double , Dim > p ) const
			{
				for( unsigned int d=0 ; d<Dim ; d++ ) if( p[d]<bBox[0][d] || p[d]>=bBox[1][d] ) return false;
				return true;
			}
			BoundingBox operator + ( BoundingBox b ) const
			{
				BoundingBox _b;
				for( unsigned int d=0 ; d<Dim ; d++ ) _b.bBox[0][d] = std::min< double >( bBox[0][d] , b.bBox[0][d] ) , _b.bBox[1][d] = std::max< double >( bBox[1][d] , b.bBox[1][d] );
				return _b;
			}
			BoundingBox& operator += ( BoundingBox b )
			{
				for( unsigned int d=0 ; d<Dim ; d++ )
				{
					bBox[0][d] = std::min< double >( bBox[0][d] , b.bBox[0][d] );
					bBox[1][d] = std::max< double >( bBox[1][d] , b.bBox[1][d] );
				}
				return *this;
			}
			BoundingBox& operator += ( Point< double , Dim > p )
			{
				for( unsigned int d=0 ; d<Dim ; d++ ) bBox[0][d] = std::min< double >( bBox[0][d] , p[d] ) , bBox[1][d] = std::max< double >( bBox[1][d] , p[d] );
				return *this;
			}
		};

		struct SimplexInfo
		{
			SimplexInfo( void ){}
			SimplexInfo( Simplex< double , Dim , K > s , Index idx ) : s(s) , idx(idx) , bBox(s) {}
			Simplex< double , Dim , K > s;		// The simplex being stored
			Index idx;							// Its (original) index
			BoundingBox bBox;					// Its bounding box

			double operator[]( unsigned int d ) const { return ( bBox.bBox[0][d] + bBox.bBox[1][d] )/2.; }
		};

		// Constructor from a simplicial mesh represented by vertices and indices
		KDSimplexTree( const Point< double , Dim > *vertices , const SimplexIndex< K , Index > *simplices , size_t sCount )
		{
			_sInfo.resize( sCount );
			ThreadPool::ParallelFor
				(
					0 , sCount ,
					[&]( unsigned int , size_t i )
					{
						Simplex< double , Dim , K > s;
						for( unsigned int d=0 ; d<=Dim ; d++ ) s[d] = vertices[ simplices[i][d] ];
						_sInfo[i] = SimplexInfo( s , i );
					}
				);
			_root = new _Node( &_sInfo[0] , sCount );
		}
		~KDSimplexTree( void ){ delete _root; }

		std::pair< Index , Point< double , Dim > > closest( Point< double , Dim > p , double maxDist2=std::numeric_limits< double >::infinity() ) const;
	protected:
		// Class for nodes in the K-D tree
		struct _Node
		{
			// A place to put all the stuff required while traversing the K-D
			// tree, so we don't have to pass tons of variables at each fcn call
			struct TraversalInfo
			{
				Point< double , Dim > source , target;
				Index closest;
				double closest_d , closest_d2;
			};

			// The node itself

			size_t n; // If this is 0, intermediate node.  If nonzero, leaf.

			struct InteriorNode
			{
				BoundingBox bBox;
				unsigned int splitaxis;
				Point< double , Dim > center;
				double radius;
				_Node *child1 , *child2;
			};

			struct LeafNode{ SimplexInfo *sInfo[MaxSimplicesPerNode]; };

			union
			{
				InteriorNode iNode;
				LeafNode lNode;
			};

			// [NOTE] This constructor may permute the entries of the array
			_Node( SimplexInfo *sInfo , size_t n );
			~_Node();

			void find_closest_to_pt( TraversalInfo &k ) const;

			void *operator new( size_t n ){ return _memPool.alloc(n); }
			void operator delete( void *p , size_t n ) { _memPool.free(p,n); }
		protected:
			static PoolAlloc _memPool;
		};

		std::vector< SimplexInfo > _sInfo;
		_Node *_root;
	};

	// Class static variable
	template< unsigned int Dim , unsigned int K , unsigned int MaxSimplicesPerNode , typename Index >
	PoolAlloc KDSimplexTree< Dim , K , MaxSimplicesPerNode , Index >::_Node::_memPool( sizeof( KDSimplexTree::_Node ) );


	// Create a KD tree from the points pointed to by the array pts
	template< unsigned int Dim , unsigned int K , unsigned int MaxSimplicesPerNode , typename Index >
	KDSimplexTree< Dim , K , MaxSimplicesPerNode , Index >::_Node::_Node( SimplexInfo *sInfo , size_t n )
	{
		// Leaf nodes
		if( n<=MaxSimplicesPerNode )
		{
			this->n = n;
			for( unsigned int i=0 ; i<n ; i++ ) lNode.sInfo[i] = sInfo + i;
			return;
		}

		// Else, interior nodes
		this->n = 0;

		// Find the bounding box
		iNode.bBox = BoundingBox();
		for( unsigned int i=0 ; i<n ; i++ ) iNode.bBox += sInfo[i].bBox;

		// Find longest axis
		Point< double , Dim > diagonal = iNode.bBox.diagonal();
		iNode.center = iNode.bBox.center();
		iNode.radius = sqrt( diagonal.squareNorm() ) / 2.;
		iNode.splitaxis = 0;
		for( unsigned int dd=1 ; dd<Dim ; dd++ ) if( diagonal[dd]>diagonal[ iNode.splitaxis ] ) iNode.splitaxis = dd;

		// Partition
		double splitval = iNode.center[ iNode.splitaxis ];
		SimplexInfo *left = sInfo , *right = sInfo + n - 1;
		while( true )
		{
			while( (size_t)(left -sInfo)<n && (*left )[ iNode.splitaxis ] <  splitval ) left ++;
			while( (size_t)(right-sInfo)>0 && (*right)[ iNode.splitaxis ] >= splitval ) right--;
			if( right < left ) break;

			if( (*left)[ iNode.splitaxis ] == (*right)[ iNode.splitaxis ] )
			{
				// Several clustered equal points - ensure even split
				left += (right - left) / 2;
				break;
			}

			std::swap( *left , *right );
			left++ , right--;
		}

		// Check for bad cases of clustered points
		if( left-sInfo == 0 || left-sInfo == n ) left = sInfo + n/2;

		// Build subtrees
		iNode.child1 = new _Node( sInfo, (unsigned int)(left-sInfo) );
		iNode.child2 = new _Node( left , (unsigned int)(n-(left-sInfo)) );
	}


	// Destroy a KD tree node
	template< unsigned int Dim , unsigned int K , unsigned int MaxSimplicesPerNode , typename Index >
	KDSimplexTree< Dim , K , MaxSimplicesPerNode , Index >::_Node::~_Node()
	{
		if( !n )
		{
			delete iNode.child1;
			delete iNode.child2;
		}
	}


	// Crawl the KD tree
	template< unsigned int Dim , unsigned int K , unsigned int MaxSimplicesPerNode , typename Index >
	void KDSimplexTree< Dim , K , MaxSimplicesPerNode , Index >::_Node::find_closest_to_pt( TraversalInfo &k ) const
	{
		// Leaf nodes
		if( n )
		{
			for( unsigned int i=0 ; i<n ; i++ )
			{
				Point< double , Dim > p = ClosestPointOnSimplex( k.source , lNode.sInfo[i]->s );
				double myd2 = ( k.source - p ).squareNorm();
				if( myd2 < k.closest_d2 )
				{
					k.closest_d2 = myd2;
					k.closest = lNode.sInfo[i]->idx;
					k.target = p;
				}
			}
			k.closest_d = sqrt( k.closest_d2 );
			return;
		}

		// Check whether to abort
		// The closest the point can come to something inside the bounding with center c and radius r is:
		//		min(p) = || p - c || - r
		if( ( iNode.bBox.center() - k.source ).squareNorm() >= ( iNode.radius + k.closest_d ) * ( iNode.radius + k.closest_d ) ) return;
		iNode.child1->find_closest_to_pt( k );
		iNode.child2->find_closest_to_pt( k );
	}


	// Return the closest point in the KD tree to p
	template< unsigned int Dim , unsigned int K , unsigned int MaxSimplicesPerNode , typename Index >
	std::pair< Index , Point< double , Dim > > KDSimplexTree< Dim , K , MaxSimplicesPerNode , Index >::closest( Point< double , Dim > p , double maxDist2 ) const
	{
		typename _Node::TraversalInfo k;
		k.source = p;
		k.closest = -1;
		k.closest_d2 = maxDist2;
		k.closest_d = sqrt( k.closest_d2 );

		_root->find_closest_to_pt(k);

		return std::pair< Index , Point< double , Dim > >( k.closest , k.target );
	}
}
#endif // KD_SIMPLEX_TREE_INCLUDE
