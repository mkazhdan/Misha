#ifndef KDTREE_H
#define KDTREE_H
/*
Szymon Rusinkiewicz
Princeton University

KDtree.h
A K-D tree for points, with limited capabilities (find nearest point to 
a given point, or to a ray). 
*/

#include <vector>
#include <cmath>
#include <string.h>
#include <algorithm>
#include "mempool.h"

namespace SzymonRusinkiewicz
{
	template< class Real , int D >
	class KDtree
	{
	private:
		class Node;
		Node *root;
		void build( const Real *ptlist , int n );

	public:
		// Compatibility function for closest-compatible-point searches
		struct CompatFunc
		{
			virtual bool operator () ( const Real *p ) const = 0;
			virtual ~CompatFunc() {}  // To make the compiler shut up
		};

		// Constructor from an array of points
		KDtree( const Real *ptlist , int n ) { build( ptlist , n ); }
		~KDtree( void );

		// The queries: returns closest point to a point or a ray,
		// provided it's within sqrt(maxdist2) and is compatible
		const Real *closest_to_pt ( const Real *p ,                   Real maxdist2 , const CompatFunc *iscompat = NULL ) const;
		const Real *closest_to_ray( const Real *p , const Real *dir , Real maxdist2 , const CompatFunc *iscompat = NULL ) const;
	};


	using std::vector;
	using std::swap;
	using std::sqrt;


	// Small utility fcns
	template< class Real > static inline Real sqr( Real  x ){ return x*x; }
	template< class Real , int D > static inline Real dist2( const Real x[D] , const Real y[D] )
	{
		Real dist2 = 0;
		for( int i=0 ; i<D ; i++ ) dist2 += sqr( x[i] - y[i] );
		return dist2;
	}

	template< class Real , int D > static inline Real dist2ray2( const Real x[D] , const Real p[D] , const Real d[D] )
	{
		Real n2 = 0 , dot = 0;
		for( int i=0 ; i<D ; i++ )
		{
			Real xp = x[i] - p[i];
			n2 += xp * xp , dot += xp * d[i];
		}
		return n2 + dot * dot;
	}


	// Class for nodes in the K-D tree
	template< class Real , int D >
	class KDtree< Real , D >::Node
	{
	private:
		static PoolAlloc memPool;

	public:
		// A place to put all the stuff required while traversing the K-D
		// tree, so we don't have to pass tons of variables at each fcn call
		struct Traversal_Info
		{
			const Real *p, *dir;
			const Real *closest;
			Real closest_d, closest_d2;
			const typename KDtree< Real , D >::CompatFunc *iscompat;
		};

		enum { MAX_PTS_PER_NODE = 7 };


		// The node itself

		int npts; // If this is 0, intermediate node.  If nonzero, leaf.

		union
		{
			struct
			{
				Real center[D];
				Real r;
				int splitaxis;
				Node *child1, *child2;
			} node;
			struct
			{
				const Real *p[MAX_PTS_PER_NODE];
			} leaf;
		};

		Node( const Real **pts , int n );
		~Node();

		void find_closest_to_pt ( Traversal_Info &k ) const;
		void find_closest_to_ray( Traversal_Info &k ) const;

		void *operator new( size_t n ) { return memPool.alloc(n); }
		void operator delete( void *p , size_t n ) { memPool.free(p,n); }
	};


	// Class static variable
	template< class Real , int D > PoolAlloc KDtree< Real , D >::Node::memPool( sizeof(KDtree::Node) );


	// Create a KD tree from the points pointed to by the array pts
	template< class Real , int D >
	KDtree< Real , D >::Node::Node( const Real **pts , int n )
	{
		// Leaf nodes
		if( n<=MAX_PTS_PER_NODE )
		{
			npts = n;
			memcpy(leaf.p, pts, n * sizeof(Real *));
			return;
		}


		// Else, interior nodes
		npts = 0;

		// Find bbox
		Real min[D] , max[D];
		for( int d=0 ; d<D ; d++ ) min[d] = max[d] = pts[0][d];
		for( int i=1 ; i<n ; i++ ) for( int d=0 ; d<D ; d++ ) min[d] = std::min< Real >( min[d] , pts[i][d] ) , max[d] = std::max< Real >( max[d] , pts[i][d] );

		// Find node center and size
		for( int d=0 ; d<D ; d++ ) node.center[d] = (Real)0.5 * ( min[d] + max[d] );
		Real d[D] , r2 = 0;
		for( int i=0 ; i<D ; i++ ) d[i] = max[i] - min[i] , r2 += d[i] * d[i];
		node.r = (Real) 0.5 * sqrt( r2 );

		// Find longest axis
		node.splitaxis = D-1;
		for( int i=0 ; i<D-1 ; i++ ) if( d[i]>d[node.splitaxis] ) node.splitaxis = i;

		// Partition
		const Real splitval = node.center[node.splitaxis];
		const Real **left = pts, **right = pts + n - 1;
		while (1) {
			while ((*left)[node.splitaxis] < splitval)
				left++;
			while ((*right)[node.splitaxis] >= splitval)
				right--;
			if (right < left)
				break;
			swap(*left, *right);
			left++; right--;
		}

		// Check for bad cases of clustered points
		if (left-pts == 0 || left-pts == n)
			left = pts + n/2;

		// Build subtrees
		node.child1 = new Node(pts, int(left-pts) );
		node.child2 = new Node(left, int(n-(left-pts)));
	}


	// Destroy a KD tree node
	template< class Real , int D >
	KDtree< Real , D >::Node::~Node()
	{
		if (!npts) {
			delete node.child1;
			delete node.child2;
		}
	}


	// Crawl the KD tree
	template< class Real , int D >
	void KDtree< Real , D >::Node::find_closest_to_pt( typename KDtree< Real , D >::Node::Traversal_Info &k) const
	{
		// Leaf nodes
		if (npts) {
			for (int i = 0; i < npts; i++) {
				Real myd2 = dist2< Real , D >(leaf.p[i], k.p);
				if ((myd2 < k.closest_d2) &&
					(!k.iscompat || (*k.iscompat)(leaf.p[i]))) {
					k.closest_d2 = myd2;
					k.closest_d = sqrt(k.closest_d2);
					k.closest = leaf.p[i];
				}
			}
			return;
		}


		// Check whether to abort
		if( dist2< Real , D >(node.center, k.p) >= sqr(node.r + k.closest_d))
			return;

		// Recursive case
		Real myd = node.center[node.splitaxis] - k.p[node.splitaxis];
		if (myd >= 0.0f) {
			node.child1->find_closest_to_pt(k);
			if (myd < k.closest_d)
				node.child2->find_closest_to_pt(k);
		} else {
			node.child2->find_closest_to_pt(k);
			if (-myd < k.closest_d)
				node.child1->find_closest_to_pt(k);
		}
	}


	// Crawl the KD tree to look for the closest point to
	// the line going through k.p in the direction k.dir
	template< class Real , int D >
	void KDtree< Real , D >::Node::find_closest_to_ray( typename KDtree< Real , D >::Node::Traversal_Info &k ) const
	{
		// Leaf nodes
		if (npts) {
			for (int i = 0; i < npts; i++) {
				Real myd2 = dist2ray2< Real , D >(leaf.p[i], k.p, k.dir);
				if ((myd2 < k.closest_d2) &&
					(!k.iscompat || (*k.iscompat)(leaf.p[i]))) {
					k.closest_d2 = myd2;
					k.closest_d = sqrt(k.closest_d2);
					k.closest = leaf.p[i];
				}
			}
			return;
		}


		// Check whether to abort
		if (dist2ray2< Real , D >(node.center, k.p, k.dir) >= sqr(node.r + k.closest_d))
			return;

		// Recursive case
		if (k.p[node.splitaxis] < node.center[node.splitaxis] ) {
			node.child1->find_closest_to_ray(k);
			node.child2->find_closest_to_ray(k);
		} else {
			node.child2->find_closest_to_ray(k);
			node.child1->find_closest_to_ray(k);
		}
	}


	// Create a KDtree from a list of points (i.e., ptlist is a list of D*n floats)
	template< class Real , int D >
	void KDtree< Real , D >::build( const Real *ptlist , int n )
	{
		vector< const Real* > pts(n);
		for( int i=0 ; i<n ; i++ ) pts[i] = ( const Real* ) ( ( ( unsigned char* ) ptlist ) + i * sizeof( Real ) * D );
		root = new Node( &( pts[0] ) , n );
	}


	// Delete a KDtree
	template< class Real , int D > KDtree< Real , D >::~KDtree( void ) { delete root; }


	// Return the closest point in the KD tree to p
	template< class Real , int D >
	const Real *KDtree< Real , D >::closest_to_pt( const Real *p , Real maxdist2 , const CompatFunc *iscompat /* = NULL */ ) const
	{
		typename Node::Traversal_Info k;

		k.p = p;
		k.iscompat = iscompat;
		k.closest = NULL;
		if( maxdist2 <= 0.0f ) maxdist2 = sqr( root->node.r );
		k.closest_d2 = maxdist2;
		k.closest_d = sqrt( k.closest_d2 );

		root->find_closest_to_pt( k );

		return k.closest;
	}


	// Return the closest point in the KD tree to the line
	// going through p in the direction dir
	template< class Real , int D >
	const Real *KDtree< Real , D >::closest_to_ray(const Real *p, const Real *dir,
		Real maxdist2,
		const CompatFunc *iscompat /* = NULL */) const
	{
		fprintf( stderr , "[ERROR] Nearest point to ray not supported for kD-trees" ) , exit( 0 );

		typename Node::Traversal_Info k;

		Real one_over_dir_len = 1.0f / sqrt(sqr(dir[0])+sqr(dir[1])+sqr(dir[2]));
		Real normalized_dir[3] = { dir[0] * one_over_dir_len, 
			dir[1] * one_over_dir_len, 
			dir[2] * one_over_dir_len };
		k.dir = normalized_dir;
		k.p = p;
		k.iscompat = iscompat;
		k.closest = NULL;
		if (maxdist2 <= 0.0f)
			maxdist2 = sqr(root->node.r);
		k.closest_d2 = maxdist2;
		k.closest_d = sqrt(k.closest_d2);

		root->find_closest_to_ray(k);

		return k.closest;
	}
}

#endif
