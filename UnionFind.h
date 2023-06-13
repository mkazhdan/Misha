#ifndef UNION_FIND_INCLUDED
#define UNION_FIND_INCLUDED
namespace Misha
{
	class UnionFind
	{
		struct Node
		{
			size_t sz;
			Node *root , *previous , *tail;
			Node( void ) { root = previous = tail = NULL , sz=0; }
		};
		Node* _nodes;
	public:
		UnionFind( void ) { _nodes = NULL; }
		UnionFind( size_t sz )
		{
			_nodes = new Node[ sz ];
			for( size_t i=0 ; i<sz ; i++ ) _nodes[i].root = _nodes[i].tail = _nodes+i , _nodes[i].previous = NULL , _nodes[i].sz=1;
		};
		~UnionFind( void ){ if( _nodes ) delete[] _nodes; }
		size_t find( size_t idx ) const { return _nodes[idx].root - _nodes; }
		size_t size( size_t idx ) const { return _nodes[ find(idx) ].sz; }
		size_t merge( size_t idx1 , size_t idx2 )
		{
			size_t _idx1 , _idx2;
			size_t sz1 = size( idx1 ) , sz2 = size( idx2 );
			if( sz1<sz2 ) _idx1 = find( idx1 ) , _idx2 = find( idx2 );
			else          _idx1 = find( idx2 ) , _idx2 = find( idx1 );
			// _idx1 points to the smaller of the two sets
			{
				Node* p = _nodes + _idx2;
				Node* c = _nodes + _idx1;
				p->sz = sz1 + sz2;
				_nodes[_idx2].tail->previous = c;
				_nodes[_idx2].tail = c->tail;
				while( c ) { c->root = p , c = c->previous; }
			}
			return sz1 + sz2;
		}
	};
};
#endif // UNION_FIND_INCLUDED