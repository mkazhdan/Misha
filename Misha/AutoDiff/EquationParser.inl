/*
Copyright (c) 2025, Michael Kazhdan
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

/////////////////
// Node::State //
/////////////////
Node::State::State( NodeType type , std::string name ) : type(type) , name(name){}

//////////
// Node //
//////////
inline Node::Node( void ) : _type( NodeType::ZERO ) {}
inline Node::Node( double c ) : _type( NodeType::CONSTANT ) , _value(c) {}
inline Node Node::Variable( unsigned int idx ){ Node n ; n._type = NodeType::VARIABLE ; n._variableIndex = idx ; return n; }
inline Node Node::DVariable( unsigned int vIndex , unsigned int dIndex )
{
	Node node;
	if( vIndex==dIndex )
	{
		node._type = NodeType::CONSTANT;
		node._value = 1.;
	}
	else node._type = NodeType::ZERO;
	return node;
}

inline Node Node::Parse( std::string eqn , const std::vector< std::string > & vars )
{
	_StateInfo stateInfo( eqn , vars );
	return _Parse( stateInfo , vars );
}
	
inline Node Node::Function( NodeType type , const Node & n )
{
	Node node;
	node._children.push_back( n );
	node._type = type;
	switch( type )
	{
	case NodeType::NEGATION:           node._function = []( const double * values ){ return      -values[0]  ; } ; break;
	case NodeType::EXPONENTIAL:        node._function = []( const double * values ){ return exp ( values[0] ); } ; break;
	case NodeType::NATURAL_LOGARITHM:  node._function = []( const double * values ){ return log ( values[0] ); } ; break;
	case NodeType::COSINE:             node._function = []( const double * values ){ return cos ( values[0] ); } ; break;
	case NodeType::SINE:               node._function = []( const double * values ){ return sin ( values[0] ); } ; break;
	case NodeType::TANGENT:            node._function = []( const double * values ){ return tan ( values[0] ); } ; break;
	case NodeType::HYPERBOLIC_COSINE:  node._function = []( const double * values ){ return cosh( values[0] ); } ; break;
	case NodeType::HYPERBOLIC_SINE:    node._function = []( const double * values ){ return sinh( values[0] ); } ; break;
	case NodeType::HYPERBOLIC_TANGENT: node._function = []( const double * values ){ return tanh( values[0] ); } ; break;
	default: MK_THROW( "Node type is not a unary operator or function: " , NodeTypeNames[ static_cast< unsigned int >(type) ] );
	}
	return node;
}

inline Node Node::Function( NodeType type , const Node & node1 , const Node & node2 )
{
	Node node;
	node._children.push_back( node1 );
	node._children.push_back( node2 );
	node._type = type;
	switch( type )
	{
	case NodeType::ADDITION:       node._function = []( const double * values ){ return values[0]+values[1]; } ; break;
	case NodeType::SUBTRACTION:    node._function = []( const double * values ){ return values[0]-values[1]; } ; break;
	case NodeType::MULTIPLICATION: node._function = []( const double * values ){ return values[0]*values[1]; } ; break;
	case NodeType::DIVISION:       node._function = []( const double * values ){ return values[0]/values[1]; } ; break;
	case NodeType::POWER:          node._function = []( const double * values ){ return pow( values[0] , values[1] ); } ; break;
	default: MK_THROW( "Node type is not a binary operator: " , NodeTypeNames[ static_cast< unsigned int >(type) ] );
	}
	return node;
}

inline Node Node::DFunction( NodeType type , const Node & n , const Node & d )
{
	switch( type )
	{
	case NodeType::NEGATION:           return - d;
	case NodeType::EXPONENTIAL:        return   d * Exp( n );
	case NodeType::NATURAL_LOGARITHM:  return   d * Pow( n , -1. );
	case NodeType::COSINE:             return - d * Sin( n );
	case NodeType::SINE:               return   d * Cos( n );
	case NodeType::TANGENT:            return   d * ( 1. + Pow( Tan(n) , 2 ) );
	case NodeType::HYPERBOLIC_COSINE:  return   d * Sinh( n );
	case NodeType::HYPERBOLIC_SINE:    return   d * Cosh( n );
	case NodeType::HYPERBOLIC_TANGENT: return   d / Pow( Cosh(n) , 2 );
	default: MK_THROW( "Node type is not a unary operator or function: " , NodeTypeNames[ static_cast< unsigned int >(type) ] );
	};

	return Node();
}

inline Node Node::DFunction( NodeType type , const Node & node1 , const Node & dNode1 , const Node & node2 , const Node & dNode2 )
{
	switch( type )
	{
	case NodeType::ADDITION:       return dNode1 + dNode2;
	case NodeType::SUBTRACTION:    return dNode1 - dNode2;
	case NodeType::MULTIPLICATION: return ( dNode1 * node2 ) + ( node1 * dNode2 );
	case NodeType::DIVISION:       return ( dNode1 / node2 ) - ( ( node1 * dNode2 ) / Pow( node2 , 2 ) );
	case NodeType::POWER:
		if( node2._type==NodeType::CONSTANT ) return node2._value * dNode1 * Pow( node1 ,node2._value-1 );
		else                                  return Pow( node1 , node2 ) * ( Log( node1 ) * dNode2 + node2 / node1 * dNode1 );
	default: MK_THROW( "Node type is not a binary operator: " , NodeTypeNames[ static_cast< unsigned int >(type) ] );
	}
}

inline Node Node::_Parse( _StateInfo stateInfo , const std::vector< std::string > & vars )
{
	auto GetConstantNode = []( std::string str )
	{
		double v;
		try
		{
			if( str=="Pi" ) v = M_PI;
			else v = std::stod( str );
		}
		catch( ... ){ MK_THROW( "Failed to convert constant to double: " , str ); }
		return Node( v );
	};
	auto GetVariableNode = [&]( std::string str )
	{
		unsigned int idx = static_cast< unsigned int >(-1);
		for( unsigned i=0 ; i<vars.size() ; i++ ) if( vars[i]==str ) idx = i;
		if( idx==static_cast< unsigned int >(-1) ) MK_THROW( "Could not find variable: " , str );
		return Variable( idx );
	};

	Node node;
	if( stateInfo.state.size()==0 ) return node;
	else if( stateInfo.state.size()==1 )
	{
		if     ( stateInfo.state[0].type==_StateInfo::NodeType::CONSTANT ) return GetConstantNode( stateInfo.state[0].name );
		else if( stateInfo.state[0].type==_StateInfo::NodeType::VARIABLE ) return GetVariableNode( stateInfo.state[0].name );
		else MK_THROW( "Expected constant or variable" );
	}
	else
	{
		if( stateInfo.state.front().type!=_StateInfo::NodeType::L_PARENTH || stateInfo.state.back().type!=_StateInfo::NodeType::R_PARENTH ) MK_THROW( "Expected parenthesized state" );
		if( stateInfo.state.size()==3 ) return _Parse( stateInfo.sub(1,2) , vars );
		else
		{
			unsigned int idx = -1;
			unsigned int pCount = 0;
			for( unsigned int i=1 ; i<stateInfo.state.size()-1 ; i++ )
			{
				if     ( stateInfo.state[i].type==_StateInfo::NodeType::L_PARENTH ) pCount++;
				else if( stateInfo.state[i].type==_StateInfo::NodeType::R_PARENTH ) pCount--;
				else if( pCount==0 )
				{
					if( stateInfo.state[i].type==_StateInfo::NodeType::FUNCTION || stateInfo.state[i].type==_StateInfo::NodeType::UNARY_OPERATOR || stateInfo.state[i].type==_StateInfo::NodeType::BINARY_OPERATOR )
					{
						if( idx!=-1 ) MK_THROW( "Expected a single operator or function" );
						else idx = i;
					}
				}
			}
			if( idx==-1 ) MK_THROW( "Could not find operator or function" );

			NodeType type = static_cast< Node::NodeType >( stateInfo.state[idx] );
			if( stateInfo.state[idx].type==_StateInfo::NodeType::FUNCTION || stateInfo.state[idx].type==_StateInfo::NodeType::UNARY_OPERATOR )
			{
				if( idx+1>stateInfo.state.size() ) MK_THROW( "Expected argument after function" );
				unsigned int begin = idx+1;

				Node childNode;
				if     ( stateInfo.state[begin].type==_StateInfo::NodeType::CONSTANT  ) childNode = GetConstantNode( stateInfo.state[begin].name );
				else if( stateInfo.state[begin].type==_StateInfo::NodeType::VARIABLE  ) childNode = GetVariableNode( stateInfo.state[begin].name );
				else if( stateInfo.state[begin].type==_StateInfo::NodeType::L_PARENTH ) childNode = _Parse( stateInfo.sub( begin , stateInfo.closingParenth( begin )+1 ) , vars );
				else MK_THROW( "Expected constant, variable, or parentheses-enclosed equation" );
				return Function( type , childNode );
			}
			else if( stateInfo.state[idx].type==_StateInfo::NodeType::BINARY_OPERATOR )
			{
				if( idx==0 ) MK_THROW( "Expected content before binary operator" );
				if( idx+1>stateInfo.state.size() ) MK_THROW( "Expected content after binary operator" );

				Node childNode1 , childNode2;
				{
					unsigned int end = idx-1;
					if     ( stateInfo.state[end].type==_StateInfo::NodeType::CONSTANT  ) childNode1 = GetConstantNode( stateInfo.state[end].name );
					else if( stateInfo.state[end].type==_StateInfo::NodeType::VARIABLE  ) childNode1 = GetVariableNode( stateInfo.state[end].name );
					else if( stateInfo.state[end].type==_StateInfo::NodeType::R_PARENTH ) childNode1 = _Parse( stateInfo.sub( stateInfo.openingParenth( end ) , end+1 ) , vars );
					else MK_THROW( "Expected constant, variable, or parenthesized expression" );
				}
				{
					unsigned int begin = idx+1;
					if     ( stateInfo.state[begin].type==_StateInfo::NodeType::CONSTANT  ) childNode2 = GetConstantNode( stateInfo.state[begin].name );
					else if( stateInfo.state[begin].type==_StateInfo::NodeType::VARIABLE  ) childNode2 = GetVariableNode( stateInfo.state[begin].name );
					else if( stateInfo.state[begin].type==_StateInfo::NodeType::L_PARENTH ) childNode2 = _Parse( stateInfo.sub( begin , stateInfo.closingParenth( begin )+1 ) , vars );
					else MK_THROW( "Expected constant, variable, or parenthesized expression" );
				}
				return Function( type , childNode1 , childNode2 );
			}
			else MK_THROW( "Expected operator or function" );
		}
	}
	return node;
}

inline std::ostream & operator << ( std::ostream & stream , const Node & node )
{
	node._insert( stream , []( unsigned int idx ){ return std::string( "x" ) + std::to_string( idx ); } );
	return stream;
}

inline std::string Node::operator()( const std::vector< std::string > & varNames ) const
{
	std::stringstream sStream;
	_insert( sStream , [&]( unsigned int idx ){ return varNames[idx]; } );
	return sStream.str();
}

inline void Node::_insert( std::ostream & stream , const std::function< std::string ( unsigned int ) > & varNameFunctor ) const
{
	switch( _type )
	{
	case Node::NodeType::ZERO:
		stream << 0;
		break;
	case Node::NodeType::CONSTANT:
		stream << _value;
		break;
	case Node::NodeType::VARIABLE:
		stream << varNameFunctor( _variableIndex );
		break;
	default:
		if( IsUnaryOperator( _type ) )
		{
			stream << ToString( _type ) << "(";
			_children[0]._insert( stream , varNameFunctor );
			stream << ")";
		}
		else if( IsBinaryOperator( _type ) )
		{
			if( _children[0]._type==Node::NodeType::CONSTANT || _children[0]._type==Node::NodeType::VARIABLE ) _children[0]._insert( stream , varNameFunctor );
			else
			{
				stream << "(";
				_children[0]._insert( stream , varNameFunctor );
				stream << ")";
			}
			stream << ToString( _type );
			if( _children[1]._type==Node::NodeType::CONSTANT || _children[1]._type==Node::NodeType::VARIABLE ) _children[1]._insert( stream , varNameFunctor );
			else
			{
				stream << "(";
				_children[1]._insert( stream , varNameFunctor );
				stream << ")";
			}
		}
		else if( IsFunction( _type ) )
		{
			stream << ToString( _type ) << "(";
			_children[0]._insert( stream , varNameFunctor );
			stream << ")";
		}
		else MK_ERROR_OUT( "Unrecognized node type: " , NodeTypeNames[ static_cast< unsigned int >( _type ) ] );
	}
}

inline unsigned int Node::hasVariable( unsigned int idx ) const
{
	unsigned int count = ( _type==NodeType::VARIABLE && _variableIndex==idx ) ? 1 : 0;
	for( unsigned int i=0 ; i<_children.size() ; i++ ) count += _children[i].hasVariable( idx );
	return count;
}

inline bool Node::isConstant( void ) const
{
	bool isConstant = _type!=NodeType::VARIABLE;
	for( unsigned int i=0 ; i<_children.size() ; i++ ) isConstant &= _children[i].isConstant();
	return isConstant;
}

inline bool Node::isProduct( void ) const
{
	bool isProduct = _type==NodeType::VARIABLE || _type==NodeType::CONSTANT || _type==NodeType::ZERO || _type==NodeType::MULTIPLICATION || _type==NodeType::NEGATION;
	for( unsigned int i=0 ; i<_children.size() ; i++ ) isProduct &= _children[i].isProduct();
	return isProduct;
}

inline bool Node::__compress( void )
{
	auto IsNegation         = []( const Node & n ){ return n._type==NodeType::NEGATION; };
	auto IsZero             = []( const Node & n ){ return n._type==NodeType::ZERO; };
	auto IsOne              = []( const Node & n ){ return n._type==NodeType::CONSTANT && n._value==1; };
	auto IsMinusOne         = []( const Node & n ){ return n._type==NodeType::CONSTANT && n._value==-1; };
	auto IsScalar           = []( const Node & n ){ return n._type==NodeType::ZERO || n._type==NodeType::CONSTANT; };
	auto IsSum              = []( const Node & n ){ return n._type==NodeType::ADDITION; };
	auto IsDifference       = []( const Node & n ){ return n._type==NodeType::SUBTRACTION; };
	auto IsProduct          = []( const Node & n ){ return n._type==NodeType::MULTIPLICATION; };
	auto IsQuotient         = []( const Node & n ){ return n._type==NodeType::DIVISION; };
	auto IsExponent         = []( const Node & n ){ return n._type==NodeType::POWER; };
	auto IsScalarProduct    = []( const Node & n ){ return n._type==NodeType::MULTIPLICATION && ( n._children[0]._type==NodeType::CONSTANT || n._children[1]._type==NodeType::CONSTANT ); };
	auto IsVariable         = []( const Node & n ){ return n._type==NodeType::VARIABLE; };
	auto IsVariableExponent = []( const Node & n ){ return n._type==NodeType::POWER && n._children[0]._type==NodeType::VARIABLE; };

	if( isConstant() && _type!=NodeType::ZERO && _type!=NodeType::CONSTANT )
	{
		_value = this->operator()( nullptr );
		_type = _value==0 ? NodeType::ZERO : NodeType::CONSTANT;
		_children.resize(0);
		return true;
	}
	if( _type==NodeType::CONSTANT && !_value )
	{
		_type = NodeType::ZERO;
		return true;
	}

	// Collapse double negation
	if( IsNegation( *this ) && IsNegation( _children[0]  ) )
	{
		Node n = _children[0]._children[0];
		*this = n;
		return true;
	}

	// Reorder product terms
	if( IsProduct( *this ) )
	{
		if( _children[0]._type!=NodeType::CONSTANT && _children[1]._type==NodeType::CONSTANT )
		{
			std::swap( _children[0] , _children[1] );
			return true;
		}
		if( IsVariable( _children[0] ) && IsVariable( _children[1] ) )
		{
			if( _children[0]._variableIndex==_children[1]._variableIndex )
			{
				*this = _children[0] ^ 2;
				return true;
			}
			else if( _children[1]._variableIndex<_children[0]._variableIndex )
			{
				std::swap( _children[0] , _children[1] );
				return true;
			}
		}
		if( IsVariableExponent( _children[0] ) && IsVariableExponent( _children[1] ) )
		{
			if( _children[0]._children[0]._variableIndex==_children[1]._children[0]._variableIndex )
			{
				*this = _children[0]._children[0] ^ (_children[0]._children[1] + _children[1]._children[1] );
				return true;
			}
			else if( _children[1]._children[0]._variableIndex<_children[0]._children[0]._variableIndex )
			{
				std::swap( _children[0] , _children[1] );
				return true;
			}
		}
		if( IsVariable( _children[0] ) && IsVariableExponent( _children[1] ) )
		{
			if( _children[0]._variableIndex==_children[1]._children[0]._variableIndex )
			{
				*this = _children[0] ^ ( _children[1]._children[1] + 1. );
				return true;
			}
			else if( _children[1]._children[0]._variableIndex<_children[0]._variableIndex )
			{
				std::swap( _children[0] , _children[1] );
				return true;
			}
		}
		if( IsVariableExponent( _children[0] ) && IsVariable( _children[1] )  )
		{
			if( _children[1]._variableIndex==_children[0]._children[0]._variableIndex )
			{
				*this = _children[1] ^ ( _children[0]._children[1] + 1. );
				return true;
			}
			else if( _children[1]._variableIndex<_children[0]._children[0]._variableIndex )
			{
				std::swap( _children[0] , _children[1] );
				return true;
			}
		}
	}

	// Give precidence to quotients over products
	if( IsProduct( *this ) && ( IsQuotient( _children[0] ) || IsQuotient( _children[1] ) ) )
	{
		if( IsQuotient( _children[0] ) && IsQuotient( _children[1] ) )
			*this = ( _children[0]._children[0] * _children[1]._children[0] ) / ( _children[0]._children[1] * _children[1]._children[1] );
		else if( IsQuotient( _children[0] ) )
			*this = ( _children[0]._children[0] * _children[1] ) / _children[0]._children[1];
		else if( IsQuotient( _children[1] ) )
			*this = ( _children[0] * _children[1]._children[0] ) / _children[1]._children[1];
		return true;
	}


	// Collapse binary operators when (at least) one of the arguments is zero/negative
	if( IsBinaryOperator( _type ) )
	{
		bool z1 = IsZero( _children[0] ) , z2 = IsZero( _children[1] );
		bool n1 = IsNegation( _children[0] ) , n2 = IsNegation( _children[1] );
		bool s1 = IsScalar( _children[0] ) , s2 = IsScalar( _children[1] );
		bool sp1 = IsScalarProduct( _children[0] ) , sp2 = IsScalarProduct( _children[1] );

		if( IsSum( *this ) && ( z1 || z2 ) )
		{
			if( z1 && z2 ) *this = Node();
			else if( z1 ){ Node n = _children[1] ; *this = n; }
			else if( z2 ){ Node n = _children[0] ; *this = n; }
			return true;
		}
		if( IsDifference( *this ) && ( z1 || z2 ) )
		{
			if( z1 && z2 ) *this = Node();
			else if( z1 )  *this = -_children[1];
			else if( z2 ){ Node n =_children[0] ; *this = n; }
			return true;
		}
		if( IsProduct( *this ) && ( z1 || z2 ) )
		{
			*this = Node();
			return true;
		}
		if( IsQuotient( *this ) && ( z1 && !z2 ) )
		{
			*this = Node();
			return true;
		}
		if( IsExponent( *this ) )
		{
			if( z1 && z2 ) ;
			else if( z1 ){ *this = Node() ; return true; }
			else if( z2 ){ *this = Node(1.) ; return true; }
		}

		if( IsSum( *this ) && ( n1 || n2 ) )
		{
			if( n1 && n2 ) *this = - ( _children[0]._children[0] + _children[1]._children[0] );
			else if( n1 )  *this = - ( _children[0]._children[0] - _children[1] );
			else if( n2 )  *this =     _children[0] - _children[1]._children[0];
			return true;
		}
		if( IsDifference( *this ) && ( n1 || n2 ) )
		{
			if( n1 && n2 ) *this = - ( _children[0]._children[0] - _children[1]._children[0] );
			else if( n1 )  *this = - ( _children[0]._children[0] + _children[1] );
			else if( n2 )  *this =     _children[0] + _children[1]._children[0];
			return true;
		}
		if( IsProduct( *this ) && ( n1 || n2 ) )
		{
			if( n1 && n2 ) *this = _children[0]._children[0] * _children[1]._children[0];
			else if( n1 )  *this = -( _children[0]._children[0] * _children[1] );
			else if( n2 )  *this = -( _children[0] * _children[1]._children[0] );
			return true;
		}
		if( IsQuotient( *this ) && ( n1 || n2 ) )
		{
			if( n1 && n2 ) *this = _children[0]._children[0] / _children[1]._children[0];
			else if( n1 )  *this = - ( _children[0]._children[0] / _children[1] );
			else if( n2 )  *this = - ( _children[0] / _children[1]._children[0] );
			return true;
		}

		if( IsScalarProduct( *this ) && ( sp1 || sp2 ) )
		{
			if( sp1 && sp2 ) MK_ERROR_OUT( "Shouldn't happen" );
			else if( sp1 ) *this = ( _children[0]._children[0]._value * _children[1]._value ) * _children[0]._children[1];
			else if( sp2 ) *this = ( _children[0]._value * _children[1]._children[0]._value ) * _children[1]._children[1];
			return true;
		}

		if( IsProduct( *this ) && ( sp1 || sp2 ) )
		{
			if( sp1 && sp2 ) *this = ( _children[0]._children[0]._value * _children[1]._children[0]._value ) * ( _children[0]._children[1] * _children[1]._children[1] );
			else if( sp1 )   *this = ( _children[0]._children[0]._value ) * ( _children[0]._children[1] * _children[1] );
			else if( sp2 )   *this = ( _children[1]._children[0]._value ) * ( _children[0] * _children[1]._children[1] );
			return true;
		}

		if( IsProduct( *this ) && IsOne( _children[0] ) )
		{
			Node n = _children[1];
			*this = n;
			return true;
		}

		if( IsProduct( *this ) && IsOne( _children[1] ) )
		{
			Node n = _children[0];
			*this = n;
			return true;
		}

		if( IsProduct( *this ) && IsMinusOne( _children[0] ) )
		{
			*this = -_children[1];
			return true;
		}

		if( IsProduct( *this ) && IsMinusOne( _children[1] ) )
		{
			*this = -_children[0];
			return true;
		}

		if( IsQuotient( *this ) && IsOne( _children[1] ) )
		{
			Node n = _children[0];
			*this = n;
			return true;
		}

		if( IsProduct( *this ) && IsMinusOne( _children[1] ) )
		{
			*this = -_children[0];
			return true;
		}

		if( IsExponent( *this ) && IsOne( _children[1] ) )
		{
			Node n = _children[0];
			*this = n;
			return true;
		}

		if( IsExponent( *this ) && IsMinusOne( _children[1] ) )
		{
			*this = 1. / _children[0];
			return true;
		}
	}
	return false;
}

inline bool Node::_compress( void )
{
	bool compressed = false;
	for( unsigned int i=0 ; i<_children.size() ; i++ ) compressed |= _children[i]._compress();
	while( __compress() ) compressed = true;
	return compressed;
}

inline void Node::_sanityCheck( void )
{
	for( unsigned int i=0 ; i<_children.size() ; i++ ) _children[i]._sanityCheck();
	if     ( IsBinaryOperator( _type ) && _children.size()!=2 ) MK_THROW( "Expected two children " , _children.size() );
	else if( IsUnaryOperator ( _type ) && _children.size()!=1 ) MK_THROW( "Expected one child "    , _children.size() );
	else if( IsFunction      ( _type ) && _children.size()!=1 ) MK_THROW( "Expected one child: "   , _children.size() );
}

inline void Node::compress( void ){ while( _compress() ); }

template< unsigned int Dim >
double Node::operator()( Point< double , Dim > p ) const { return operator()( &p[0] ); }

inline double Node::operator()( const double * values ) const
{
	double _values[2];
	switch( _type )
	{
	case Node::NodeType::ZERO:     return 0;
	case Node::NodeType::CONSTANT: return _value;
	case Node::NodeType::VARIABLE: return values[ _variableIndex ];
	default:
		if( IsUnaryOperator( _type ) || IsFunction( _type ) )
		{
			_values[0] = _children[0]( values );
			return _function( _values );
		}
		else if( IsBinaryOperator( _type ) )
		{
			_values[0] = _children[0]( values );
			_values[1] = _children[1]( values );
			return _function( _values );
		}
		else
		{
			MK_THROW( "Urecognized node type: " , NodeTypeNames[ static_cast< unsigned int >( _type ) ] );
			return 0;
		}
	}
}

inline Node Node::d( unsigned int dIndex ) const
{
	Node node;
	if     ( _type==NodeType::ZERO || _type==NodeType::CONSTANT ) ;
	else if( _type==NodeType::VARIABLE ) node = DVariable( _variableIndex , dIndex );
	else if( IsUnaryOperator( _type ) || IsFunction( _type ) ) node = DFunction( _type , _children[0] , _children[0].d( dIndex ) );
	else if( IsBinaryOperator( _type ) ) node = DFunction( _type , _children[0] , _children[0].d(dIndex) , _children[1] , _children[1].d(dIndex) );
	else MK_THROW( "Unrecognized type: " , NodeTypeNames[ static_cast< unsigned int >( _type ) ] );
	node.compress();
	return node;
}

inline unsigned int Node::size( void ) const
{
	unsigned int sz = 1;
	for( unsigned int i=0 ; i<_children.size() ; i++ ) sz += _children[i].size();
	return sz;
}

/////////////////////////////
// Node::_StateInfo::State //
/////////////////////////////
Node::_StateInfo::State::State( NodeType type , std::string name ) : type(type) , name(name){}

Node::_StateInfo::State::operator Node::NodeType( void ) const
{
	Node::NodeType _type;
	switch( type )
	{
	case NodeType::VARIABLE: _type = Node::NodeType::VARIABLE  ; break;
	case NodeType::CONSTANT:  _type = Node::NodeType::CONSTANT ; break;
	case NodeType::UNARY_OPERATOR:
		if( name=="-" ) _type = Node::NodeType::NEGATION;
		else MK_THROW( "Unrecognized unary operator: " , name );
		break;
	case NodeType::BINARY_OPERATOR:
		if     ( name=="+" ) _type = Node::NodeType::ADDITION;
		else if( name=="-" ) _type = Node::NodeType::SUBTRACTION;
		else if( name=="*" ) _type = Node::NodeType::MULTIPLICATION;
		else if( name=="/" ) _type = Node::NodeType::DIVISION;
		else if( name=="^" ) _type = Node::NodeType::POWER;
		else MK_THROW( "Unrecognized binary operator: " , name );
		break;
	case NodeType::FUNCTION:
		if     ( name=="exp"  ) _type = Node::NodeType::EXPONENTIAL;
		else if( name=="cos"  ) _type = Node::NodeType::COSINE;
		else if( name=="sin"  ) _type = Node::NodeType::SINE;
		else if( name=="tan"  ) _type = Node::NodeType::TANGENT;
		else if( name=="cosh" ) _type = Node::NodeType::HYPERBOLIC_COSINE;
		else if( name=="sinh" ) _type = Node::NodeType::HYPERBOLIC_SINE;
		else if( name=="tanh" ) _type = Node::NodeType::HYPERBOLIC_TANGENT;
		else MK_THROW( "Unrecognized function: " , name );
		break;
	default:
		MK_THROW( "Could not convert node type: " , NodeTypeNames[ static_cast< int >( type ) ] );
	}
	return _type;
}

inline Node::_StateInfo::_StateInfo( void ){}

inline Node::_StateInfo::_StateInfo( std::string eqn , const std::vector< std::string > & vars )
{
	auto PrintType = []( const std::vector< NodeType > & state , const std::vector< std::string > & tokens )
		{
			for( unsigned int i=0 ; i<state.size() ; i++ )
				std::cout << "\t" << NodeTypeNames[ static_cast< int >( state[i] ) ] << " : " << tokens[i] << std::endl;
		};
	auto PrintState = []( const std::vector< State > & state , bool showType )
		{
			for( unsigned int i=0 ; i<state.size() ; i++ )
				if( showType ) std::cout << "\t" << NodeTypeNames[ static_cast< int >( state[i].type ) ] << " : " << state[i].name << std::endl;
				else           std::cout << " " << state[i].name;
			if( !showType ) std::cout << std::endl;
		};


	// Add spaces around operators and parentheses
	{
		std::string _eqn;
		for( unsigned int i=0 ; i<eqn.size() ; i++ )
		{
			if( IsOperator( eqn[i] ) || IsParenth( eqn[i] ) )
			{
				_eqn.push_back( ' ' );
				_eqn.push_back( eqn[i] );
				_eqn.push_back( ' ' );
			}
			else _eqn.push_back( eqn[i] );
		}
		eqn = _eqn;
	}

	// Tokenize on spaces
	std::vector< std::string > tokens;
	{
		auto RemoveLeadingWhiteSpace = []( const std::string str )
			{
				unsigned int idx = 0;
				while( idx<str.size() && std::isspace( str[idx] ) ) idx++;
				return str.substr( idx );
			};
		eqn = RemoveLeadingWhiteSpace( eqn );

		while( eqn.size() )
		{
			std::string token;
			for( unsigned int i=0 ; i<eqn.size() && !std::isspace( eqn[i] ) ; i++ ) token += eqn[i];
			if( !token.size() ) MK_THROW( "Expected non-empty token" );
			tokens.push_back( token );
			eqn = RemoveLeadingWhiteSpace( eqn.substr( token.size() ) );
		}
	}

	// Identify what we can
	std::vector< NodeType > _state( tokens.size() , NodeType::UNKNOWN );
	{
		int pCount = 0;
		for( unsigned int i=0 ; i<tokens.size() ; i++ )
		{
			if( IsOperator( tokens[i] ) ) _state[i] = NodeType::BINARY_OPERATOR;
			else if( tokens[i]=="(" || tokens[i]==")" )
			{
				if     ( tokens[i]=="(" ) pCount++ , _state[i] = NodeType::L_PARENTH;
				else if( tokens[i]==")" ) pCount-- , _state[i] = NodeType::R_PARENTH;
				if( pCount<0 ){ MK_THROW( "Negative parentheses count: " , eqn ); }
			}
		}
		if( pCount ){ MK_THROW( "Non-zero parentheses count" ); }
	}

	for( unsigned int i=0 ; i<tokens.size() ; i++ ) if( _state[i]==NodeType::UNKNOWN ) for( unsigned int j=0 ; j<vars.size() ; j++ ) if( tokens[i]==vars[j] )
		if( ( i==0 || _state[i-1]==NodeType::BINARY_OPERATOR || _state[i-1]==NodeType::L_PARENTH ) && ( i+1==tokens.size() || _state[i+1]==NodeType::BINARY_OPERATOR || _state[i+1]==NodeType::R_PARENTH ) )
		{
			_state[i] = NodeType::VARIABLE;
			break;
		}

	for( unsigned int i=0 ; i<_state.size() ; i++ )
	{
		if( _state[i]==NodeType::UNKNOWN )
		{
			try
			{
				if( tokens[i]=="Pi" ) state.emplace_back( NodeType::CONSTANT , tokens[i] );
				else
				{
					double foo = std::stod( tokens[i] );
					state.emplace_back( NodeType::CONSTANT , tokens[i] );
				}
			}
			catch( ... ){ state.emplace_back( NodeType::FUNCTION , tokens[i] ); }
		}
		else state.emplace_back( _state[i] , tokens[i] );
	}

	for( unsigned int i=0 ; i<state.size() ; i++ )
		if( state[i].type==NodeType::FUNCTION && ( i+1>=state.size() || state[i+1].type!=NodeType::L_PARENTH ) )
			MK_THROW( "Expected left parenth after function: " , i );

	addFunctionParenths();
	for( unsigned int i=0 ; i<UnaryOperators.size() ; i++ ) addUnaryOperatorParenths( UnaryOperators[i] );
	for( unsigned int i=0 ; i<BinaryOperators.size() ; i++ ) addBinaryOperatorParenths( BinaryOperators[i] );

	for( unsigned int i=0 ; i<state.size() ; i++ ) if( state[i].type==NodeType::UNKNOWN ) MK_THROW( "Failed to parse: " , i );
}

inline unsigned int Node::_StateInfo::openingParenth( unsigned int idx ) const
{
	if( state[idx].type!=NodeType::R_PARENTH ) MK_THROW( "Expected right parenth: " , idx );
	int pCount = 1;
	for( unsigned int i=idx ; i>0 ; i-- )
	{
		if     ( state[i-1].type==NodeType::L_PARENTH ) pCount--;
		else if( state[i-1].type==NodeType::R_PARENTH ) pCount++;
		if( pCount==0 ) return i-1;
	}
	MK_THROW( "Could not find opening parenth: " , idx );
}

inline unsigned int Node::_StateInfo::closingParenth( unsigned int idx ) const
{
	if( state[idx].type!=NodeType::L_PARENTH ) MK_THROW( "Expected left parenth: " , idx );
	int pCount = 1;
	for( unsigned int i=idx+1 ; i<state.size() ; i++ )
	{
		if     ( state[i].type==NodeType::L_PARENTH ) pCount++;
		else if( state[i].type==NodeType::R_PARENTH ) pCount--;
		if( pCount==0 ) return i;
	}
	MK_THROW( "Could not find closing parenth: " , idx );
	return static_cast< unsigned int >(-1);
}

inline Node::_StateInfo Node::_StateInfo::sub( unsigned int begin , unsigned end ) const
{
	_StateInfo stateInfo;
	for( unsigned int i=begin ; i<end ; i++ ) stateInfo.state.emplace_back( state[i] );
	return stateInfo;
}

inline void Node::_StateInfo::addFunctionParenths( void )
{
	for( unsigned int i=0 ; i<state.size() ; i++ )
		if( state[i].type==NodeType::FUNCTION && ( i==0 || ( i>0 && state[i-1].type!=NodeType::L_PARENTH ) ) )
		{
			unsigned int idx = closingParenth( i+1 );

			std::vector< State > _state;
			for( unsigned int j=0 ; j<i ; j++ ) _state.push_back( state[j] );
			_state.emplace_back( NodeType::L_PARENTH , "(" );
			for( unsigned int j=i ; j<idx ; j++ ) _state.push_back( state[j] );
			_state.emplace_back( NodeType::R_PARENTH , ")" );
			for( unsigned int j=idx ; j<state.size() ; j++ ) _state.push_back( state[j] );
			state = _state;
			addFunctionParenths();
		}
}

inline void Node::_StateInfo::addUnaryOperatorParenths( const std::vector< std::string > & ops )
{
	for( unsigned int i=0 ; i<state.size() ; i++ )
	{
		if( state[i].type==NodeType::BINARY_OPERATOR && IsOperator( state[i].name , ops ) )
		{
			if( i && ( state[i-1].type==NodeType::VARIABLE || state[i-1].type==NodeType::CONSTANT || state[i-1].type==NodeType::R_PARENTH ) ) continue;
			unsigned int end;

			if( state[i+1].type==NodeType::VARIABLE || state[i+1].type==NodeType::CONSTANT ) end = i+1;
			else if( state[i+1].type==NodeType::L_PARENTH ) end = closingParenth( i+1 );
			else MK_THROW( "Expected variable or left parenth" );

			if( i==0 || state[i-1].type!=NodeType::L_PARENTH || closingParenth(i-1)!=end+1 )
			{
				std::vector< State > _state;
				for( unsigned int j=0 ; j<i ; j++ ) _state.push_back( state[j] );
				_state.emplace_back( NodeType::L_PARENTH , "(" );
				for( unsigned int j=i ; j<=end ; j++ ) _state.push_back( state[j] );
				_state.emplace_back( NodeType::R_PARENTH , ")" );
				for( unsigned int j=end+1 ; j<state.size() ; j++ ) _state.push_back( state[j] );
				_state[i+1].type = NodeType::UNARY_OPERATOR;
				state = _state;
				addUnaryOperatorParenths( ops );
			}
			else state[i].type = NodeType::UNARY_OPERATOR;
		}
	}
}

inline void Node::_StateInfo::addBinaryOperatorParenths( const std::vector< std::string > & ops )
{
	for( unsigned int i=0 ; i<state.size() ; i++ ) if( state[i].type==NodeType::BINARY_OPERATOR && IsOperator( state[i].name , ops ) )
	{
		unsigned int begin , end;

		if( state[i-1].type==NodeType::VARIABLE || state[i-1].type==NodeType::CONSTANT ) begin = i-1;
		else if( state[i-1].type==NodeType::R_PARENTH ) begin = openingParenth( i-1 );
		else MK_THROW( "Expected variable or right parenth" );

		if( state[i+1].type==NodeType::VARIABLE || state[i+1].type==NodeType::CONSTANT ) end = i+1;
		else if( state[i+1].type==NodeType::L_PARENTH ) end = closingParenth( i+1 );
		else MK_THROW( "Expected variable or left parenth" );

		if( begin==0 || state[begin-1].type!=NodeType::L_PARENTH || closingParenth( begin-1 )!=end+1 )
		{
			std::vector< State > _state;
			for( unsigned int j=0 ; j<begin ; j++ ) _state.push_back( state[j] );
			_state.emplace_back( NodeType::L_PARENTH , "(" );
			for( unsigned int j=begin ; j<=end ; j++ ) _state.push_back( state[j] );
			_state.emplace_back( NodeType::R_PARENTH , ")" );
			for( unsigned int j=end+1 ; j<state.size() ; j++ ) _state.push_back( state[j] );
			state = _state;
			addBinaryOperatorParenths( ops );
		}
	}
}

inline bool Node::_StateInfo::IsParenth( char c ){ return c=='(' || c==')'; }

inline bool Node::_StateInfo::IsOperator( char c ){ return IsOperator( std::string( 1 , c ) ); }
inline bool Node::_StateInfo::IsOperator( std::string str )
{
	return IsUnaryOperator( str ) || IsBinaryOperator( str );
}

inline bool Node::_StateInfo::IsOperator( char c , const std::vector< std::string > & ops ){ return IsOperator( std::string( 1 , c ) , ops ); }
inline bool Node::_StateInfo::IsOperator( std::string str , const std::vector< std::string > & ops )
{
	for( unsigned int i=0 ; i<ops.size() ; i++ ) if( str==ops[i] ) return true;
	return false;
}

inline bool Node::_StateInfo::IsUnaryOperator( char c ){ return IsUnaryOperator( std::string( 1 , c ) ); }
inline bool Node::_StateInfo::IsUnaryOperator( std::string str )
{
	for( unsigned int i=0 ; i<UnaryOperators.size() ; i++ ) if( IsOperator( str , UnaryOperators[i] ) ) return true;
	return false;
}

inline bool Node::_StateInfo::IsBinaryOperator( char c ){ return IsBinaryOperator( std::string( 1 , c ) ); }
inline bool Node::_StateInfo::IsBinaryOperator( std::string str )
{
	for( unsigned int i=0 ; i<BinaryOperators.size() ; i++ ) if( IsOperator( str , BinaryOperators[i] ) ) return true;
	return false;
}

inline Node & operator += ( Node & n , const Node & _n ){ return n = n + _n; }
inline Node & operator -= ( Node & n , const Node & _n ){ return n = n - _n; }
inline Node & operator += ( Node & n , double s ){ return n = n + s; }
inline Node & operator -= ( Node & n , double s ){ return n = n - s; }
inline Node & operator *= ( Node & n , const Node & _n ){ return n = n * _n; }
inline Node & operator /= ( Node & n , const Node & _n ){ return n = n / _n; }
inline Node & operator *= ( Node & n , double s ){ return n = n * s; }
inline Node & operator /= ( Node & n , double s ){ return n = n / s; }
inline Node operator - ( const Node & node ){ return Node::Function( Node::NodeType::NEGATION , node ); }
inline Node operator + ( const Node & node1 , const Node & node2 ){ return Node::Function( Node::NodeType::ADDITION , node1 , node2 ); }
inline Node operator + ( const Node & node , double s ){ return Node::Function( Node::NodeType::ADDITION , node , Node( s ) ); }
inline Node operator + ( double s , const Node & node ){ return Node::Function( Node::NodeType::ADDITION , Node( s ) , node ); }
inline Node operator - ( const Node & node1 , const Node & node2 ){ return Node::Function( Node::NodeType::SUBTRACTION , node1 , node2 ); }
inline Node operator - ( const Node & node , double s ){ return Node::Function( Node::NodeType::SUBTRACTION , node , Node( s ) ); }
inline Node operator - ( double s , const Node & node ){ return Node::Function( Node::NodeType::SUBTRACTION , Node( s ) , node ); }
inline Node operator * ( const Node & node1 , const Node & node2 ){ return Node::Function( Node::NodeType::MULTIPLICATION , node1 , node2 ); }
inline Node operator * ( const Node & node , double s ){ return Node::Function( Node::NodeType::MULTIPLICATION , node , Node( s ) ); }
inline Node operator * ( double s , const Node & node ){ return Node::Function( Node::NodeType::MULTIPLICATION , Node( s ) , node ); }
inline Node operator / ( const Node & node1 , const Node & node2 ){ return Node::Function( Node::NodeType::DIVISION , node1 , node2 ); }
inline Node operator / ( const Node & node , double s ){ return Node::Function( Node::NodeType::DIVISION , node , Node( s ) ); }
inline Node operator / ( double s , const Node & node ){ return Node::Function( Node::NodeType::DIVISION , Node( s ) , node ); }
inline Node operator ^ ( const Node & node1 , const Node & node2 ){ return Node::Function( Node::NodeType::POWER , node1 , node2 ); }
inline Node operator ^ ( const Node & node , double s ){ return Node::Function( Node::NodeType::POWER , node , Node( s ) ); }
inline Node operator ^ ( double s , const Node & node ){ return Node::Function( Node::NodeType::POWER , Node( s ) , node ); }

inline Node Pow( const Node & n1 , const Node & n2 ){ return Node::Function( Node::NodeType::POWER , n1 , n2 ); }
inline Node Pow( const Node & n , double s ){ return Node::Function( Node::NodeType::POWER , n , s ); }
inline Node Pow( double s , const Node & n ){ return Node::Function( Node::NodeType::POWER , s , n ); }
inline Node Exp( const Node & n ){ return Node::Function( Node::NodeType::EXPONENTIAL , n ); }
inline Node Log( const Node & n ){ return Node::Function( Node::NodeType::NATURAL_LOGARITHM , n ); }
inline Node Cos( const Node & n ){ return Node::Function( Node::NodeType::COSINE , n ); }
inline Node Sin( const Node & n ){ return Node::Function( Node::NodeType::SINE , n ); }
inline Node Tan( const Node & n ){ return Node::Function( Node::NodeType::TANGENT , n ); }
inline Node Cosh( const Node & n ){ return Node::Function( Node::NodeType::HYPERBOLIC_COSINE , n ); }
inline Node Sinh( const Node & n ){ return Node::Function( Node::NodeType::HYPERBOLIC_SINE , n ); }
inline Node Tanh( const Node & n ){ return Node::Function( Node::NodeType::HYPERBOLIC_TANGENT , n ); }
inline Node Sqrt( const Node & n ){ return n ^ 0.5; }
