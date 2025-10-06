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

inline Node Node::Parse( std::string eqn , const std::vector< std::string > & vars )
{
	_StateInfo stateInfo( eqn , vars );
	return _Parse( stateInfo , vars );
}
	
inline Node Node::_GetConstant( std::string str )
{
#if 1 // NEW_CODE
	if( str=="Pi" ) return _GetConstant( M_PI );
	else return _GetConstant( std::stod(str) );
#else // !NEW_CODE
	return _GetConstant( std::stod(str) );
#endif // NEW_CODE
}

inline Node Node::_GetConstant( double value )
{
	Node node;
	node._type = _StateInfo::ConvertNodeType( _StateInfo::NodeType::CONSTANT );
	node._value = value;
	return node;
}

inline Node Node::_GetVariable( std::string str , const std::vector< std::string > & vars )
{
	auto VariableIndex = [&]( const std::string & name )
		{
			for( unsigned int i=0 ; i<vars.size() ; i++ ) if( name==vars[i] ) return i;
			MK_THROW( "Could not find variable in variable list: " , name );
			return static_cast< unsigned int >(-1);
		};

	Node node;
	node._type = _StateInfo::ConvertNodeType( _StateInfo::NodeType::VARIABLE );
	node._variableIndex = VariableIndex( str );
	return node;
}

inline Node Node::_GetFunction( std::string str , const Node & n )
{
	Node node;
	node._children.push_back( n );
	node._functionName = str;
	if     ( str=="-"    ) node._function = []( const double * values ){ return      -values[0]  ; } , node._type = NodeType::UNARY_OPERATOR;
	else if( str=="exp"  ) node._function = []( const double * values ){ return exp ( values[0] ); } , node._type = NodeType::FUNCTION;
	else if( str=="log"  ) node._function = []( const double * values ){ return log ( values[0] ); } , node._type = NodeType::FUNCTION;
	else if( str=="cos"  ) node._function = []( const double * values ){ return cos ( values[0] ); } , node._type = NodeType::FUNCTION;
	else if( str=="sin"  ) node._function = []( const double * values ){ return sin ( values[0] ); } , node._type = NodeType::FUNCTION;
	else if( str=="tan"  ) node._function = []( const double * values ){ return tan ( values[0] ); } , node._type = NodeType::FUNCTION;
	else if( str=="cosh" ) node._function = []( const double * values ){ return cosh( values[0] ); } , node._type = NodeType::FUNCTION;
	else if( str=="sinh" ) node._function = []( const double * values ){ return sinh( values[0] ); } , node._type = NodeType::FUNCTION;
	else if( str=="tanh" ) node._function = []( const double * values ){ return tanh( values[0] ); } , node._type = NodeType::FUNCTION;
	else MK_THROW( "Failed to parse function: " , str );
	return node;
}

inline Node Node::_GetFunction( std::string str , const Node & node1 , const Node & node2 )
{
	Node node;
	node._children.push_back( node1 );
	node._children.push_back( node2 );
	node._functionName = str;
	if     ( str=="+" ) node._function = []( const double * values ){ return values[0]+values[1]; }          , node._type = NodeType::BINARY_OPERATOR;
	else if( str=="-" ) node._function = []( const double * values ){ return values[0]-values[1]; }          , node._type = NodeType::BINARY_OPERATOR;
	else if( str=="*" ) node._function = []( const double * values ){ return values[0]*values[1]; }          , node._type = NodeType::BINARY_OPERATOR;
	else if( str=="/" ) node._function = []( const double * values ){ return values[0]/values[1]; }          , node._type = NodeType::BINARY_OPERATOR;
	else if( str=="^" ) node._function = []( const double * values ){ return pow( values[0] , values[1] ); } , node._type = NodeType::BINARY_OPERATOR;
	else MK_THROW( "Failed to parse function: " , str );
	return node;
}

inline Node Node::operator + ( const Node & node ) const { return _GetFunction( "+" , *this , node ); }
inline Node Node::operator - ( const Node & node ) const { return _GetFunction( "-" , *this , node ); }
inline Node Node::operator * ( const Node & node ) const { return _GetFunction( "*" , *this , node ); }
inline Node Node::operator / ( const Node & node ) const { return _GetFunction( "/" , *this , node ); }
inline Node Node::operator * ( double s ) const { return _GetFunction( "*" , *this , _GetConstant( s ) ); }
inline Node Node::operator / ( double s ) const { return _GetFunction( "*" , *this , _GetConstant( 1./s ) ); }

inline Node Node::_GetDConstant( double )
{
	Node node;
	node._type = NodeType::ZERO;
	return node;
}

inline Node Node::_GetDVariable( unsigned int vIndex , unsigned int dIndex )
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

inline Node Node::_GetDFunction( std::string fName, const Node & n , const Node & d )
{
	if     ( fName=="-"   ) return _GetFunction( "-" , d );
	else if( fName=="exp" ) return _GetFunction( "*" , d , _GetFunction( "exp" , n ) );
	else if( fName=="log" ) return _GetFunction( "*" , d , _GetFunction( "^" , n , _GetConstant( -1. ) ) );
	else if( fName=="cos" ) return _GetFunction( "*" , d , _GetFunction( "-" , _GetFunction( "sin" , n ) ) );
	else if( fName=="sin" ) return _GetFunction( "*" , d , _GetFunction( "cos" , n ) );
	// d(sin/cos) = d(sin) / cos + sin * d(1/cos)
	//            = cos / cos - sin * 1/cos^2 * d(cos)
	//            = 1 + sin^2 / cos^2
	//            = 1 + tan^2
	else if( fName=="tan" ) return _GetFunction( "*" , d , _GetFunction( "+" , _GetConstant( 1. ) , _GetFunction( "^" , _GetFunction( "tan" , n ) , _GetConstant( 2. ) ) ) );
	else if( fName=="cosh" ) return _GetFunction( "*" , d , _GetFunction( "sinh" , n ) );
	else if( fName=="sinh" ) return _GetFunction( "*" , d , _GetFunction( "cosh" , n ) );
	else if( fName=="tanh" ) return _GetFunction( "/" , d , _GetFunction( "^" , _GetFunction( "cosh" , n ) , _GetConstant( 2 ) ) );
	else MK_THROW( "Failed to parse function: " , fName );
	return Node();
}

inline Node Node::_GetDFunction( std::string fName , const Node & node1 , const Node & dNode1 , const Node & node2 , const Node & dNode2 )
{
	if     ( fName=="+" ) return _GetFunction( "+" , dNode1 , dNode2 );
	else if( fName=="-" ) return _GetFunction( "-" , dNode1 , dNode2 );
	else if( fName=="*" ) return _GetFunction( "+" , _GetFunction( "*" , dNode1 , node2 ) , _GetFunction( "*" , node1 , dNode2 ) );
	else if( fName=="/" ) return _GetFunction( "-" , _GetFunction( "/" , dNode1 , node2 ) , _GetFunction( "*" , dNode2 , _GetFunction( "/" , node1 , _GetFunction( "^" , node2 , _GetConstant( 2. ) ) ) ) );
	// d( f^g ) = d( e^( log f * g ) )
	//          = e^( log f * g ) * d( log f * g )
	//          = f^g * ( log f * d(g) + d( log f ) * g )
	//          = f^g * ( log f * d(g) + 1/f * d(f) * g )
	else if( fName=="^" )
	{
		if( node2._type==NodeType::CONSTANT ) return _GetFunction( "*" , _GetConstant( node2._value ) , _GetFunction( "*" , _GetFunction( "^" , node1 , _GetConstant( node2._value-1 ) ) , dNode1 ) );
		else return _GetFunction( "*" , _GetFunction( "^" , node1 , node2 ) , _GetFunction( "+" , _GetFunction( "*" , _GetFunction( "log" , node1 ) , dNode2 ) , _GetFunction( "*" , dNode1 , _GetFunction( "/" , node2 , node1 ) ) ) );
	}
	else MK_THROW( "Failed to parse function: " , fName );
	return Node();

}

inline Node Node::_Parse( _StateInfo stateInfo , const std::vector< std::string > & vars )
{
	Node node;
	if( stateInfo.state.size()==0 ) return node;
	else if( stateInfo.state.size()==1 )
	{
		if     ( stateInfo.state[0].type==_StateInfo::NodeType::CONSTANT ) return _GetConstant( stateInfo.state[0].name );
		else if( stateInfo.state[0].type==_StateInfo::NodeType::VARIABLE ) return _GetVariable( stateInfo.state[0].name , vars );
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

			node._type = _StateInfo::ConvertNodeType( stateInfo.state[idx].type );
			if( stateInfo.state[idx].type==_StateInfo::NodeType::FUNCTION || stateInfo.state[idx].type==_StateInfo::NodeType::UNARY_OPERATOR )
			{
				if( idx+1>stateInfo.state.size() ) MK_THROW( "Expected argument after function" );
				unsigned int begin = idx+1;

				Node childNode;
				if     ( stateInfo.state[begin].type==_StateInfo::NodeType::CONSTANT  ) childNode = _GetConstant( stateInfo.state[begin].name );
				else if( stateInfo.state[begin].type==_StateInfo::NodeType::VARIABLE  ) childNode = _GetVariable( stateInfo.state[begin].name , vars );
				else if( stateInfo.state[begin].type==_StateInfo::NodeType::L_PARENTH ) childNode = _Parse( stateInfo.sub( begin , stateInfo.closingParenth( begin )+1 ) , vars );
				else MK_THROW( "Expected constant, variable, or parentheses-enclosed equation" );
				return _GetFunction( stateInfo.state[idx].name , childNode );
			}
			else if( stateInfo.state[idx].type==_StateInfo::NodeType::BINARY_OPERATOR )
			{
				if( idx==0 ) MK_THROW( "Expected content before binary operator" );
				if( idx+1>stateInfo.state.size() ) MK_THROW( "Expected content after binary operator" );

				Node childNode1 , childNode2;
				{
					unsigned int end = idx-1;
					if     ( stateInfo.state[end].type==_StateInfo::NodeType::CONSTANT  ) childNode1 = _GetConstant( stateInfo.state[end].name );
					else if( stateInfo.state[end].type==_StateInfo::NodeType::VARIABLE  ) childNode1 = _GetVariable( stateInfo.state[end].name , vars );
					else if( stateInfo.state[end].type==_StateInfo::NodeType::R_PARENTH ) childNode1 = _Parse( stateInfo.sub( stateInfo.openingParenth( end ) , end+1 ) , vars );
					else MK_THROW( "Expected constant, variable, or parenthesized expression" );
				}
				{
					unsigned int begin = idx+1;
					if     ( stateInfo.state[begin].type==_StateInfo::NodeType::CONSTANT  ) childNode2 = _GetConstant( stateInfo.state[begin].name );
					else if( stateInfo.state[begin].type==_StateInfo::NodeType::VARIABLE  ) childNode2 = _GetVariable( stateInfo.state[begin].name, vars );
					else if( stateInfo.state[begin].type==_StateInfo::NodeType::L_PARENTH ) childNode2 = _Parse( stateInfo.sub( begin , stateInfo.closingParenth( begin )+1 ) , vars );
					else MK_THROW( "Expected constant, variable, or parenthesized expression" );
				}
				return _GetFunction( stateInfo.state[idx].name , childNode1 , childNode2 );
			}
			else MK_THROW( "Expected operator or function" );
		}
	}
	return node;
}

std::ostream & operator << ( std::ostream & stream , const Node & node )
{
	switch( node._type )
	{
	case Node::NodeType::ZERO:            return stream << 0;
	case Node::NodeType::CONSTANT:        return stream << node._value;
	case Node::NodeType::VARIABLE:        return stream << "x" << node._variableIndex;
	case Node::NodeType::UNARY_OPERATOR:  return stream << node._functionName << "(" << node._children[0] << ")";
	case Node::NodeType::BINARY_OPERATOR:
		if( node._children[0]._type==Node::NodeType::CONSTANT || node._children[0]._type==Node::NodeType::VARIABLE ) stream << node._children[0];
		else stream << "(" << node._children[0] << ")";
		stream << node._functionName;
		if( node._children[1]._type==Node::NodeType::CONSTANT || node._children[1]._type==Node::NodeType::VARIABLE ) return stream << node._children[1];
		else return stream << "(" << node._children[1] << ")";
	case Node::NodeType::FUNCTION:        return stream << node._functionName << "(" << node._children[0] << ")";
	default:
		MK_THROW( "Unrecognized type: " , static_cast< int >( node._type ) );
		return stream;
	}
}

inline bool Node::isConstant( void ) const
{
	bool isConstant = true;
	if( _children.size() ) for( unsigned int i=0 ; i<_children.size() ; i++ ) isConstant &= _children[i].isConstant();
	else isConstant &= _type!=NodeType::VARIABLE;
	return isConstant;
}

inline void Node::compress( void )
{
	for( unsigned int i=0 ; i<_children.size() ; i++ ) _children[i].compress();
#if 1 // NEW_CODE
	if( _type==NodeType::BINARY_OPERATOR && _functionName=="*" && ( _children[0]._type==NodeType::ZERO || _children[1]._type==NodeType::ZERO ) )
	{
		_type = NodeType::ZERO;
		_children.resize(0);
	}
#endif // NEW_CODE
	if( isConstant() )
	{
		_value = this->operator()( nullptr );
		_type = NodeType::CONSTANT;
#if 1 // NEW_CODE
		_children.resize(0);
#endif // NEW_CODE
	}
	if( _type==NodeType::CONSTANT && !_value ) _type = NodeType::ZERO;


	// Clean up addition/multiplication/exponentation by 0
	if( _type==NodeType::BINARY_OPERATOR )
	{
		if( _children[0]._type==NodeType::ZERO )
		{
			if( _functionName=="*" || _functionName=="/" || _functionName=="^" )
			{
				_children.clear();
				_type = NodeType::ZERO;
			}
			else if( _functionName=="+" )
			{
				_type          = _children[1]._type;
				_functionName  = _children[1]._functionName;
				_function      = _children[1]._function;
				_value         = _children[1]._value;
				_variableIndex = _children[1]._variableIndex;
				std::swap( _children , _children[1]._children );
			}
			else if( _functionName=="-" )
			{
				_type = NodeType::UNARY_OPERATOR;
				_functionName = "-";
				_function = []( const double * values ){ return -values[0]; };
				std::swap( _children[0] , _children[1] );
				_children.pop_back();
			}
		}
		else if( _children[1]._type==NodeType::ZERO )
		{
			if( _functionName=="*" )
			{
				_children.clear();
				_type = NodeType::ZERO;
			}
			else if( _functionName=="^" )
			{
				_children.clear();
				_type = NodeType::CONSTANT;
				_value = 1.;
			}
			else if( _functionName=="+" || _functionName=="-" )
			{
				_type          = _children[0]._type;
				_functionName  = _children[0]._functionName;
				_function      = _children[0]._function;
				_value         = _children[0]._value;
				_variableIndex = _children[0]._variableIndex;
				std::swap( _children , _children[0]._children );
			}
		}
	}

	// Clean up multiplication/exponentation by 1/-1
	if( _type==NodeType::BINARY_OPERATOR )
	{
		if( _children[0]._type==NodeType::CONSTANT && ( _children[0]._value==1 || _children[0]._value==-1 ) )
		{
			if( _functionName=="*" )
				if( _children[0]._value==1 )
				{
					_type          = _children[1]._type;
					_functionName  = _children[1]._functionName;
					_function      = _children[1]._function;
					_value         = _children[1]._value;
					_variableIndex = _children[1]._variableIndex;
					std::swap( _children , _children[1]._children );
				}
				else if( _children[0]._value==-1 )
				{
					_type = NodeType::UNARY_OPERATOR;
					_functionName = "-";
					_function = []( const double * values ){ return -values[0]; };
					std::swap( _children[0] , _children[1] );
					_children.pop_back();
				}
		}
		else if( _children[1]._type==NodeType::CONSTANT && ( _children[1]._value==1 || _children[1]._value==-1 ) )
		{
			if( _functionName=="*" || _functionName=="/" )
			{
				if( _children[1]._value==1 )
				{
					_type          = _children[0]._type;
					_functionName  = _children[0]._functionName;
					_function      = _children[0]._function;
					_value         = _children[0]._value;
					_variableIndex = _children[0]._variableIndex;
					std::swap( _children , _children[0]._children );
				}
				else if( _children[1]._value==-1 )
				{
					_type = NodeType::UNARY_OPERATOR;
					_functionName = "-";
					_function = []( const double * values ){ return -values[0]; };
					_children.pop_back();
				}
			}
			else if( _functionName=="^" && _children[1]._value==1 )
			{
				_type          = _children[0]._type;
				_functionName  = _children[0]._functionName;
				_function      = _children[0]._function;
				_value         = _children[0]._value;
				_variableIndex = _children[0]._variableIndex;
				std::swap( _children , _children[0]._children );
			}
		}
	}
}

template< unsigned int Dim >
double Node::operator()( Point< double , Dim > p ) const { return operator()( &p[0] ); }

inline double Node::operator()( const double * values ) const
{
	double _values[2];
	switch( _type )
	{
	case Node::NodeType::ZERO:
		return 0;
	case Node::NodeType::CONSTANT:
		return _value;
	case Node::NodeType::VARIABLE:
		return values[ _variableIndex ];
	case Node::NodeType::UNARY_OPERATOR:
	case Node::NodeType::FUNCTION:
		_values[0] = _children[0]( values );
		return _function( _values );
	case Node::NodeType::BINARY_OPERATOR:
		_values[0] = _children[0]( values );
		_values[1] = _children[1]( values );
		return _function( _values );
	default:
		MK_THROW( "Unrecognized type" );
		return 0.;
	}
}

inline Node Node::d( unsigned int dIndex ) const
{
	Node node;
	if     ( _type==NodeType::ZERO ) ;
	else if( _type==NodeType::CONSTANT ) node = _GetDConstant( _value );
	else if( _type==NodeType::VARIABLE ) node = _GetDVariable( _variableIndex , dIndex );
	else if( _type==NodeType::UNARY_OPERATOR || _type==NodeType::FUNCTION ) node = _GetDFunction( _functionName , _children[0] , _children[0].d( dIndex ) );
	else if( _type==NodeType::BINARY_OPERATOR ) node = _GetDFunction( _functionName , _children[0] , _children[0].d(dIndex) , _children[1] , _children[1].d(dIndex) );
	else MK_THROW( "Unrecognized type" );
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

//////////////////////
// Node::_StateInfo //
//////////////////////
inline Node::NodeType Node::_StateInfo::ConvertNodeType( NodeType type )
{
	Node::NodeType _type;
	switch( type )
	{
	case NodeType::VARIABLE:        _type = Node::NodeType::VARIABLE        ; break;
	case NodeType::CONSTANT:        _type = Node::NodeType::CONSTANT        ; break;
	case NodeType::UNARY_OPERATOR:  _type = Node::NodeType::UNARY_OPERATOR  ; break;
	case NodeType::BINARY_OPERATOR: _type = Node::NodeType::BINARY_OPERATOR ; break;
	case NodeType::FUNCTION:        _type = Node::NodeType::FUNCTION        ; break;
	default: MK_THROW( "Could not convert node type: " , NodeTypeNames[ static_cast< int >( type ) ] );
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
#if 1 // NEW_CODE
				if( tokens[i]=="Pi" ) state.emplace_back( NodeType::CONSTANT , tokens[i] );
				else
				{
					double foo = std::stod( tokens[i] );
					state.emplace_back( NodeType::CONSTANT , tokens[i] );
				}
#else // !NEW_CODE
				double foo = std::stod( tokens[i] );
				state.emplace_back( NodeType::CONSTANT , tokens[i] );
#endif // NEW_CODE
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
