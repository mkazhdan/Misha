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


#include <vector>
#include <string>
#include <functional>
#include <Misha/Exceptions.h>
#include <Misha/Geometry.h>

#define NEW_EQUATION_PARSER

namespace MishaK
{
	namespace EquationParser
	{
		// Interior nodes store functions/operators and the nodes to which they are applied
		// Leaf nodes store variables/constants
		struct Node
		{
			enum struct NodeType
			{
				ZERO ,
				CONSTANT ,
				VARIABLE ,
				UNARY_OPERATOR ,
				BINARY_OPERATOR ,
				FUNCTION
			};

			static inline const std::string NodeTypeNames[] =
			{
				"zero" ,
				"constant" ,
				"variable" ,
				"operator (unary)" ,
				"operator (binary)" ,
				"function" ,
			};

			struct State
			{
				NodeType type;
				std::string name;

				State( NodeType type=NodeType::CONSTANT , std::string name="" );
			};


			Node( void );

			static Node Parse( std::string eqn , const std::vector< std::string > & vars );

			double operator()( const double * values ) const;
			template< unsigned int Dim >
			double operator()( Point< double , Dim > p ) const;
			Node d( unsigned int dIndex ) const;

			Node operator + ( const Node & n ) const;
			Node operator - ( const Node & n ) const;
			Node operator * ( const Node & n ) const;
			Node operator / ( const Node & n ) const;
			Node operator * ( double s ) const;
			Node operator / ( double s ) const;

			bool isConstant( void ) const;
			void compress( void );

			unsigned int size( void ) const;

			friend std::ostream & operator << ( std::ostream & stream , const Node & node );

		protected:
			std::vector< Node > _children;

			NodeType _type;
			// These should be a union...
			std::string _functionName;
			std::function< double ( const double * ) > _function;
			double _value;
			unsigned int _variableIndex;

			struct _StateInfo
			{
				friend struct Node;

				enum struct NodeType
				{
					VARIABLE ,
					CONSTANT ,
					UNARY_OPERATOR ,
					BINARY_OPERATOR ,
					FUNCTION ,
					UNKNOWN ,
					L_PARENTH ,
					R_PARENTH
				};

				static inline const std::string NodeTypeNames[] =
				{
					"variable" ,
					"constant" ,
					"operator (unary)" ,
					"operator (binary)" ,
					"function" ,
					"unknown" ,
					"parentheses (left)" ,
					"parentheses (right)"
				};

				static Node::NodeType ConvertNodeType( NodeType type );

				struct State
				{
					NodeType type;
					std::string name;

					State( NodeType type=NodeType::UNKNOWN , std::string name="" );
				};

				inline static const std::vector< std::vector< std::string > > UnaryOperators = { { "-" } };
				inline static const std::vector< std::vector< std::string > > BinaryOperators = { { "^" } , { "*" , "/" } , { "+" , "-" } };

				std::vector< State > state;

				_StateInfo( void );
				_StateInfo( std::string eqn , const std::vector< std::string > & vars );

				unsigned int openingParenth( unsigned int idx ) const;
				unsigned int closingParenth( unsigned int idx ) const;

				_StateInfo sub( unsigned int begin , unsigned end ) const;

				void addFunctionParenths( void );
				void addUnaryOperatorParenths( const std::vector< std::string > & ops );
				void addBinaryOperatorParenths( const std::vector< std::string > & ops );

				static bool IsParenth( char c );
				static bool IsOperator( char c );
				static bool IsOperator( std::string str );
				static bool IsOperator( char c , const std::vector< std::string > & ops );
				static bool IsOperator( std::string str , const std::vector< std::string > & ops );
				static bool IsUnaryOperator( char c );
				static bool IsUnaryOperator( std::string str );
				static bool IsBinaryOperator( char c );
				static bool IsBinaryOperator( std::string str );
			};

			static Node _GetConstant( double value );
			static Node _GetConstant( std::string str );
			static Node _GetVariable( std::string str , const std::vector< std::string > & vars );
			static Node _GetFunction( std::string , const Node & node );
			static Node _GetFunction( std::string , const Node & node1 , const Node & node2 );

			static Node _GetDConstant( double );
			static Node _GetDVariable( unsigned int vIndex , unsigned int dIndex );
			static Node _GetDFunction( std::string fName , const Node & node , const Node & dNode );
			static Node _GetDFunction( std::string fName  , const Node & node1 , const Node & dNode1 , const Node & node2 , const Node & dNode2 );

			static Node _Parse( _StateInfo stateInfo , const std::vector< std::string > & vars );
		};

		Node operator * ( double s , const Node & n ) { return n * s; }

#include "EquationParser.inl"
	}
}