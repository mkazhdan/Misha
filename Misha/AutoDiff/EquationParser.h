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
				NEGATION ,
				ADDITION ,
				SUBTRACTION ,
				MULTIPLICATION ,
				DIVISION ,
				POWER ,
				EXPONENTIAL ,
				NATURAL_LOGARITHM ,
				COSINE ,
				SINE ,
				TANGENT ,
				HYPERBOLIC_COSINE ,
				HYPERBOLIC_SINE ,
				HYPERBOLIC_TANGENT
			};

			static inline const std::string NodeTypeNames[] =
			{
				"zero" ,
				"constant" ,
				"variable" ,
				"negation" ,
				"addition" ,
				"subtraction" ,
				"multiplication" ,
				"division" ,
				"power" ,
				"exponential" ,
				"logarithm" ,
				"cosine" ,
				"sine" ,
				"tangent" ,
				"hyperbolic cosine" ,
				"hyperbolic sine" ,
				"hyperbolic tangent"
			};

			static bool IsUnaryOperator( NodeType type ){ return type==NodeType::NEGATION; }
			static bool IsBinaryOperator( NodeType type )
			{
				switch( type )
				{
				case NodeType::ADDITION:
				case NodeType::SUBTRACTION:
				case NodeType::MULTIPLICATION:
				case NodeType::DIVISION:
				case NodeType::POWER:
					return true;
				default: return false;
				}
			}
			static bool IsFunction( NodeType type )
			{
				switch( type )
				{
				case NodeType::EXPONENTIAL:
				case NodeType::COSINE:
				case NodeType::SINE:
				case NodeType::TANGENT:
				case NodeType::HYPERBOLIC_COSINE:
				case NodeType::HYPERBOLIC_SINE:
				case NodeType::HYPERBOLIC_TANGENT:
					return true;
				default: return false;
				}
			}

			static std::string ToString( NodeType type )
			{
				switch( type )
				{
				case NodeType::NEGATION:           return "-";
				case NodeType::ADDITION:           return "+";
				case NodeType::SUBTRACTION:        return "-";
				case NodeType::MULTIPLICATION:     return "*";
				case NodeType::DIVISION:           return "/";
				case NodeType::POWER:              return "^";
				case NodeType::EXPONENTIAL:        return "exp";
				case NodeType::COSINE:             return "cos";
				case NodeType::SINE:               return "sin";
				case NodeType::TANGENT:            return "tan";
				case NodeType::HYPERBOLIC_COSINE:  return "cosh";
				case NodeType::HYPERBOLIC_SINE:    return "sinh";
				case NodeType::HYPERBOLIC_TANGENT: return "tanh";
				default:
					MK_THROW( "Unreognized node type: " , NodeTypeNames[ static_cast< unsigned int >(type) ] );
					return "";
				}
			}

			struct State
			{
				NodeType type;
				std::string name;

				State( NodeType type=NodeType::CONSTANT , std::string name="" );
			};


			Node( void );
			Node( double c );
			static Node Variable( unsigned int idx );
			static Node DVariable( unsigned int idx , unsigned int dIdx );
			static Node Function( NodeType type , const Node & node );
			static Node Function( NodeType type , const Node & node1 , const Node & node2 );
			static Node DFunction( NodeType type , const Node & node , const Node & dNode );
			static Node DFunction( NodeType type , const Node & node1 , const Node & dNode1 , const Node & node2 , const Node & dNode2 );

			static Node Parse( std::string eqn , const std::vector< std::string > & vars );

			double operator()( const double * values ) const;
			template< unsigned int Dim >
			double operator()( Point< double , Dim > p ) const;
			Node d( unsigned int dIndex ) const;

			bool isConstant( void ) const;
			bool isProduct( void ) const;
			unsigned int hasVariable( unsigned int idx ) const;
			void compress( void );

			unsigned int size( void ) const;

			std::string operator()( const std::vector< std::string > &varNames ) const;

			friend std::ostream & operator << ( std::ostream & stream , const Node & node );

		protected:
			std::vector< Node > _children;

			NodeType _type;
			// These should be a union...
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

				struct State
				{
					NodeType type;
					std::string name;

					State( NodeType type=NodeType::UNKNOWN , std::string name="" );
					explicit operator Node::NodeType( void ) const;
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

			static Node _Parse( _StateInfo stateInfo , const std::vector< std::string > & vars );

			bool _compress( void );
			bool __compress( void );
			void _sanityCheck( void );
			void _insert( std::ostream & stream , const std::function< std::string ( unsigned int ) > & varName ) const;

			friend Node & operator += ( Node & n , const Node & _n );
			friend Node & operator -= ( Node & n , const Node & _n );
			friend Node & operator *= ( Node & n , const Node & _n );
			friend Node & operator /= ( Node & n , const Node & _n );
			friend Node & operator += ( Node & n , double s );
			friend Node & operator -= ( Node & n , double s );
			friend Node & operator *= ( Node & n , double s );
			friend Node & operator /= ( Node & n , double s );
			friend Node operator - ( const Node & n );
			friend Node operator + ( const Node & n1 , const Node & n2 );
			friend Node operator + ( const Node & n , double s );
			friend Node operator + ( double s , const Node & n );
			friend Node operator - ( const Node & n1 , const Node & n2 );
			friend Node operator - ( const Node & n , double s );
			friend Node operator - ( double s , const Node & n );
			friend Node operator * ( const Node & n1 , const Node & n2 );
			friend Node operator * ( const Node & n , double s );
			friend Node operator * ( double s , const Node & n );
			friend Node operator / ( const Node & n1 , const Node & n2 );
			friend Node operator / ( const Node & n , double s );
			friend Node operator / ( double s , const Node & n );
			friend Node operator ^ ( const Node & n1 , const Node & n2 );
			friend Node operator ^ ( const Node & n , double s );
			friend Node operator ^ ( double s , const Node & n );
			friend Node Exp( const Node & n );
			friend Node Pow( const Node & n1 , const Node & n2 );
			friend Node Pow( const Node & n , double s );
			friend Node Pow( double s , const Node & n );
			friend Node Log( const Node & n );
			friend Node Cos( const Node & n );
			friend Node Sin( const Node & n );
			friend Node Tan( const Node & n );
			friend Node Cosh( const Node & n );
			friend Node Sinh( const Node & n );
			friend Node Tanh( const Node & n );
			friend Node Sqrt( const Node & n );
		};

		Node & operator += ( Node & n , const Node & _n );
		Node & operator -= ( Node & n , const Node & _n );
		Node & operator *= ( Node & n , const Node & _n );
		Node & operator /= ( Node & n , const Node & _n );
		Node & operator += ( Node & n , double s );
		Node & operator -= ( Node & n , double s );
		Node & operator *= ( Node & n , double s );
		Node & operator /= ( Node & n , double s );

		Node operator - ( const Node & n );
		Node operator + ( const Node & n1 , const Node & n2 );
		Node operator + ( const Node & n , double s );
		Node operator + ( double s , const Node & n );
		Node operator - ( const Node & n1 , const Node & n2 );
		Node operator - ( const Node & n , double s );
		Node operator - ( double s , const Node & n );
		Node operator * ( const Node & n1 , const Node & n2 );
		Node operator * ( const Node & n , double s );
		Node operator * ( double s , const Node & n );
		Node operator / ( const Node & n1 , const Node & n2 );
		Node operator / ( const Node & n , double s );
		Node operator / ( double s , const Node & n );
		Node operator ^ ( const Node & n1 , const Node & n2 );
		Node operator ^ ( const Node & n , double s );
		Node operator ^ ( double s , const Node & n );
		Node Exp( const Node & n );
		Node Pow( const Node & n1 , const Node & n2 );
		Node Pow( const Node & n , double s );
		Node Pow( double s , const Node & n );
		Node Log( const Node & n );
		Node Cos( const Node & n );
		Node Sin( const Node & n );
		Node Tan( const Node & n );
		Node Cosh( const Node & n );
		Node Sinh( const Node & n );
		Node Tanh( const Node & n );
		Node Sqrt( const Node & n );

#include "EquationParser.inl"
	}
}