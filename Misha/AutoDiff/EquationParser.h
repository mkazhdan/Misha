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
#include <optional>
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
			// Constructs an equation-tree evaluating to zero
			Node( void );

			// Constructs an equation-tree evaluating to the idx-th value
			Node( unsigned int idx );

			// Constructs an equation-tree for the specified equation
			Node( std::string eqn , const std::vector< std::string > & vars );

			// Evaluates the the equation-tree at a prescribed set of values
			double operator()( const double * values ) const;
			template< unsigned int Dim > double operator()( Point< double , Dim > p ) const;

			// Returns of the equation-tree the derivative with respect to the dIndex-th variable
			Node d( unsigned int dIndex ) const;

			// Attempts to simplify the equation-tree
			void compress( void );

			// Returns the number of nodes in the equation-tree
			unsigned int size( void ) const;

			// Converts the equation-tree to an equation string (with the prescribed variable names)
			std::string operator()( const std::vector< std::string > &varNames ) const;

			// Writes the equation-tree, as an equation, to the stream
			friend std::ostream & operator << ( std::ostream & stream , const Node & node );

			// Operators for comparing two equation-trees
			bool operator <  ( const Node & n ) const;
			bool operator == ( const Node & n ) const;
			bool operator != ( const Node & n ) const;

			// Confirms that the compressed equation-tree and the input equation-tree give the same values when evaluated on random input
			static double SanityCheckCompression( const Node & node , unsigned int evaluationCount=1 , double evaluationRadius=1. );

			// Confirms that converting the node to a string and parsing back to an equation-tree reproduces the input
			static void SanityCheckParsing( const Node & node );

			// Functionality for stepping/printing along a path through the equation-tree
			static void Trace( const Node & node );

		protected:
			enum struct _NodeType;
			struct _StateInfo;

			std::vector< Node > _children;
			_NodeType _type;
			// These should be a union...
			std::function< double ( const std::vector< double > & ) > _function;
			double _value;
			unsigned int _variableIndex;


			enum struct _NodeType
			{
				// Constant terminals
				ZERO ,
				CONSTANT ,
				// Variable terminals
				VARIABLE ,
				// N-ary
				ADDITION ,
				MULTIPLICATION ,
				// Binary
				POWER ,
				// Functions
				EXPONENTIAL ,
				NATURAL_LOGARITHM ,
				COSINE ,
				SINE ,
				TANGENT ,
				HYPERBOLIC_COSINE ,
				HYPERBOLIC_SINE ,
				HYPERBOLIC_TANGENT
			};

			static inline const std::string _NodeTypeNames[] =
			{
				"zero" ,
				"constant" ,
				"variable" ,
				"addition" ,
				"multiplication" ,
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

			static bool _IsTerminal( _NodeType type );
			static bool _IsBinaryOperator( _NodeType type );
			static bool _IsNAryOperator( _NodeType type );
			static bool _IsSymmetricOperator( _NodeType type );
			static bool _IsOperator( _NodeType type );
			static bool _IsFunction( _NodeType type );
			static std::string _ToString( _NodeType type );

			static Node _Constant( double c );
			static Node _Function( _NodeType type , const Node & node );
			static Node _Function( _NodeType type , const Node & node1 , const Node & node2 );
			static Node _Function( _NodeType type , const std::vector< Node > & nodes );
			static Node _DVariable( unsigned int idx , unsigned int dIdx );
			static Node _DFunction( _NodeType type , const Node & node , const Node & dNode );
			static Node _DFunction( _NodeType type , const Node & node1 , const Node & dNode1 , const Node & node2 , const Node & dNode2 );
			static Node _DFunction( _NodeType type , const std::vector< Node > & nodes , const std::vector< Node > & dNodes );

			static Node _Parse( _StateInfo stateInfo , const std::vector< std::string > & vars );

			bool _preCompress( void );
			bool __preCompress( void );
			void _sort( void );
			bool _postCompress( void );
			bool __postCompress( void );
			void _sanityCheck( void ) const;
			unsigned int _maxVarIndex( void ) const;
			bool _isDivisible( const Node & node ) const;
			bool _divide( const Node & node );
			bool _isNegative( void ) const;
			bool _isReciprocal( void ) const;
			bool _isConstant( void ) const;
			bool _hasNumerator( const Node & node ) const;
			bool _hasDenominator( const Node & node ) const;
			bool _removeNumerator( const Node & node );
			bool _removeDenominator( const Node & node );
			static void _Insert( std::ostream & stream , const Node & node , const std::function< std::string ( unsigned int ) > & varName , bool processSign );
			static void __Insert( std::ostream & stream , const Node & node , const std::function< std::string ( unsigned int ) > & varName );


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
			friend Node Pow( const Node & n1 , const Node & n2 );
			friend Node Pow( const Node & n , double s );
			friend Node Pow( double s , const Node & n );
			friend Node Exp( const Node & n );
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
		Node Pow( const Node & n1 , const Node & n2 );
		Node Pow( const Node & n , double s );
		Node Pow( double s , const Node & n );
		Node Exp( const Node & n );
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