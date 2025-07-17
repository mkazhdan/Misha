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

inline std::string DataTypeName( DataType dataType )
{
	switch( dataType )
	{
	case Binary: return std::string( "bool" );
	case UChar:  return std::string( "unsigned char" );
	case Short:  return std::string( "short" );
	case UShort: return std::string( "unsigned short" );
	case Int:    return std::string( "int" );
	case Float:  return std::string( "float" );
	case Double: return std::string( "double" );
	default:     return std::string( "unknown" );
	}
}

inline unsigned int Dimension( std::string hdr_file )
{
	nifti_1_header hdr;
	FILE * fp;


	/********** open and read header */
	fp = fopen( hdr_file.c_str() , "rb" );
	if( fp==NULL )
	{
		fprintf( stderr , "\nError opening header file %s\n" , hdr_file.c_str() );
		exit(1);
	}
	if( fread( &hdr , MIN_HEADER_SIZE , 1 , fp )!=1 )
	{
		fprintf( stderr , "\nError reading header file %s\n" , hdr_file.c_str() );
		exit(1);
	}

	return static_cast< unsigned int >( hdr.dim[0] );
}

template< unsigned int Dim >
long SetHeader( std::string hdr_file , Header< Dim > & header )
{
	nifti_1_header hdr;
	FILE * fp;


	/********** open and read header */
	fp = fopen( hdr_file.c_str() , "rb" );
	if( fp==NULL )
	{
		fprintf( stderr , "\nError opening header file %s\n" , hdr_file.c_str() );
		exit(1);
	}
	if( fread( &hdr , MIN_HEADER_SIZE , 1 , fp )!=1 )
	{
		fprintf( stderr , "\nError reading header file %s\n" , hdr_file.c_str() );
		fclose( fp );
		exit(1);
	}

	if( hdr.dim[0]!=Dim )
	{
		fprintf( stderr , "\nError dimensions don't match %d != %d\n" , hdr.dim[0] , static_cast< int >( Dim ) );
		fclose( fp );
		exit(1);
	}

	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		header.res[d] = hdr.dim[d+1];
		header.scale[d] = hdr.pixdim[d+1];
	}

	switch( hdr.datatype )
	{
	case DT_UNKNOWN:       header.dataType = DataType::Unknown ; break;
	case DT_BINARY:        header.dataType = DataType::Binary  ; break;
	case DT_UNSIGNED_CHAR: header.dataType = DataType::UChar   ; break;
	case DT_SIGNED_SHORT:  header.dataType = DataType::Short   ; break;
	case DT_UINT16:        header.dataType = DataType::UShort  ; break;
	case DT_SIGNED_INT:    header.dataType = DataType::Int     ; break;
	case DT_FLOAT:         header.dataType = DataType::Float   ; break;
	case DT_DOUBLE:        header.dataType = DataType::Double  ; break;
	default:
		fprintf( stderr , "\nError parsing data-type %d\n" , hdr.datatype );
		exit(1);
	}


//	fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",hdr.scl_slope,hdr.scl_inter);

	fclose(fp);

	return static_cast< long >(hdr.vox_offset);
}

template< unsigned int Dim , DataType DT >
typename DataTypeInfo< DT >::Type * ReadData( std::string fileName , Header< Dim > & header )
{
	using Type = typename DataTypeInfo< DT >::Type;

	FILE * fp;
	long offset = SetHeader( fileName , header );
	Type * data;

	if( header.dataType!=DT )
	{
		fprintf( stderr , "\nData types don't match %d != %d\n" , header.dataType , DT );
		exit(1);
	}


	fp = fopen( fileName.c_str() , "rb" );
	if( fp==NULL )
	{
		fprintf( stderr , "\nError opening data file %s\n" , fileName.c_str() );
		exit(1);
	}

	if( fseek( fp, offset , SEEK_SET )!=0 )
	{
		fprintf(stderr, "\nError doing fseek() to %ld in data file %s\n" , offset , fileName.c_str() );
		fclose( fp );
		exit(1);
	}

	size_t sz = 1;
	for( unsigned int d=0 ; d<Dim ; d++ ) sz *= header.res[d];

	data = new Type[ sz ];
	if( !data )
	{
		fprintf( stderr , "\nError allocating data buffer for %s\n" , fileName.c_str() );
		fclose( fp );
		exit(1);
	}

	size_t _sz = fread( data , sizeof(Type) , sz , fp );
	if( _sz != sz )
	{
		fprintf( stderr , "\nError reading volume 1 from %s: %lld != %lld\n" , fileName.c_str() , sz , _sz );
		fclose( fp );
		exit(1);
	}

	fclose(fp);

	return data;
}
