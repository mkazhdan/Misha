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
#ifndef NIFTI_INCLUDED
#define NIFTI_INCLUDED

// Code modified from https://github.com/NIFTI-Imaging/nifti_clib/blob/master/real_easy/stand_alone_app/nifti1_read_write.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include "NIFTI/nifti1.h"

namespace MishaK
{
	namespace NIFTI
	{
		static const unsigned int MIN_HEADER_SIZE = 348;
		static const unsigned int NII_HEADER_SIZE = 352;

		enum DataType
		{
			Unknown ,
			Binary ,
			UChar ,
			UShort ,
			Short ,
			Int ,
			Float ,
			Double
		};

		std::string DataTypeName( DataType dataType );

		template< unsigned int > struct DataTypeInfo;

		template<> struct DataTypeInfo< DataType::Unknown >{ unsigned int niftiType = DT_UNKNOWN       ; };
		template<> struct DataTypeInfo< DataType::Binary  >{ unsigned int niftiType = DT_BINARY        ; using Type = bool; };
		template<> struct DataTypeInfo< DataType::UChar   >{ unsigned int niftiType = DT_UNSIGNED_CHAR ; using Type = unsigned char; };
		template<> struct DataTypeInfo< DataType::UShort  >{ unsigned int niftiType = DT_UINT16        ; using Type = unsigned short; };
		template<> struct DataTypeInfo< DataType::Short   >{ unsigned int niftiType = DT_SIGNED_SHORT  ; using Type = short; };
		template<> struct DataTypeInfo< DataType::Int     >{ unsigned int niftiType = DT_SIGNED_INT    ; using Type = int; };
		template<> struct DataTypeInfo< DataType::Float   >{ unsigned int niftiType = DT_FLOAT         ; using Type = float; };
		template<> struct DataTypeInfo< DataType::Double  >{ unsigned int niftiType = DT_DOUBLE        ; using Type = double; };

		template< unsigned int Dim >
		struct Header
		{
			unsigned int res[Dim];
			DataType dataType;
			float scale[Dim];
		};

		unsigned int Dimension( std::string fileName );

		template< unsigned int Dim >
		long SetHeader( std::string fileName , Header< Dim > & header );

		template< unsigned int Dim , DataType DT >
		typename DataTypeInfo< DT >::Type * ReadData( std::string fileName , Header< Dim > & header );

#include "NIFTI.inl"
	};
}
#endif // NIFTI_INCLUDED
