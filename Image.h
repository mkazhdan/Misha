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

#ifndef IMAGE_INCLUDED
#define IMAGE_INCLUDED

#if defined( NEW_CMD_LINE_PARSER ) || 1
#include <string.h>
#endif // NEW_CMD_LINE_PARSER
#include "Miscellany.h"

namespace MishaK
{
	struct ImageReader
	{
		virtual unsigned int nextRow( unsigned char* row ) = 0;
		static unsigned char* Read( std::string fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
		{
			ImageReader* reader = Get( fileName , width , height , channels );
			unsigned char* pixels = new unsigned char[ width*height*channels ];
			for( unsigned int j=0 ; j<height ; j++ ) reader->nextRow( pixels + j*width*channels );
			delete reader;
			return pixels;
		}
		static unsigned char* ReadColor( std::string fileName , unsigned int& width , unsigned int& height )
		{
			unsigned int channels;
			ImageReader* reader = Get( fileName , width , height , channels );
			if( channels!=1 && channels!=3 ) ERROR_OUT( "Only one- or three-channel input supported: " , channels ) , exit( 0 );
			unsigned char* pixels = new unsigned char[ width*height*3 ];
			unsigned char* pixelRow = new unsigned char[ width*channels];
			for( unsigned int j=0 ; j<height ; j++ )
			{
				reader->nextRow( pixelRow );
				if     ( channels==3 ) memcpy( pixels+j*width*3 , pixelRow , sizeof(unsigned char)*width*3 );
				else if( channels==1 ) for( unsigned int i=0 ; i<width ; i++ ) for( unsigned int c=0 ; c<3 ; c++ ) pixels[j*width*3+i*3+c] = pixelRow[i];
			}
			delete[] pixelRow;
			delete reader;
			return pixels;
		}

		static ImageReader* Get( std::string fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
		static bool ValidExtension( std::string ext );
		static void GetInfo( std::string fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
		virtual ~ImageReader( void ){ }

	};

	struct ImageWriter
	{
		virtual unsigned int nextRow( const unsigned char* row ) = 0;
		virtual unsigned int nextRows( const unsigned char* rows , unsigned int rowNum ) = 0;
		static bool Write( std::string fileName , const unsigned char* pixels , unsigned int width , unsigned int height , int channels , int quality=100 )
		{
			ImageWriter* writer = Get( fileName , width , height , channels , quality );
			if( !writer ) return false;
			for( unsigned int j=0 ; j<height ; j++ ) writer->nextRow( pixels + j*width*channels );
			delete writer;
			return true;
		}
		static ImageWriter* Get( std::string fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality=100 );
		virtual ~ImageWriter( void ){ }
	};
}

#include "PNG.h"
#include "JPEG.h"
#include "BMP.h"
#include "PBM.h"
#include "CmdLineParser.h"


namespace MishaK
{

	inline std::string GetImageExtension( std::string imageName ){ return GetFileExtension( imageName ); }

	inline bool ImageReader::ValidExtension( std::string ext )
	{
		for( unsigned int i=0 ; i<ext.size() ; i++ ) ext[i] = std::tolower( ext[i] );

		if     ( ext==std::string( "jpeg" ) || ext==std::string( "jpg" ) ) return true;
		else if( ext==std::string( "png" )                               ) return true;
		else if( ext==std::string( "igrid" )                             ) return true;
		return false;
	}

	inline ImageReader* ImageReader::Get( std::string fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
	{
		ImageReader *reader = NULL;
		std::string ext = GetImageExtension( fileName );
		if     ( ext==std::string( "jpeg" ) || ext==std::string( "jpg" ) ) reader = new JPEGReader( fileName , width , height , channels );
		else if( ext==std::string( "png" )                               ) reader = new  PNGReader( fileName , width , height , channels );
		else if( ext==std::string( "bmp" )                               ) reader = new  BMPReader( fileName , width , height , channels );
		else if( ext==std::string( "pbm" )                               ) reader = new  PBMReader( fileName , width , height , channels );

		return reader;
	}
	inline void ImageReader::GetInfo( std::string fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
	{
		std::string ext = GetImageExtension( fileName );
		if     ( ext==std::string( "jpeg" ) || ext==std::string( "jpg" ) ) JPEGReader::GetInfo( fileName , width , height , channels );
		else if( ext==std::string( "png" )                               )  PNGReader::GetInfo( fileName , width , height , channels );
		else if( ext==std::string( "bmp" )                               )  BMPReader::GetInfo( fileName , width , height , channels );
		else if( ext==std::string( "pbm" )                               )  PBMReader::GetInfo( fileName , width , height , channels );
	}
	inline ImageWriter* ImageWriter::Get( std::string fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality )
	{
		ImageWriter* writer = NULL;
		std::string ext = GetImageExtension( fileName );
		if     ( ext==std::string( "jpeg" ) || ext==std::string( "jpg" ) ) writer = new JPEGWriter( fileName , width , height , channels , quality );
		else if( ext==std::string( "png" )                               ) writer = new  PNGWriter( fileName , width , height , channels , quality );
		else if( ext==std::string( "bmp" )                               ) writer = new  BMPWriter( fileName , width , height , channels , quality );
		else if( ext==std::string( "pbm" )                               ) writer = new  PBMWriter( fileName , width , height , channels , quality );

		return writer;
	}
}
#endif // IMAGE_INCLUDED
