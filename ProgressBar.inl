/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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
// ProgressBar //
/////////////////
inline double ProgressBar::Time( void )
{
#ifdef WIN32
	struct _timeb t;
	_ftime64_s(&t);
	return double(t.time)+double(t.millitm)/1000.0;
#else // WIN32
	struct timeval t;
	gettimeofday(&t,NULL);
	return t.tv_sec+(double)t.tv_usec/1000000;
#endif // WIN32
}

inline ProgressBar::ProgressBar( int bins , size_t total , const char* header , bool outputOnDestruction ) : _outputOnDestruction(outputOnDestruction)
{
	_startTime = Time();
	_bins = bins;
	_total = total;
	_header = header;
	_idx = 0;
	_previousTime = -1;
}

inline void ProgressBar::update( bool output )
{
	if( output ) print( );
	AddAtomic( _idx , (size_t)1 );
}
inline void ProgressBar::print( void )
{
	int currentBin = int( (_idx*_bins) / (_total-1 ) );
	double currentTime = Time() - _startTime;
	if( int( currentTime*10 ) != int( _previousTime*10 ) )
	{
		printf( "\r[" );
		for( int i=0 ; i<currentBin ; i++ ) printf( "." );
		for( int i=currentBin ; i<_bins ; i++ ) printf( " " );
		printf( "] %s: %.1f (s)\r" , _header , currentTime );
		_previousTime = currentTime;
	}
}
inline ProgressBar::~ProgressBar( void )
{
	if( _outputOnDestruction )
	{
		double currentTime = Time() - _startTime;
		printf( "[" );
		for( int i=0 ; i<_bins ; i++ ) printf( "." );
		printf( "] %s: %.1f (s)\n" , _header , currentTime );
	}
}