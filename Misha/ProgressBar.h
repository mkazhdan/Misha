/*
Copyright (c) 2009, Michael Kazhdan and Ming Chuang
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
#ifndef PROGRESS_BAR_INCLUDED
#define PROGRESS_BAR_INCLUDED

#if 1 // NEW_CODE
#else // !NEW_CODE
#include <sys/timeb.h>
#ifndef WIN32
#include <sys/time.h>
#endif // WIN32
#endif // NEW_CODE
#include "MultiThreading.h"
#include "Atomic.h"
#include "Miscellany.h"

namespace MishaK
{
	class ProgressBar
	{
		int _bins;
		size_t _total;
		size_t _idx;
		const char* _header;
		double _startTime , _previousTime;
		bool _outputOnDestruction , _showPeakMemory;
	public:
		ProgressBar( int bins , size_t total , const char *header , bool outputOnDestruction=true , bool showPeakMemory=false );
		~ProgressBar( void );
		void update( bool output=true );
		void print( void );

#if 1 // NEW_CODE
#else // !NEW_CODE
		static double Time( void );
#endif // NEW_CODE
	};

#include "ProgressBar.inl"
}
#endif // PROGRESS_BAR_INCLUDED
