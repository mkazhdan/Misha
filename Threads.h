/*
Copyright (c) 2019, Michael Kazhdan
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

#ifndef MISHA_THREADS_INCLUDED
#define MISHA_THREADS_INCLUDED

#include <thread>
#include <vector>

namespace Misha
{
	// The assumption is that Kernel is the type of a function taking two arguments, the first is the index of the thread and the second is index of the iteration
	struct Thread
	{
		static unsigned int Threads;

		template< typename Kernel , typename ... Args >
		static void parallel_for( size_t start , size_t end , Kernel kernel , Args ... args )
		{
			size_t range = end - start;

			auto _functor = [&]( unsigned int t )
			{
				size_t _start = start + ( range * t ) / Threads;
				size_t _end = start + ( range *(t+1) ) / Threads;
				for( size_t i=_start ; i<_end ; i++ ) kernel( t , i , args ... );
			};

			std::vector< std::thread > threads( Threads );

			for( unsigned int t=0 ; t<Threads ; t++ ) threads[t] = std::thread( [&]( unsigned int _t ){ _functor(_t); } , t );

			for( unsigned int t=0 ; t<Threads ; t++ ) threads[t].join();
		}
	};
	unsigned int Thread::Threads = std::thread::hardware_concurrency();
}
#endif // MISHA_THREADS_INCLUDED
