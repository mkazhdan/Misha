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

#ifndef FACTOR_INCLUDED
#define FACTOR_INCLUDED

#include <math.h>

#define PI 3.1415926535897932384
#define SQRT_3 1.7320508075688772935

namespace MishaK
{
	double ArcTan2(const double& y,const double& x);
	double Angle(const double in[2]);
	void Sqrt(const double in[2],double out[2]);
	void Add(const double in1[2],const double in2[2],double out[2]);
	void Subtract(const double in1[2],const double in2[2],double out[2]);
	void Multiply(const double in1[2],const double in2[2],double out[2]);
	void Divide(const double in1[2],const double in2[2],double out[2]);

	int Factor(double a1,double a0,double roots[1][2],const double& EPS);
	int Factor(double a2,double a1,double a0,double roots[2][2],const double& EPS);
	int Factor(double a3,double a2,double a1,double a0,double roots[3][2],const double& EPS);
	int Factor(double a4,double a3,double a2,double a1,double a0,double roots[4][2],const double& EPS);

	int Solve(const double* eqns,const double* values,double* solutions,const int& dim);
#include "Factor.inl"
}
#endif // FACTOR_INCLUDED
