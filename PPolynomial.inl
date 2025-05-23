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

////////////////////////
// StartingPolynomial //
////////////////////////
template<int Degree>
template<int Degree2>
StartingPolynomial<Degree+Degree2> StartingPolynomial<Degree>::operator * (const StartingPolynomial<Degree2>& p) const{
	StartingPolynomial<Degree+Degree2> sp;
	if(start>p.start){sp.start=start;}
	else{sp.start=p.start;}
	sp.p=this->p*p.p;
	return sp;
}
template<int Degree>
StartingPolynomial<Degree> StartingPolynomial<Degree>::scale(const double& s) const{
	StartingPolynomial q;
	q.start=start*s;
	q.p=p.scale(s);
	return q;
}
template<int Degree>
StartingPolynomial<Degree> StartingPolynomial<Degree>::shift(const double& s) const{
	StartingPolynomial q;
	q.start=start+s;
	q.p=p.shift(s);
	return q;
}


template<int Degree>
int StartingPolynomial<Degree>::operator < (const StartingPolynomial<Degree>& sp) const{
	if(start<sp.start){return 1;}
	else{return 0;}
}
template<int Degree>
int StartingPolynomial<Degree>::Compare(const void* v1,const void* v2){
	double d=((StartingPolynomial*)(v1))->start-((StartingPolynomial*)(v2))->start;
	if		(d<0)	{return -1;}
	else if	(d>0)	{return  1;}
	else			{return  0;}
}

/////////////////
// PPolynomial //
/////////////////
template<int Degree>
PPolynomial<Degree>::PPolynomial(void){
	polyCount=0;
	polys=NULL;
}
template<int Degree>
PPolynomial<Degree>::PPolynomial(const PPolynomial<Degree>& p){
	polyCount=0;
	polys=NULL;
	set(p.polyCount);
	memcpy(polys,p.polys,sizeof(StartingPolynomial<Degree>)*p.polyCount);
}

template<int Degree>
PPolynomial<Degree>::~PPolynomial(void){
	if(polyCount){free(polys);}
	polyCount=0;
	polys=NULL;
}
template<int Degree>
void PPolynomial<Degree>::set(StartingPolynomial<Degree>* sps,const int& count){
	int i,c=0;
	set(count);
	qsort(sps,count,sizeof(StartingPolynomial<Degree>),StartingPolynomial<Degree>::Compare);
	for(i=0;i<count;i++){
		if(!c || sps[i].start!=polys[c-1].start){polys[c++]=sps[i];}
		else{polys[c-1].p+=sps[i].p;}
	}
	reset(c);
}
template <int Degree>
int PPolynomial<Degree>::size(void) const{return int(sizeof(StartingPolynomial<Degree>)*polyCount);}

template<int Degree>
void PPolynomial<Degree>::set(const size_t &size){
	if(polyCount){free(polys);}
	polyCount=0;
	polys=NULL;
	polyCount=size;
	if(size){
		polys=(StartingPolynomial<Degree>*)malloc(sizeof(StartingPolynomial<Degree>)*size);
		memset(polys,0,sizeof(StartingPolynomial<Degree>)*size);
	}
}
template<int Degree>
void PPolynomial<Degree>::reset(const size_t& newSize){
	polyCount=newSize;
	polys=(StartingPolynomial<Degree>*)realloc(polys,sizeof(StartingPolynomial<Degree>)*newSize);
}

template<int Degree>
PPolynomial<Degree>& PPolynomial<Degree>::operator = (const PPolynomial<Degree>& p){
	set(p.polyCount);
	memcpy(polys,p.polys,sizeof(StartingPolynomial<Degree>)*p.polyCount);
	return *this;
}

template<int Degree>
template<int Degree2>
PPolynomial<Degree>& PPolynomial<Degree>::operator  = (const PPolynomial<Degree2>& p){
	set(p.polyCount);
	for(int i=0;i<int(polyCount);i++){
		polys[i].start=p.polys[i].start;
		polys[i].p=p.polys[i].p;
	}
	return *this;
}

template<int Degree>
double PPolynomial<Degree>::operator ()(const double& t) const{
	double v=0;
	for(int i=0;i<int(polyCount) && t>polys[i].start;i++) v+=polys[i].p(t);
	return v;
}

template<int Degree>
double PPolynomial<Degree>::integral(const double& tMin,const double& tMax) const{
	int m=1;
	double start,end,s,v=0;
	start=tMin;
	end=tMax;
	if(tMin>tMax){
		m=-1;
		start=tMax;
		end=tMin;
	}
	for(int i=0;i<int(polyCount) && polys[i].start<end;i++){
		if(start<polys[i].start){s=polys[i].start;}
		else{s=start;}
		v+=polys[i].p.integral(s,end);
	}
	return v*m;
}
template<int Degree>
double PPolynomial<Degree>::Integral(void) const{return integral(polys[0].start,polys[polyCount-1].start);}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator + (const PPolynomial<Degree>& p) const{
	PPolynomial q;
	int i,j;
	size_t idx=0;
	q.set(polyCount+p.polyCount);
	i=j=-1;

	while(idx<q.polyCount){
		if		(j>=int(p.polyCount)-1)				{q.polys[idx]=  polys[++i];}
		else if	(i>=int(  polyCount)-1)				{q.polys[idx]=p.polys[++j];}
		else if(polys[i+1].start<p.polys[j+1].start){q.polys[idx]=  polys[++i];}
		else										{q.polys[idx]=p.polys[++j];}
//		if(idx && polys[idx].start==polys[idx-1].start)	{polys[idx-1].p+=polys[idx].p;}
//		else{idx++;}
		idx++;
	}
	return q;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator - (const PPolynomial<Degree>& p) const{
	PPolynomial q;
	int i,j;
	size_t idx=0;
	q.set(polyCount+p.polyCount);
	i=j=-1;

	while(idx<q.polyCount){
		if		(j>=int(p.polyCount)-1)				{q.polys[idx]=  polys[++i];}
		else if	(i>=int(  polyCount)-1)				{q.polys[idx].start=p.polys[++j].start;q.polys[idx].p=p.polys[j].p*(-1.0);}
		else if(polys[i+1].start<p.polys[j+1].start){q.polys[idx]=  polys[++i];}
		else										{q.polys[idx].start=p.polys[++j].start;q.polys[idx].p=p.polys[j].p*(-1.0);}
//		if(idx && polys[idx].start==polys[idx-1].start)	{polys[idx-1].p+=polys[idx].p;}
//		else{idx++;}
		idx++;
	}
	return q;
}
template<int Degree>
PPolynomial<Degree>& PPolynomial<Degree>::addScaled(const PPolynomial<Degree>& p,const double& scale){
	int i,j;
	StartingPolynomial<Degree>* oldPolys=polys;
	size_t idx=0,cnt=0,oldPolyCount=polyCount;
	polyCount=0;
	polys=NULL;
	set(oldPolyCount+p.polyCount);
	i=j=-1;
	while(cnt<polyCount){
		if		(j>=int( p.polyCount)-1)				{polys[idx]=oldPolys[++i];}
		else if	(i>=int(oldPolyCount)-1)				{polys[idx].start= p.polys[++j].start;polys[idx].p=p.polys[j].p*scale;}
		else if	(oldPolys[i+1].start<p.polys[j+1].start){polys[idx]=oldPolys[++i];}
		else											{polys[idx].start= p.polys[++j].start;polys[idx].p=p.polys[j].p*scale;}
		if(idx && polys[idx].start==polys[idx-1].start)	{polys[idx-1].p+=polys[idx].p;}
		else{idx++;}
		cnt++;
	}
	free(oldPolys);
	reset(idx);
	return *this;
}
template<int Degree>
template<int Degree2>
PPolynomial<Degree+Degree2> PPolynomial<Degree>::operator * (const PPolynomial<Degree2>& p) const{
	PPolynomial<Degree+Degree2> q;
	StartingPolynomial<Degree+Degree2> *sp;
	int i,j,spCount=int(polyCount*p.polyCount);

	sp=(StartingPolynomial<Degree+Degree2>*)malloc(sizeof(StartingPolynomial<Degree+Degree2>)*spCount);
	for(i=0;i<int(polyCount);i++){
		for(j=0;j<int(p.polyCount);j++){
			sp[i*p.polyCount+j]=polys[i]*p.polys[j];
		}
	}
	q.set(sp,spCount);
	free(sp);
	return q;
}
template<int Degree>
template<int Degree2>
PPolynomial<Degree+Degree2> PPolynomial<Degree>::operator * (const Polynomial::Polynomial1D< Degree2 >& p) const{
	PPolynomial<Degree+Degree2> q;
	q.set(polyCount);
	for(int i=0;i<int(polyCount);i++){
		q.polys[i].start=polys[i].start;
		q.polys[i].p=polys[i].p*p;
	}
	return q;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::scale(const double& s) const{
	PPolynomial q;
	q.set(polyCount);
	for(size_t i=0;i<polyCount;i++){q.polys[i]=polys[i].scale(s);}
	return q;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::shift(const double& s) const{
	PPolynomial q;
	q.set(polyCount);
	for(size_t i=0;i<polyCount;i++){q.polys[i]=polys[i].shift(s);}
	return q;
}
template<int Degree>
PPolynomial< (Degree>0 ? Degree-1 : 0) > PPolynomial<Degree>::derivative(void) const{
	PPolynomial< (Degree>0 ? Degree-1 : 0) > q;
	q.set(polyCount);
	for(size_t i=0;i<polyCount;i++){
		q.polys[i].start=polys[i].start;
		q.polys[i].p=polys[i].p.derivative();
	}
	return q;
}
template<int Degree>
PPolynomial<Degree+1> PPolynomial<Degree>::integral(void) const{
	int i;
	PPolynomial<Degree+1> q;
	q.set(polyCount);
	for(i=0;i<int(polyCount);i++){
		q.polys[i].start=polys[i].start;
		q.polys[i].p=polys[i].p.integral();
		q.polys[i].p-=q.polys[i].p(q.polys[i].start);
	}
	return q;
}
template<int Degree>
PPolynomial<Degree>& PPolynomial<Degree>::operator  += (const double &s){ polys[0].p+=s ; return *this; }
template<int Degree>
PPolynomial<Degree>& PPolynomial<Degree>::operator  -= (const double &s){ polys[0].p-=s ; return *this; }
template<int Degree>
PPolynomial<Degree>& PPolynomial<Degree>::operator  *= (const double &s){
	for(int i=0;i<int(polyCount);i++){polys[i].p*=s;}
	return *this;
}
template<int Degree>
PPolynomial<Degree>& PPolynomial<Degree>::operator  /= (const double &s){
	for(size_t i=0;i<polyCount;i++){polys[i].p/=s;}
	return *this;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator + (const double& s) const{
	PPolynomial q=*this;
	q+=s;
	return q;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator - (const double& s) const{
	PPolynomial q=*this;
	q-=s;
	return q;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator * (const double& s) const{
	PPolynomial q=*this;
	q*=s;
	return q;
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator / (const double& s) const{
	PPolynomial q=*this;
	q/=s;
	return q;
}

template<int Degree>
void PPolynomial<Degree>::printnl(void) const{
	Polynomial::Polynomial1D<Degree> p;

	if(!polyCount){
		Polynomial::Polynomial1D<Degree> p;
		printf("[-Infinity,Infinity]\n");
	}
	else{
		for(size_t i=0;i<polyCount;i++){
			printf("[");
			if		(polys[i  ].start== DBL_MAX){printf("Infinity,");}
			else if	(polys[i  ].start==-DBL_MAX){printf("-Infinity,");}
			else								{printf("%f,",polys[i].start);}
			if(i+1==polyCount)					{printf("Infinity]\t");}
			else if (polys[i+1].start== DBL_MAX){printf("Infinity]\t");}
			else if	(polys[i+1].start==-DBL_MAX){printf("-Infinity]\t");}
			else								{printf("%f]\t",polys[i+1].start);}
			p=p+polys[i].p;
			p.printnl();
		}
	}
	printf("\n");
}
template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::ConstantFunction(const double& radius){
	if(Degree<0){
		fprintf(stderr,"Could not set degree %d polynomial as constant\n",Degree);
		exit(0);
	}
	PPolynomial q;
	q.set(2);

	q.polys[0].start=-radius;
	q.polys[1].start= radius;

	q.polys[0].p.coefficients[0]= 1.0;
	q.polys[1].p.coefficients[0]=-1.0;
	return q;
}

template<>
PPolynomial<0> PPolynomial<0>::GaussianApproximation(const double& width)
{
	return ConstantFunction(width);
}

template<int Degree>
PPolynomial<Degree> PPolynomial<Degree>::GaussianApproximation(const double& width){return PPolynomial<Degree-1>::GaussianApproximation().MovingAverage(width);}
template<int Degree>
PPolynomial<Degree+1> PPolynomial<Degree>::MovingAverage(const double& radius){
	PPolynomial<Degree+1> A;
	Polynomial::Polynomial1D< Degree+1 > p;
	StartingPolynomial<Degree+1>* sps;

	sps=(StartingPolynomial<Degree+1>*)malloc(sizeof(StartingPolynomial<Degree+1>)*polyCount*2);

	for(int i=0;i<int(polyCount);i++){
		sps[2*i  ].start=polys[i].start-radius;
		sps[2*i+1].start=polys[i].start+radius;
		p=polys[i].p.integral()-polys[i].p.integral()(polys[i].start);
		sps[2*i  ].p=p.shift(-radius);
		sps[2*i+1].p=p.shift( radius)*-1;
	}
	A.set(sps,int(polyCount*2));
	free(sps);
	return A*1.0/(2*radius);
}

template<int Degree>
void PPolynomial<Degree>::getSolutions(const double& c,std::vector<double>& roots,const double& EPS,const double& min,const double& max) const{
	Polynomial::Polynomial1D< Degree > p;
	std::vector<double> tempRoots;

	p.setZero();
	for(size_t i=0;i<polyCount;i++){
		p+=polys[i].p;
		if(polys[i].start>max){break;}
		if(i<polyCount-1 && polys[i+1].start<min){continue;}
		p.getSolutions(c,tempRoots,EPS);
		for(size_t j=0;j<tempRoots.size();j++){
			if(tempRoots[j]>polys[i].start && (i+1==polyCount || tempRoots[j]<=polys[i+1].start)){
				if(tempRoots[j]>min && tempRoots[j]<max){roots.push_back(tempRoots[j]);}
			}
		}
	}
}

template<int Degree>
void PPolynomial<Degree>::write(FILE* fp,const int& samples,const double& min,const double& max) const{
	fwrite(&samples,sizeof(int),1,fp);
	for(int i=0;i<samples;i++){
		double x=min+i*(max-min)/(samples-1);
		float v=(*this)(x);
		fwrite(&v,sizeof(float),1,fp);
	}
}
