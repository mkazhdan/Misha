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

inline void DrawSphere( Point3D< double > p , double r , unsigned int res )
{
	unsigned int thetaRes = 2*res , phiRes = res;

	auto DrawSpherePoint = [&]( unsigned int x , unsigned int y )
		{
			double theta = ( 2. * M_PI * x ) / thetaRes;
			double phi = ( M_PI * y ) / phiRes;
			Point3D< double > _v = Point3D< double >( cos( theta ) * sin( phi ) , cos( phi ) , sin( theta ) * sin( phi ) );
			Point3D< double > v = p + _v * r;
			glNormal3d( _v[0] , _v[1] , _v[2] );
			glVertex3d( v[0] , v[1] , v[2] );
		};

	glBegin( GL_TRIANGLES );
	{
		unsigned int j=0;
		for( unsigned int i=0 ; i<2*thetaRes ; i++ )
		{
			DrawSpherePoint(i+1,j+1);
			DrawSpherePoint(i,j+1);
			DrawSpherePoint(i,j);
		}
	}
	for( unsigned int j=1 ; j<phiRes-1 ; j++ )
	{
		for( unsigned int i=0 ; i<2*thetaRes ; i++ )
		{
			DrawSpherePoint(i,j);
			DrawSpherePoint(i+1,j);
			DrawSpherePoint(i+1,j+1);

			DrawSpherePoint(i+1,j+1);
			DrawSpherePoint(i,j+1);
			DrawSpherePoint(i,j);
		}
	}
	{
		unsigned int j=phiRes-1;
		for( unsigned int i=0 ; i<2*thetaRes ; i++ )
		{
			DrawSpherePoint(i,j);
			DrawSpherePoint(i+1,j);
			DrawSpherePoint(i+1,j+1);
		}
	}
	glEnd();
}
inline void DrawDisk( Point3D< double > p , Point3D< double > dir , double r , unsigned int res )
{
	unsigned int thetaRes = 2*res;

	OrthogonalFrame< 3 > frame( &dir , 1 );

	auto DrawDiskPoint = [&]( unsigned int x )
		{
			double theta = ( 2. * M_PI * x ) / thetaRes;
			Point3D< double > _v = frame[1] * cos(theta) + frame[2] * sin(theta);
			Point3D< double > v = p + _v * r;
			glNormal3d( _v[0] , _v[1] , _v[2] );
			glVertex3d( v[0] , v[1] , v[2] );
		};

	glBegin( GL_POLYGON );
	for( unsigned int i=0 ; i<2*thetaRes ; i++ ) DrawDiskPoint(i);
	glEnd();
}

inline void DrawCylinder( Point3D< double > p1 , Point3D< double > p2 , double r , unsigned int res )
{
	unsigned int thetaRes = 2*res;

	Point< double , 3 > v = p2 - p1;
	OrthogonalFrame< 3 > frame( &v , 1 );

	auto DrawCylinderPoint = [&]( unsigned int x , unsigned int y )
	{
		double theta = ( 2. * M_PI * x ) / thetaRes;
		Point3D< double > _v = frame[1] * cos(theta) + frame[2] * sin(theta);
		Point3D< double > v = p2*y + p1*(1-y) + _v * r;
		glNormal3d( _v[0] , _v[1] , _v[2] );
		glVertex3d( v[0] , v[1] , v[2] );
	};

	glBegin( GL_QUADS );
	for( unsigned int i=0 ; i<2*thetaRes ; i++ )
	{
		DrawCylinderPoint(i  ,0);
		DrawCylinderPoint(i+1,0);
		DrawCylinderPoint(i+1,1);
		DrawCylinderPoint(i  ,1);
	}
	glEnd();
}

inline void DrawCone( Point3D< double > p1 , Point3D< double > p2 , double r , unsigned int res )
{
	unsigned int thetaRes = 2*res;

	Point< double , 3 > v = p2 - p1;
	OrthogonalFrame< 3 > frame( &v , 1 );

	auto DrawConePoint = [&]( unsigned int x , unsigned int y )
	{
		double theta = ( 2. * M_PI * x ) / thetaRes;
		Point3D< double > _v = frame[1] * cos(theta) + frame[2] * sin(theta);
		Point3D< double > v = p2*y + p1*(1-y) + _v * r * (1-y);
		glNormal3d( _v[0] , _v[1] , _v[2] );
		glVertex3d( v[0] , v[1] , v[2] );
	};

	glBegin( GL_TRIANGLES );
	for( unsigned int i=0 ; i<2*thetaRes ; i++ )
	{
		DrawConePoint(i  ,0);
		DrawConePoint(i+1,0);
		DrawConePoint(i  ,1);
	}
	glEnd();

	glBegin( GL_POLYGON );
	for( unsigned int i=0 ; i<2*thetaRes ; i++ ) DrawConePoint(i,0);
	glEnd();
}

inline void DrawArrow( Point3D< double > p1 , Point3D< double > p2  , double cylinderRadius , double coneRadius , double coneLength , unsigned int res )
{
#if 1
	DrawDisk( p1 , p1-p2 , cylinderRadius , res );
	DrawCylinder( p1 , p2 , cylinderRadius , res );
	Point< double , 3 > v = ( p2 - p1 );
	v /= Point< double , 3 >::Length( v );
	DrawCone( p2 , p2 + v*coneLength , coneRadius , res );	
#else
	DrawCylinder( p1 , p1+(p2-p1)*(1-tipFraction ) , cylinderRadius , res );	
	DrawCone( p1+(p2-p1)*(1-tipFraction ) , p2 , coneRadius , res );	
#endif
}

inline void DrawCurve( const std::vector< Point3D< double > > & curve , double r , unsigned int res )
{
	auto GetDir = [&]( unsigned int idx )
		{
			Point< double , 3 > dir;
			if( idx>0 )              dir += curve[idx] - curve[idx-1];
			if( idx+1<curve.size() ) dir += curve[idx+1] - curve[idx];
			return dir / Point< double , 3 >::Length( dir );
		};

	unsigned int thetaRes = 2*res;

	if( curve.size()<2 ){ MK_WARN_ONCE( "Cannot draw a curve with less than two vertices: " , curve.size() ); }
	else
	{
		Point< double , 3 > startDir = GetDir( 0 );
		Point< double , 3 > startFrame[2];
		while( Point< double , 3 >::SquareNorm( startFrame[0] )<1e-8 )
		{
			startFrame[0] = RandomSpherePoint< double , 3 >();
			startFrame[0] -= startDir * Point< double , 3 >::Dot( startDir , startFrame[0] );
		}
		startFrame[0] /= Point< double , 3 >::Length( startFrame[0] );
		startFrame[1] = Point< double , 3 >::CrossProduct( startDir , startFrame[0] );

		for( unsigned int i=0 ; i<curve.size()-1 ; i++ )
		{
			Point< double , 3 > endDir = GetDir( i+1 );
			SquareMatrix< double , 3 > R = SimplexProcessing::RodriguesRotation( startDir , endDir );
			Point< double , 3 > endFrame[] = { R * startFrame[0] , R * startFrame[1] };

			auto DrawPoint = [&]( bool front , unsigned int x )
				{
					double theta = ( 2. * M_PI * x ) / thetaRes;
					Point3D< double> v , n;
					if( front )
					{
						n = startFrame[0] * cos(theta) + startFrame[1] * sin(theta);
						v = curve[i] + n * r;
					}
					else
					{
						n = endFrame[0] * cos(theta) + endFrame[1] * sin(theta);
						v = curve[i+1] + n * r;
					}
					glNormal3d( n[0] , n[1] , n[2] );
					glVertex3d( v[0] , v[1] , v[2] );
				};

			glBegin( GL_QUADS );
			for( unsigned int y=0 ; y<2*thetaRes ; y++ )
			{
				DrawPoint( true  , y   );
				DrawPoint( true  , y+1 );
				DrawPoint( false , y+1 );
				DrawPoint( false , y   );
			}
			glEnd();

			startDir = endDir;
			startFrame[0] = endFrame[0];
			startFrame[1] = endFrame[1];
		}
	}
}