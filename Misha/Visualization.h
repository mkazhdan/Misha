/*
Copyright (c) 2018, Michael Kazhdan
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

#ifndef VISUALIZATION_INCLUDED
#define VISUALIZATION_INCLUDED

#define NEW_VISUALIZATION_CODE
#define NEW_CALL_BACK

#include <iomanip>
#include <algorithm>
#include <sstream>
#include <functional>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <vector>
#include <sys/timeb.h>
#include <functional>
#include <algorithm>
#include "Misha/Image.h"
#include "Misha/Array.h"
#include "Misha/Exceptions.h"


namespace MishaK
{
	static const int KEY_UPARROW    = 101;
	static const int KEY_DOWNARROW	= 103;
	static const int KEY_LEFTARROW	= 100;
	static const int KEY_RIGHTARROW	= 102;
	static const int KEY_PGUP		= 104;
	static const int KEY_PGDN		= 105;
	static const int KEY_CTRL_C     =   3;
	static const int KEY_BACK_SPACE =   8;
	static const int KEY_ENTER      =  13;
	static const int KEY_ESC        =  27;

#ifdef VERBOSE_MESSAGING
	inline void AssertOpenGLState( const char *fileName , int line , const char *functionName )
	{
		GLenum errorCode;
		if( ( errorCode=glGetError() )!=GL_NO_ERROR )
		{
			std::string buffer = MakeMessageString( "[OPEN_GL ERROR]" , fileName , line , functionName , "%s (%d)" , gluErrorString( errorCode ) , errorCode );
			std::cerr << buffer << std::endl;
			exit( 0 );
		}
	}
#ifndef ASSERT_OPEN_GL_STATE
#define ASSERT_OPEN_GL_STATE( ... ) MishaK::AssertOpenGLState( __FILE__ , __LINE__ , __FUNCTION__ )
#endif // ASSERT_OPEN_GL_STATE
#else // !VERBOSE_MESSAGING
	inline void AssertOpenGLState( const char *functionName )
	{
		GLenum errorCode;
		if( ( errorCode=glGetError() )!=GL_NO_ERROR )
		{
			std::string buffer = MakeMessageString( "[OPEN_GL ERROR]" , functionName , "%s (%d)" , gluErrorString( errorCode ) , errorCode );
			std::cerr << buffer << std::endl;
			exit( 0 );
		}
	}
#ifndef ASSERT_OPEN_GL_STATE
#define ASSERT_OPEN_GL_STATE( ... ) Misha::AssertOpenGLState( __FUNCTION__ )
#endif // ASSERT_OPEN_GL_STATE
#endif // VERBOSE_MESSAGING

	double Time( void )
	{
#ifdef WIN32
		struct _timeb t;
		_ftime(&t);
		return double(t.time)+double(t.millitm)/1000.0;
#else // WIN32
		struct timeval t;
		gettimeofday(&t,NULL);
		return t.tv_sec+(double)t.tv_usec/1000000;
#endif // WIN32
	}


	struct Font
	{
		static const Font Fonts[];

#ifdef NEW_VISUALIZATION_CODE
		enum Type
		{
			FONT_8_BY_13 ,
			FONT_9_BY_15 ,
			FONT_HELVETICA_10 ,
			FONT_HELVETICA_12 ,
			FONT_HELVETICA_18 ,
			FONT_TIMES_ROMAN_10 ,
			FONT_TIMES_ROMAN_24 ,
			FONT_COUNT
		};

		Font( void ) : _font(nullptr) , _fontHeight(0){}

		const std::string & name( void ) const { return _fontName; }
		unsigned int height( void ) const { return _fontHeight ;}
		void * operator()( void ) const { return _font; }
	protected:
		std::string _fontName;
		void * _font;
		unsigned int _fontHeight;


		Font( std::string fn , void * f , unsigned int fh ) : _fontName(fn) , _font(f) , _fontHeight(fh) {}
#else // !NEW_VISUALIZATION_CODE
		std::string fontName;
		void * font;
		unsigned int fontHeight;

		Font( std::string fn , void * f , unsigned int fh ) : fontName(fn) , font(f) , fontHeight(fh) {}
		static unsigned int FontNum();
#endif // NEW_VISUALIZATION_CODE
	};

	const Font Font::Fonts[] =
	{
		Font( "Bitmap 8x13" , GLUT_BITMAP_8_BY_13 , 13 ) ,
		Font( "Bitmap 9x15" , GLUT_BITMAP_9_BY_15 , 15 ) ,
		Font( "Helvetica 10" , GLUT_BITMAP_HELVETICA_10 , 10 ) ,
		Font( "Helvetica 12" , GLUT_BITMAP_HELVETICA_12 , 12 ) ,
		Font( "Helvetica 18" , GLUT_BITMAP_HELVETICA_18 , 18 ) ,
		Font( "Times-Roman 10" , GLUT_BITMAP_TIMES_ROMAN_10 , 10 ) ,
		Font( "Times-Roman 24" , GLUT_BITMAP_TIMES_ROMAN_24 , 24 ) ,
	};
#ifdef NEW_VISUALIZATION_CODE
#else // !NEW_VISUALIZATION_CODE
	unsigned int Font::FontNum( void ){ return sizeof( Fonts ) / sizeof( Font ); }
#endif // NEW_VISUALIZATION_CODE


#ifdef NEW_VISUALIZATION_CODE
	using CallBackFunction = std::function< void ( std::string ) >;

	struct KeyboardCallBack
	{
		template< typename T >
		static CallBackFunction GetCallBackFunction( T * obj , void (T::*memberFunctionPointer)( std::string ) )
		{
			return [obj,memberFunctionPointer]( std::string str ){ return (obj->*memberFunctionPointer)( str ); };
		}

		struct Modifiers
		{
			bool alt , ctrl;
			Modifiers( bool a=false , bool c=false ) : alt(a) , ctrl(c) {};
			bool operator == ( const Modifiers &m ) const { return alt==m.alt && ctrl==m.ctrl; }
			bool operator != ( const Modifiers &m ) const { return alt!=m.alt || ctrl!=m.ctrl; }
		};
		char key;
		std::string prompt;
		std::string description;
		CallBackFunction callBackFunction;
		Modifiers modifiers;
		KeyboardCallBack( char key , Modifiers modifiers , std::string description ,                      CallBackFunction callBackFunction );
		KeyboardCallBack( char key , Modifiers modifiers , std::string description , std::string prompt , CallBackFunction callBackFunction );
		template< typename T >
		KeyboardCallBack( char key , Modifiers modifiers , std::string description ,                      T * obj , void (T::*memberFunctionPointer)( std::string ) );
		template< typename T >
		KeyboardCallBack( char key , Modifiers modifiers , std::string description , std::string prompt , T * obj , void (T::*memberFunctionPointer)( std::string ) );
	};
#endif // NEW_VISUALIZATION_CODE

	template< typename DerivedViewableType >
	struct Viewable
	{
#ifdef NEW_VISUALIZATION_CODE
#else // !NEW_VISUALIZATION_CODE
		using CallBackFunction = std::function< void ( DerivedViewableType * , const char * ) >;
#endif // NEW_VISUALIZATION_CODE

#ifdef NEW_VISUALIZATION_CODE
#else // !NEW_VISUALIZATION_CODE
		struct KeyboardCallBack
		{
			struct Modifiers
			{
				bool alt , ctrl;
				Modifiers( bool a=false , bool c=false ) : alt(a) , ctrl(c) {};
				bool operator == ( const Modifiers &m ) const { return alt==m.alt && ctrl==m.ctrl; }
				bool operator != ( const Modifiers &m ) const { return alt!=m.alt || ctrl!=m.ctrl; }
			};
			char key;
			char prompt[1024];
			char description[1024];
			CallBackFunction callBackFunction;
			DerivedViewableType* viewable;
			Modifiers modifiers;
			KeyboardCallBack( DerivedViewableType * viewable , char key , Modifiers modifiers , const char * description ,                       CallBackFunction callBackFunction );
			KeyboardCallBack( DerivedViewableType * viewable , char key , Modifiers modifiers , const char * description , const char * prompt , CallBackFunction callBackFunction );
		};
#endif // NEW_VISUALIZATION_CODE

		struct Viewer
		{
			static DerivedViewableType *viewable;
#ifdef NEW_VISUALIZATION_CODE
			static void Run( DerivedViewableType * viewable , int argc , char * argv[] , std::string windowName="" );
#else // !NEW_VISUALIZATION_CODE
			static void Run( DerivedViewableType * viewable , int argc , char * argv[] , const char * windowName="" );
#endif // NEW_VISUALIZATION_CODE
			static void Idle             ( void );
			static void KeyboardFunc     ( unsigned char key , int x , int y );
			static void KeyboardUpFunc   ( unsigned char key , int x , int y );
			static void SpecialFunc      ( int key, int x, int y );
			static void SpecialUpFunc    ( int key, int x, int y );
			static void Display          ( void );
			static void Reshape          ( int w , int h );
			static void MouseFunc        ( int button , int state , int x , int y );
			static void MouseWheelFunc   ( int button , int state , int x , int y );
			static void MotionFunc       ( int x , int y );
			static void PassiveMotionFunc( int x , int y );
		};

	protected:
		const double _MIN_FPS_TIME = 0.5;
		double _lastFPSTime;
		int _lastFPSCount;
		double _fps;
		int _currentFrame , _totalFrames;
		bool _exitAfterSnapshot , _exitAfterVideo;
	public:
#ifdef NEW_VISUALIZATION_CODE
		unsigned int screenWidth , screenHeight;
		Font font , promptFont;
#else // !NEW_VISUALIZATION_CODE
		int screenWidth , screenHeight;
		void *font , *promptFont;
		int fontHeight , promptFontHeight;
#endif // NEW_VISUALIZATION_CODE
		bool showHelp , showInfo , showFPS;
		CallBackFunction promptCallBack;
#ifdef NEW_VISUALIZATION_CODE
		std::string promptString;
		unsigned int originalPromptLength;
		std::string snapshotName;
		std::string videoHeader;
#else // !NEW_VISUALIZATION_CODE
		char promptString[1024];
		int promptLength;
		char* snapshotName;
		char* videoHeader;
#endif // NEW_VISUALIZATION_CODE
		bool flushImage;
		std::function< void ( void ) > quitFunction;

		std::vector< KeyboardCallBack > callBacks;
#ifdef NEW_VISUALIZATION_CODE
		std::vector< std::string > info;
#else // !NEW_VISUALIZATION_CODE
		std::vector< char * > info;
#endif // NEW_VISUALIZATION_CODE
		Viewable( void );
		void setFont( unsigned int idx );
#ifdef NEW_VISUALIZATION_CODE
		void setSnapshot( std::string sName , bool exitAfterSnapshot=true );
		void setVideo( std::string sName , int frames , bool exitAfterVideo=true );
#else // !NEW_VISUALIZATION_CODE
		void setSnapshot( const char* sName , bool exitAfterSnapshot=true );
		void setVideo( const char* sName , int frames , bool exitAfterVideo=true );
#endif // NEW_VISUALIZATION_CODE
		virtual void display( void ) {}
		virtual void idle( void ) {}
		virtual void keyboardFunc( unsigned char key , int x , int y ) {}
		virtual void keyboardUpFunc( unsigned char key , int x , int y ) {}
		virtual void specialFunc( int key, int x, int y ) {}
		virtual void specialUpFunc( int key, int x, int y ) {}
		virtual void mouseFunc( int button , int state , int x , int y ) {}
		virtual void mouseWheelFunc( int button , int state , int x , int y ) {}
		virtual void motionFunc( int x , int y ) {}
		virtual void passiveMotionFunc( int x , int y ) {}
		virtual void reshape( int w , int h )
		{
			screenWidth = w , screenHeight = h;
			glViewport( 0 , 0 , screenWidth , screenHeight );
		}

		void Idle             ( void );
		void KeyboardFunc     ( unsigned char key , int x , int y );
		void KeyboardUpFunc   ( unsigned char key , int x , int y );
		void SpecialFunc      ( int key, int x, int y );
		void SpecialUpFunc    ( int key, int x, int y );
		void Display          ( void );
		void Reshape          ( int w , int h );
		void MouseFunc        ( int button , int state , int x , int y );
		void MouseWheelFunc   ( int button , int state , int x , int y );
		void MotionFunc       ( int x , int y );
		void PassiveMotionFunc( int x , int y );

#ifdef NEW_VISUALIZATION_CODE
		void           exitCallBack( std::string ){ exit( 0 ); }
		void           quitCallBack( std::string ){ quitFunction() , exit( 0 ); }
		void      toggleFPSCallBack( std::string ){ showFPS  = !showFPS ; }
		void     toggleHelpCallBack( std::string ){ showHelp = !showHelp; }
		void     toggleInfoCallBack( std::string ){ showInfo = !showInfo; }
		void setFrameBufferCallBack( std::string );

#ifdef NEW_VISUALIZATION_CODE
		static void WriteLeftString( int x , int y , const Font & font , std::string );
		static unsigned int StringWidth( const Font & font , std::string );
#else // !NEW_VISUALIZATION_CODE
		static void WriteLeftString( int x , int y , void * font , std::string );
		static unsigned int StringWidth( void * font , std::string );
#endif // NEW_VISUALIZATION_CODE
		void writeLeftString( int x , int y , std::string ) const;
		void writeRightString( int x , int y , std::string ) const;
		void writeCenterString( int x , int y , std::string ) const;
		unsigned int stringWidth( std::string ) const;

		void saveFrameBuffer( std::string fileName , int whichBuffer=GL_BACK );
		void setPromptCallBack( std::string prompt , CallBackFunction callBackFunction );

		void addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , std::string description ,                      CallBackFunction callBackFunction );
		void addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , std::string description , std::string prompt , CallBackFunction callBackFunction );
		void addCallBack( char key , std::string description ,                                                                       CallBackFunction callBackFunction );
		void addCallBack( char key , std::string description , std::string prompt ,                                                  CallBackFunction callBackFunction );
		void addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , std::string description ,                      void (DerivedViewableType::*memberFunctionPointer)( std::string ) );
		void addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , std::string description , std::string prompt , void (DerivedViewableType::*memberFunctionPointer)( std::string ) );
		void addCallBack( char key , std::string description ,                                                                       void (DerivedViewableType::*memberFunctionPointer)( std::string ) );
		void addCallBack( char key , std::string description , std::string prompt ,                                                  void (DerivedViewableType::*memberFunctionPointer)( std::string ) );
#else // !NEW_VISUALIZATION_CODE
		static void       ExitCallBack( DerivedViewableType *   , const char * ){ exit( 0 ); }
		static void       QuitCallBack( DerivedViewableType * v , const char * ){ v->quitFunction() , exit( 0 ); }
		static void  ToggleFPSCallBack( DerivedViewableType * v , const char * ){ v->showFPS  = !v->showFPS ; }
		static void ToggleHelpCallBack( DerivedViewableType * v , const char * ){ v->showHelp = !v->showHelp; }
		static void ToggleInfoCallBack( DerivedViewableType * v , const char * ){ v->showInfo = !v->showInfo; }
		static void SetFrameBufferCallBack( DerivedViewableType* v , const char* prompt );

		static void WriteLeftString( int x , int y , void * font , const char * format , ... );
		static int StringWidth( void * font , const char * format , ... );
		void writeLeftString( int x , int y , const char * format , ... ) const;
		void writeRightString( int x , int y , const char * format , ... ) const;
		void writeCenterString( int x , int y , const char * format , ... ) const;
		int stringWidth( const char * format , ... ) const;

		void saveFrameBuffer( const char* fileName , int whichBuffer=GL_BACK );
		void setPromptCallBack( const char *prompt , CallBackFunction callBackFunction );

		void addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , const char * description ,                      CallBackFunction callBackFunction );
		void addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , const char * description , const char *prompt , CallBackFunction callBackFunction );
		void addCallBack( char key , const char *description ,                      CallBackFunction callBackFunction );
		void addCallBack( char key , const char *description , const char *prompt , CallBackFunction callBackFunction );
#endif // NEW_VISUALIZATION_CODE
	};

	//////////////////////
	// KeyboardCallBack //
	//////////////////////
#ifdef NEW_VISUALIZATION_CODE
	KeyboardCallBack::KeyboardCallBack( char key , Modifiers modifiers , std::string description , CallBackFunction  callBackFunction )
		: modifiers(modifiers) , key(key) , description(description) , callBackFunction(callBackFunction)
	{}

	KeyboardCallBack::KeyboardCallBack( char key , Modifiers modifiers , std::string description , std::string prompt , CallBackFunction  callBackFunction )
		: modifiers(modifiers) , key(key) , description(description) , prompt(prompt) , callBackFunction(callBackFunction)
	{}

	template< typename T >
	KeyboardCallBack::KeyboardCallBack( char key , Modifiers modifiers , std::string description , T * obj , void (T::*memberFunctionPointer)( std::string ) )
		: modifiers(modifiers) , key(key) , description(description) , callBackFunction( GetCallBackFunction(obj,memberFunctionPointer) )
	{}

	template< typename T >
	KeyboardCallBack::KeyboardCallBack( char key , Modifiers modifiers , std::string description , std::string prompt , T * obj , void (T::*memberFunctionPointer)( std::string ) )
		: modifiers(modifiers) , key(key) , description(description) , prompt(prompt) , callBackFunction( GetCallBackFunction(obj,memberFunctionPointer) )
	{}
#else // !NEW_VISUALIZATION_CODE
	template< typename DerivedViewableType >
	Viewable< DerivedViewableType >::KeyboardCallBack::KeyboardCallBack( DerivedViewableType* viewable , char key , Modifiers modifiers , const char * description , CallBackFunction  callBackFunction )
	{
		this->modifiers = modifiers;
		this->viewable = viewable;
		this->key = key;
		strcpy( this->description , description );
		prompt[0] = 0;
		this->callBackFunction = callBackFunction;
	}

	template< typename DerivedViewableType >
	Viewable< DerivedViewableType >::KeyboardCallBack::KeyboardCallBack( DerivedViewableType* viewable , char key , Modifiers modifiers , const char * description , const char * prompt , CallBackFunction  callBackFunction )
	{
		this->modifiers = modifiers;
		this->viewable = viewable;
		this->key = key;
		strcpy( this->description , description );
		strcpy( this->prompt , prompt );
		this->callBackFunction = callBackFunction;
	}
#endif // NEW_VISUALIZATION_CODE

	//////////////////////
	// Viewable::Viewer //
	//////////////////////

	template< typename DerivedViewableType > DerivedViewableType* Viewable< DerivedViewableType >::Viewer::viewable = NULL;

	template< typename DerivedViewableType >
#ifdef NEW_VISUALIZATION_CODE
	void Viewable< DerivedViewableType >::Viewer::Run( DerivedViewableType * v , int argc , char * argv[] , std::string windowName )
#else // !NEW_VISUALIZATION_CODE
	void Viewable< DerivedViewableType >::Viewer::Run( DerivedViewableType * v , int argc , char * argv[] , const char * windowName )
#endif // NEW_VISUALIZATION_CODE
	{
		viewable = v;
		glutInitDisplayMode( GLUT_DEPTH | GLUT_RGBA | GLUT_DOUBLE );
		glutInitWindowSize( viewable->screenWidth , viewable->screenHeight );
		glutInit( &argc , argv );
#ifdef NEW_VISUALIZATION_CODE
		glutCreateWindow( windowName.c_str() );
#else // !NEW_VISUALIZATION_CODE
		glutCreateWindow( windowName );
#endif // NEW_VISUALIZATION_CODE

		if( glewInit()!=GLEW_OK ) MK_ERROR_OUT( "glewInit failed" );

		glutIdleFunc         ( Idle );
		glutDisplayFunc      ( Display );
		glutReshapeFunc      ( Reshape );
		glutMouseFunc        ( MouseFunc );
		glutMotionFunc       ( MotionFunc );
		glutPassiveMotionFunc( PassiveMotionFunc );
		glutKeyboardFunc     ( KeyboardFunc );
		glutKeyboardUpFunc   ( KeyboardUpFunc );
		glutSpecialFunc      ( SpecialFunc );
		glutSpecialUpFunc    ( SpecialUpFunc );
		glutMouseWheelFunc   ( MouseWheelFunc );

		glutMainLoop();
	}

	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::Idle( void ){ viewable->Idle(); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::KeyboardFunc( unsigned char key , int x , int y ){ viewable->KeyboardFunc( key , x , y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::KeyboardUpFunc( unsigned char key , int x , int y ){ viewable->KeyboardUpFunc( key , x , y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::SpecialFunc( int key , int x , int y ){ viewable->SpecialFunc( key , x ,  y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::SpecialUpFunc( int key , int x , int y ){ viewable->SpecialUpFunc( key , x ,  y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::Display( void ){ viewable->Display(); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::Reshape( int w , int h ){ viewable->Reshape( w , h ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::MouseFunc( int button , int state , int x , int y ){ viewable->MouseFunc( button , state , x , y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::MouseWheelFunc( int button , int state , int x , int y ){ viewable->MouseWheelFunc( button , state , x , y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::MotionFunc( int x , int y ){ viewable->MotionFunc( x , y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::PassiveMotionFunc( int x , int y ){ viewable->PassiveMotionFunc( x , y ); }

	//////////////
	// Viewable //
	//////////////
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Reshape( int w , int h ){ reshape( w , h ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::MouseFunc( int button , int state , int x , int y ){ mouseFunc( button , state , x , y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::MouseWheelFunc( int button , int state , int x , int y ){ mouseWheelFunc( button , state , x , y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::MotionFunc( int x , int y ){ motionFunc( x , y );}
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::PassiveMotionFunc( int x , int y ){ passiveMotionFunc( x , y );}
	template< typename DerivedViewableType > 
	void Viewable< DerivedViewableType >::Idle( void )
	{
#ifdef NEW_VISUALIZATION_CODE
		if( snapshotName.size() )
#else // !NEW_VISUALIZATION_CODE
		if( snapshotName )
#endif // NEW_VISUALIZATION_CODE
		{
			if( flushImage )
			{
				flushImage = false;
				glutPostRedisplay();
				return;
			}
			else
			{
				saveFrameBuffer( snapshotName , GL_FRONT );
#ifdef NEW_VISUALIZATION_CODE
				snapshotName.clear();
#else // !NEW_VISUALIZATION_CODE
				delete[] snapshotName;
				snapshotName = NULL;
#endif // NEW_VISUALIZATION_CODE
				if( _exitAfterSnapshot ) exit( 0 );
			}
		}
#ifdef NEW_VISUALIZATION_CODE
		else if( videoHeader.size() && _currentFrame<_totalFrames )
		{
			std::stringstream ss;
//			ss << videoHeader << std::setfill('0') << std::setw(4) << _currentFrame << ".png";
			ss << videoHeader << std::setfill('0') << std::setw(4) << _currentFrame << ".jpg";
			saveFrameBuffer( ss.str() , GL_FRONT );
			_currentFrame++;
			if( _currentFrame==_totalFrames && _exitAfterVideo ) exit( 0 );
		}
#else // !NEW_VISUALIZATION_CODE
		else if( videoHeader && _currentFrame<_totalFrames )
		{
			char snapshotName[512];
			sprintf( snapshotName , "%s.%04d.jpg" , videoHeader , _currentFrame );
//			sprintf( snapshotName , "%s.%04d.png" , videoHeader , _currentFrame );
			saveFrameBuffer( snapshotName , GL_FRONT );
			_currentFrame++;
			if( _currentFrame==_totalFrames && _exitAfterVideo ) exit( 0 );
		}
#endif // NEW_VISUALIZATION_CODE
		idle();
	}

#ifdef NEW_VISUALIZATION_CODE
	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::setPromptCallBack( std::string prompt , CallBackFunction callBackFunction )
	{
		promptString = prompt + std::string( ": " );
		originalPromptLength = static_cast< unsigned int >( promptString.size() );
		promptCallBack = callBackFunction;
	}
#else // !NEW_VISUALIZATION_CODE
	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::setPromptCallBack( const char *prompt , CallBackFunction callBackFunction )
	{
		sprintf( promptString , "%s: " , prompt );
		promptLength = (int)strlen( promptString );
		promptCallBack = callBackFunction;
	}
#endif // NEW_VISUALIZATION_CODE

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::KeyboardFunc( unsigned char key , int x , int y )
	{
		if( promptCallBack )
		{
#ifdef NEW_VISUALIZATION_CODE
			size_t len = promptString.size();
			if( key==KEY_BACK_SPACE )
			{
				if( len>originalPromptLength ) promptString.pop_back();
			}
#else // !NEW_VISUALIZATION_CODE
			size_t len = strlen( promptString );
			if( key==KEY_BACK_SPACE )
			{
				if( len>promptLength ) promptString[len-1] = 0;
			}
#endif // NEW_VISUALIZATION_CODE
			else if( key==KEY_ENTER )
			{
#ifdef NEW_VISUALIZATION_CODE
				promptCallBack( promptString.substr( originalPromptLength ) );
				promptString.clear();
				originalPromptLength = 0;
#else // !NEW_VISUALIZATION_CODE
				promptCallBack( (DerivedViewableType*)this , promptString+promptLength );
				promptString[0] = 0;
				promptLength = 0;
#endif // NEW_VISUALIZATION_CODE
				promptCallBack = nullptr;
			}
			else if( key==KEY_CTRL_C )
			{
#ifdef NEW_VISUALIZATION_CODE
				promptString.clear();
				originalPromptLength = 0;
#else // !NEW_VISUALIZATION_CODE
				promptString[0] = 0;
				promptLength = 0;
#endif // NEW_VISUALIZATION_CODE
				promptCallBack = nullptr;
			}
			else if( key>=32 && key<=126 ) // ' ' to '~'
			{
#ifdef NEW_VISUALIZATION_CODE
				promptString.push_back( key );
#else // !NEW_VISUALIZATION_CODE
				promptString[ len ] = key;
				promptString[ len+1 ] = 0;
#endif // NEW_VISUALIZATION_CODE
			}
			glutPostRedisplay();
			return;
		}
		switch( key )
		{
		case KEY_CTRL_C:
			exit( 0 );
			break;
		default:
		{
			int m = glutGetModifiers();
			typename KeyboardCallBack::Modifiers modifiers( m & GLUT_ACTIVE_ALT , m & GLUT_ACTIVE_CTRL );
			for( unsigned int i=0 ; i<callBacks.size() ; i++ ) if( callBacks[i].key==key && callBacks[i].modifiers==modifiers )
			{
#ifdef NEW_VISUALIZATION_CODE
				if( callBacks[i].prompt.size() ) setPromptCallBack( callBacks[i].prompt , callBacks[i].callBackFunction );
				else callBacks[i].callBackFunction( std::string() );
#else // !NEW_VISUALIZATION_CODE
				if( strlen( callBacks[i].prompt ) ) setPromptCallBack( callBacks[i].prompt , callBacks[i].callBackFunction );
				else callBacks[i].callBackFunction( (DerivedViewableType*)this , NULL );
#endif // NEW_VISUALIZATION_CODE
				break;
			}
		}
		}
		keyboardFunc( key , x , y );
		glutPostRedisplay();
	}

	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::KeyboardUpFunc( unsigned char key , int x , int y ){ keyboardUpFunc( key , x , y ); }

	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::SpecialFunc( int key , int x , int y ){ specialFunc( key , x , y ); }

	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::SpecialUpFunc( int key , int x , int y ){ specialUpFunc( key , x , y ); }

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::Display( void )
	{
		glClearColor( 1 , 1 , 1 , 1 );
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

		display();

		_lastFPSCount++;
		double t = Time();
		if( t-_lastFPSTime > _MIN_FPS_TIME )
		{
			_fps = (double)_lastFPSCount / (t-_lastFPSTime);
			_lastFPSCount = 0;
			_lastFPSTime = t;
		}
#ifdef NEW_VISUALIZATION_CODE
		if( showFPS )
		{
			std::stringstream ss;
			Miscellany::StreamFloatPrecision sfp( ss , 2 );
			ss << screenWidth << " x " << screenHeight << " @ " << _fps;
#ifdef NEW_VISUALIZATION_CODE
			writeRightString( 5 , screenHeight - font.height() - 5 , ss.str() );
#else // !NEW_VISUALIZATION_CODE
			writeRightString( 5 , screenHeight - fontHeight - 5 , ss.str() );
#endif // NEW_VISUALIZATION_CODE
		}
#else // !NEW_VISUALIZATION_CODE
		if( showFPS ) writeRightString( 5 , screenHeight - fontHeight - 5 , "%d x %d @ %.2f" , screenWidth , screenHeight , _fps );
#endif // NEW_VISUALIZATION_CODE

		GLboolean writeMask;
		glGetBooleanv( GL_DEPTH_WRITEMASK , &writeMask );
		glDepthMask( GL_FALSE );

		glDisable( GL_LIGHTING );
#ifdef NEW_VISUALIZATION_CODE
		unsigned int offset = font.height()/2;
#else // !NEW_VISUALIZATION_CODE
		int offset = fontHeight/2;
#endif // NEW_VISUALIZATION_CODE
		if( showHelp )
		{
			auto CallBackString = [&]( int i )
			{
				std::stringstream stream;
				stream << "\'" << callBacks[i].key << "\'";
				if( callBacks[i].modifiers.ctrl ) stream << "+[CTRL]";
				if( callBacks[i].modifiers.alt ) stream << "+[ALT]";
				stream << ": " << callBacks[i].description;
				return stream.str();
			};
			{
				GLint vp[4];
				glGetIntegerv( GL_VIEWPORT , vp );

				glMatrixMode( GL_PROJECTION );
				glPushMatrix();
				glLoadIdentity();
				glOrtho( vp[0] , vp[2] , vp[1] , vp[3] , 0 , 1 );

				glMatrixMode( GL_MODELVIEW );
				glPushMatrix();
				glLoadIdentity();

				int x=0 , y = offset;
#ifdef NEW_VISUALIZATION_CODE
				for( unsigned int i=0 ; i<callBacks.size() ; i++ ) if( callBacks[i].description.size() ) x = std::max< int >( x , stringWidth( CallBackString(i).c_str() ) ) , y += font.height() + offset;
#else // !NEW_VISUALIZATION_CODE
				for( unsigned int i=0 ; i<callBacks.size() ; i++ ) if( strlen( callBacks[i].description ) ) x = std::max< int >( x , stringWidth( CallBackString(i).c_str() ) ) , y += fontHeight + offset;
#endif // NEW_VISUALIZATION_CODE

				GLint srcAlpha , dstAlpha;
				glGetIntegerv( GL_BLEND_SRC_ALPHA , &srcAlpha );
				glGetIntegerv( GL_BLEND_DST_ALPHA , &dstAlpha );

				glEnable( GL_BLEND );
				glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
				glBegin( GL_QUADS );
				glColor4f( 1.f , 1.f , 1.f , 0.5f );
				glVertex2i( screenWidth-5 , 0 ) , glVertex2i( screenWidth-(x+15) , 0 ) , glVertex2i( screenWidth-(x+15) , y ) , glVertex2i( screenWidth-5 , y );
				glEnd();
				glBlendFunc( srcAlpha , dstAlpha );
				glDisable( GL_BLEND );
				glDisable( GL_DEPTH_TEST );
				glLineWidth( 2.f );
				glBegin( GL_LINE_LOOP );
				glColor4f( 0.f , 0.f , 0.f , 1.f );
				glVertex2i( screenWidth-5 , 0 ) , glVertex2i( screenWidth-(x+15) , 0 ) , glVertex2i( screenWidth-(x+15) , y ) , glVertex2i( screenWidth-5 , y );
				glEnd();

				glMatrixMode( GL_PROJECTION );
				glPopMatrix();

				glMatrixMode( GL_MODELVIEW );
				glPopMatrix();
			}

			{
				int y = offset , width = 0;

#ifdef NEW_VISUALIZATION_CODE
				for( unsigned int i=0 ; i<callBacks.size() ; i++ ) if( callBacks[i].description.size() )
					width = std::max< int >( width , stringWidth( CallBackString(i).c_str() ) );
				for( unsigned int i=0 ; i<callBacks.size() ; i++ ) if( callBacks[i].description.size() )
					writeLeftString( screenWidth - 10 - width , y , CallBackString(i).c_str() ) , y += font.height() + offset;
#else // !NEW_VISUALIZATION_CODE
				for( unsigned int i=0 ; i<callBacks.size() ; i++ ) if( strlen( callBacks[i].description ) )
					width = std::max< int >( width , stringWidth( CallBackString(i).c_str() ) );
				for( unsigned int i=0 ; i<callBacks.size() ; i++ ) if( strlen( callBacks[i].description ) )
					writeLeftString( screenWidth - 10 - width , y , CallBackString(i).c_str() ) , y += fontHeight + offset;
#endif // NEW_VISUALIZATION_CODE
			}
		}
		if( showInfo && info.size() )
		{
			{
				GLint vp[4];
				glGetIntegerv( GL_VIEWPORT , vp );

				glMatrixMode( GL_MODELVIEW );
				glPushMatrix();
				glLoadIdentity();
				glMatrixMode( GL_PROJECTION );
				glPushMatrix();
				glLoadIdentity();
				glOrtho( vp[0] , vp[2] , vp[1] , vp[3] , 0 , 1 );
				int x=0 , y = offset;
#ifdef NEW_VISUALIZATION_CODE
				for( unsigned int i=0 ; i<info.size() ; i++ ) if( info[i].size() ) x = std::max< int >( x , glutBitmapLength( font() , reinterpret_cast< const unsigned char * >( info[i].c_str() ) ) ) , y += font.height() + offset;
#else // !NEW_VISUALIZATION_CODE
				for( int i=0 ; i<info.size() ; i++ ) if( strlen( info[i] ) ) x = std::max< int >( x , glutBitmapLength( font , (unsigned char *) info[i] ) ) , y += fontHeight + offset;
#endif // NEW_VISUALIZATION_CODE
				glEnable( GL_BLEND );
				GLint srcAlpha , dstAlpha;
				glGetIntegerv( GL_BLEND_SRC_ALPHA , &srcAlpha );
				glGetIntegerv( GL_BLEND_DST_ALPHA , &dstAlpha );
				glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
				glBegin( GL_QUADS );
				glColor4f( 1.f , 1.f , 1.f , 0.5f );
				glVertex2i( 5 , 0 ) , glVertex2i( x+15 , 0 ) , glVertex2i( x+15 , y ) , glVertex2i( 5 , y );
				glEnd();
				glBlendFunc( srcAlpha , dstAlpha );
				glDisable( GL_BLEND );
				glDisable( GL_DEPTH_TEST );
				glLineWidth( 2.f );
				glBegin( GL_LINE_LOOP );
				glColor4f( 0.f , 0.f , 0.f , 1.f );
				glVertex2i( 5 , 0 ) , glVertex2i( x+15 , 0 ) , glVertex2i( x+15 , y ) , glVertex2i( 5 , y );
				glEnd();

				glMatrixMode( GL_PROJECTION );
				glPopMatrix();

				glMatrixMode( GL_MODELVIEW );
				glPopMatrix();
			}
			{
				int y = offset;
#ifdef NEW_VISUALIZATION_CODE
				for( unsigned int i=0 ; i<info.size() ; i++ ) if( info[i].size() ) writeLeftString( 10 , y , info[i] ) , y += font.height() + offset;
#else // !NEW_VISUALIZATION_CODE
				for( int i=0 ; i<info.size() ; i++ ) if( strlen( info[i] ) )
					writeLeftString( 10 , y , "%s" , info[i] ) , y += fontHeight + offset;
#endif // NEW_VISUALIZATION_CODE
			}
		}

#ifdef NEW_VISUALIZATION_CODE
		if( promptString.size() )
#else // !NEW_VISUALIZATION_CODE
		if( strlen( promptString ) )
#endif // NEW_VISUALIZATION_CODE
		{
#ifdef NEW_VISUALIZATION_CODE
			Font _font = font;
			font = promptFont;
#else // !NEW_VISUALIZATION_CODE
			void * _font = font;
			int _fontHeight = fontHeight;
			font = promptFont;
			fontHeight = promptFontHeight;
#endif // NEW_VISUALIZATION_CODE

			int sw = StringWidth ( font , promptString );
			glColor4f( 1.f , 1.f , 1.f , 0.5 );
			glEnable( GL_BLEND );
			GLint srcAlpha , dstAlpha;
			glGetIntegerv( GL_BLEND_SRC_ALPHA , &srcAlpha );
			glGetIntegerv( GL_BLEND_DST_ALPHA , &dstAlpha );
			glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
			glBegin( GL_QUADS );
			{
				glVertex2i(     0 , screenHeight );
				glVertex2i( sw+20 , screenHeight );
#ifdef NEW_VISUALIZATION_CODE
				glVertex2i( sw+20 , screenHeight-font.height()*2 );
				glVertex2i(     0 , screenHeight-font.height()*2 );
#else // !NEW_VISUALIZATION_CODE
				glVertex2i( sw+20 , screenHeight-fontHeight*2 );
				glVertex2i(     0 , screenHeight-fontHeight*2 );
#endif // NEW_VISUALIZATION_CODE
			}
			glEnd();
			glBlendFunc( srcAlpha , dstAlpha );
			glDisable( GL_BLEND );
			glColor4f( 0.f , 0.f , 0.f , 1.f );
			glLineWidth( 2.f );
			glBegin( GL_LINE_LOOP );
			{
				glVertex2i(     0 , screenHeight );
				glVertex2i( sw+20 , screenHeight );
#ifdef NEW_VISUALIZATION_CODE
				glVertex2i( sw+20 , screenHeight-font.height()*2 );
				glVertex2i(     0 , screenHeight-font.height()*2 );
#else // !NEW_VISUALIZATION_CODE
				glVertex2i( sw+20 , screenHeight-fontHeight*2 );
				glVertex2i(     0 , screenHeight-fontHeight*2 );
#endif // NEW_VISUALIZATION_CODE
			}
			glEnd();
#ifdef NEW_VISUALIZATION_CODE
			writeLeftString( 10 , screenHeight-font.height()-font.height()/2 , promptString.c_str() );
			font = _font;
#else // !NEW_VISUALIZATION_CODE
			writeLeftString( 10 , screenHeight-fontHeight-fontHeight/2 , promptString );
			font = _font;
			fontHeight = _fontHeight;
#endif // NEW_VISUALIZATION_CODE
		}
		if( writeMask ) glDepthMask( GL_TRUE );
		glutSwapBuffers();
	}

#ifdef NEW_VISUALIZATION_CODE
	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::WriteLeftString( int x , int y , const Font & font , std::string str )
#else // !NEW_VISUALIZATION_CODE
	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::WriteLeftString( int x , int y , void * font , const char * format , ... )
#endif // NEW_VISUALIZATION_CODE
	{
#ifdef NEW_VISUALIZATION_CODE
#else // !NEW_VISUALIZATION_CODE
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
#endif // NEW_VISUALIZATION_CODE

		GLint vp[4];
		glGetIntegerv( GL_VIEWPORT , vp );

		glMatrixMode( GL_PROJECTION );
		glPushMatrix();
		glLoadIdentity();
		glOrtho( vp[0] , vp[2] , vp[1] , vp[3] , 0 , 1 );

		glMatrixMode( GL_MODELVIEW );
		glPushMatrix();
		glLoadIdentity();

		GLint matrixMode;
		glGetIntegerv( GL_MATRIX_MODE , &matrixMode );
		int depth = glIsEnabled( GL_DEPTH_TEST );
		int lighting = glIsEnabled( GL_LIGHTING );
		glDisable( GL_DEPTH_TEST );
		glDisable( GL_LIGHTING );
		glColor4f( 0 , 0 , 0 , 1 );
		glRasterPos2i( x , y );
#ifdef NEW_VISUALIZATION_CODE
		for( unsigned int i=0 ; i<str.size() ; i++ ) glutBitmapCharacter( font() , str[i] );
#else // !NEW_VISUALIZATION_CODE
		int len = int( strlen( str ) );
		for( int i=0 ; i<len ; i++ ) glutBitmapCharacter( font , str[i] );
#endif // NEW_VISUALIZATION_CODE
		if( depth ) glEnable( GL_DEPTH_TEST );
		if( lighting ) glEnable( GL_LIGHTING );

		glMatrixMode( GL_PROJECTION );
		glPopMatrix();

		glMatrixMode( GL_MODELVIEW );
		glPopMatrix();

		glMatrixMode( matrixMode );
	}

#ifdef NEW_VISUALIZATION_CODE
	template< typename DerivedViewableType >
	unsigned int Viewable< DerivedViewableType >::StringWidth( const Font & font , std::string str )
	{
		return glutBitmapLength( font() , reinterpret_cast< const unsigned char * >( str.c_str() ) );
	}

	template< typename DerivedViewableType >
	unsigned int Viewable< DerivedViewableType >::stringWidth( std::string str ) const
	{
		return glutBitmapLength( font() , reinterpret_cast< const unsigned char * >( str.c_str() ) );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::writeLeftString( int x , int y , std::string str ) const
	{
		WriteLeftString( x , y , font , str );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::writeRightString( int x , int y , std::string str ) const
	{
		WriteLeftString( screenWidth-x-stringWidth( str ) , y , font  ,str );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::writeCenterString( int x , int y , std::string str ) const
	{
		WriteLeftString( x-stringWidth( str )/2 , y , font , str );
	}
#else // !NEW_VISUALIZATION_CODE
	template< typename DerivedViewableType >
	int Viewable< DerivedViewableType >::StringWidth( void * font , const char * format , ... )
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		return glutBitmapLength( font , (unsigned char *) str );
	}

	template< typename DerivedViewableType >
	int Viewable< DerivedViewableType >::stringWidth( const char * format , ... ) const
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		return glutBitmapLength( font , (unsigned char *) str );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::writeLeftString( int x , int y , const char * format , ... ) const
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		WriteLeftString( x , y , font , str );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::writeRightString( int x , int y , const char * format , ... ) const
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		WriteLeftString( screenWidth-x-glutBitmapLength( font , (unsigned char *) str ) , y , font  ,str );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::writeCenterString( int x , int y , const char * format , ... ) const
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		WriteLeftString( x-glutBitmapLength( font , (unsigned char *) str )/2 , y , font  ,str );
	}
#endif // NEW_VISUALIZATION_CODE


	template< typename DerivedViewableType >
#ifdef NEW_VISUALIZATION_CODE
	void Viewable< DerivedViewableType >::saveFrameBuffer( std::string fileName , int whichBuffer )
#else // !NEW_VISUALIZATION_CODE
	void Viewable< DerivedViewableType >::saveFrameBuffer( const char* fileName , int whichBuffer )
#endif // NEW_VISUALIZATION_CODE
	{
		Pointer( float ) pixels = AllocPointer< float >( sizeof(float) * 3 * screenWidth * screenHeight );
		Pointer( unsigned char ) _pixels = AllocPointer< unsigned char >( sizeof(unsigned char) * 3 * screenWidth * screenHeight );
		glReadBuffer( whichBuffer );
		glReadPixels( 0 , 0 , screenWidth , screenHeight , GL_RGB , GL_FLOAT , pixels );
		for( int j=0 ; j<screenHeight ; j++ ) for( int i=0 ; i<screenWidth ; i++ ) for( int c=0 ; c<3 ; c++ )
		{
			int ii = int( pixels[ c + i * 3 + ( screenHeight - 1 - j ) * screenWidth * 3 ]*256 );
			if( ii<  0 ) ii =   0;
			if( ii>255 ) ii = 255;
			_pixels[ c + i * 3 + j * screenWidth * 3 ] = (unsigned char)ii;
		}
		FreePointer( pixels );
		ImageWriter< 8 >::Write( fileName , _pixels , screenWidth , screenHeight , 3 , 95 );
		FreePointer( _pixels );
	}

#ifdef NEW_VISUALIZATION_CODE
	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , std::string description , CallBackFunction callBackFunction )
	{
		callBacks.push_back( KeyboardCallBack( key , modifiers , description , callBackFunction ) );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , std::string description , std::string prompt , CallBackFunction callBackFunction )
	{
		callBacks.push_back( KeyboardCallBack( key , modifiers , description , prompt , callBackFunction ) );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , std::string description , CallBackFunction callBackFunction )
	{
		callBacks.push_back( KeyboardCallBack( key , KeyboardCallBack::Modifiers() , description , callBackFunction ) );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , std::string description , std::string prompt , CallBackFunction callBackFunction )
	{
		callBacks.push_back( KeyboardCallBack( key , KeyboardCallBack::Modifiers() , description , prompt , callBackFunction ) );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , std::string description , void (DerivedViewableType::*memberFunctionPointer)( std::string ) )
	{
		callBacks.push_back( KeyboardCallBack( key , modifiers , description , reinterpret_cast< DerivedViewableType * >( this ) , memberFunctionPointer ) );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , std::string description , std::string prompt , void (DerivedViewableType::*memberFunctionPointer)( std::string ) )
	{
		callBacks.push_back( KeyboardCallBack( key , modifiers , description , prompt , reinterpret_cast< DerivedViewableType * >( this ) , memberFunctionPointer ) );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , std::string description , void (DerivedViewableType::*memberFunctionPointer)( std::string ) )
	{
		callBacks.push_back( KeyboardCallBack( key , KeyboardCallBack::Modifiers() , description , reinterpret_cast< DerivedViewableType * >( this ) , memberFunctionPointer ) );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , std::string description , std::string prompt , void (DerivedViewableType::*memberFunctionPointer)( std::string ) )
	{
		callBacks.push_back( KeyboardCallBack( key , KeyboardCallBack::Modifiers() , description , prompt , reinterpret_cast< DerivedViewableType * >( this ) , memberFunctionPointer ) );
	}
#else // !NEW_VISUALIZATION_CODE
	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , const char *description , CallBackFunction callBackFunction )
	{
		callBacks.push_back( KeyboardCallBack( this , key , modifiers , description , callBackFunction ) );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , const char *description , const char *prompt , CallBackFunction callBackFunction )
	{
		callBacks.push_back( KeyboardCallBack( this , key , modifiers , description , prompt , callBackFunction ) );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , const char *description , CallBackFunction callBackFunction )
	{
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , key , KeyboardCallBack::Modifiers() , description , callBackFunction ) );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , const char *description , const char *prompt , CallBackFunction callBackFunction )
	{
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , key , KeyboardCallBack::Modifiers() , description , prompt , callBackFunction ) );
	}
#endif // NEW_VISUALIZATION_CODE

	template< typename DerivedViewableType >
	Viewable< DerivedViewableType >::Viewable( void )
	{
		quitFunction = []( void ){};
#ifdef NEW_VISUALIZATION_CODE
		addCallBack( KEY_ESC    , typename KeyboardCallBack::Modifiers() , ""                                  , &Viewable::quitCallBack );
		addCallBack( KEY_CTRL_C , typename KeyboardCallBack::Modifiers() , ""                                  , &Viewable::exitCallBack );
		addCallBack( 'F'        , typename KeyboardCallBack::Modifiers() , "toggle fps"                        , &Viewable::toggleFPSCallBack );
		addCallBack( 'H'        , typename KeyboardCallBack::Modifiers() , "toggle help"                       , &Viewable::toggleHelpCallBack );
		addCallBack( 'I'        , typename KeyboardCallBack::Modifiers() , "toggle info"                       , &Viewable::toggleInfoCallBack );
		addCallBack( 'i'        , typename KeyboardCallBack::Modifiers() , "save frame buffer" , "Ouput image" , &Viewable::setFrameBufferCallBack );
#else // !NEW_VISUALIZATION_CODE
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , KEY_ESC    , KeyboardCallBack::Modifiers() , "" , QuitCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , KEY_CTRL_C , KeyboardCallBack::Modifiers() , "" , ExitCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , 'F' , KeyboardCallBack::Modifiers() , "toggle fps"  , ToggleFPSCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , 'H' , KeyboardCallBack::Modifiers() , "toggle help" , ToggleHelpCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , 'I' , KeyboardCallBack::Modifiers() , "toggle info" , ToggleInfoCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , 'i' , KeyboardCallBack::Modifiers() , "save frame buffer" , "Ouput image" , SetFrameBufferCallBack ) );
		snapshotName = videoHeader = NULL;
#endif // NEW_VISUALIZATION_CODE
		flushImage = false;
		showHelp = showInfo = showFPS = true;
		_exitAfterSnapshot = _exitAfterVideo = false;
		screenWidth = screenHeight = 512;
#ifdef NEW_VISUALIZATION_CODE
		font = Font::Fonts[ Font::Type::FONT_HELVETICA_12 ];
		promptFont = Font::Fonts[ Font::Type::FONT_TIMES_ROMAN_24 ];
#else // !NEW_VISUALIZATION_CODE
		font = GLUT_BITMAP_HELVETICA_12;
		fontHeight = 12;
		promptFont = GLUT_BITMAP_TIMES_ROMAN_24;
		promptFontHeight = 24;
#endif // NEW_VISUALIZATION_CODE
		promptCallBack = NULL;
#ifdef NEW_VISUALIZATION_CODE
		originalPromptLength = 0;
#else // !NEW_VISUALIZATION_CODE
		strcpy( promptString , "" );
		promptLength = 0;
#endif // NEW_VISUALIZATION_CODE

		_lastFPSTime = Time();
		_lastFPSCount = 0;
		_fps = 0;
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::setFont( unsigned int idx )
	{
#ifdef NEW_VISUALIZATION_CODE
		if( idx>=Font::Type::FONT_COUNT ) MK_WARN( "Font index out of bounds: " , idx , " < " , Font::Type::FONT_COUNT );
		else font = Font::Fonts[idx];
#else // !NEW_VISUALIZATION_CODE
		if( idx>=Font::FontNum() ) MK_WARN( "Font index out of bounds: " , idx , " < " , Font::FontNum() );
		else font = Font::Fonts[idx].font , fontHeight = Font::Fonts[idx].fontHeight;
#endif // NEW_VISUALIZATION_CODE
	}


	template< typename DerivedViewableType >
#ifdef NEW_VISUALIZATION_CODE
	void Viewable< DerivedViewableType >::setSnapshot( std::string sName , bool exitAfterSnapshot )
#else // !NEW_VISUALIZATION_CODE
	void Viewable< DerivedViewableType >::setSnapshot( const char* sName , bool exitAfterSnapshot )
#endif // NEW_VISUALIZATION_CODE
	{
		_exitAfterSnapshot = exitAfterSnapshot;
#ifdef NEW_VISUALIZATION_CODE
		snapshotName = sName;
#else // !NEW_VISUALIZATION_CODE
		snapshotName = new char[ strlen( sName ) + 1 ];
		strcpy( snapshotName , sName );
#endif // NEW_VISUALIZATION_CODE
		showHelp = showInfo = showFPS = false;
		flushImage = true;
	}

	template< typename DerivedViewableType >
#ifdef NEW_VISUALIZATION_CODE
	void Viewable< DerivedViewableType >::setVideo( std::string vHeader , int frames , bool exitAfterVideo )
#else // !NEW_VISUALIZATION_CODE
	void Viewable< DerivedViewableType >::setVideo( const char* vHeader , int frames , bool exitAfterVideo )
#endif // NEW_VISUALIZATION_CODE
	{
		_exitAfterVideo = exitAfterVideo;
#ifdef NEW_VISUALIZATION_CODE
		videoHeader = vHeader;
#else // !NEW_VISUALIZATION_CODE
		videoHeader = new char[ strlen( vHeader ) + 1 ];
		strcpy( videoHeader , vHeader );
#endif // NEW_VISUALIZATION_CODE
		_currentFrame = 0;
		_totalFrames = frames;
	}

	template< typename DerivedViewableType >
#ifdef NEW_VISUALIZATION_CODE
	void Viewable< DerivedViewableType >::setFrameBufferCallBack( std::string prompt )
#else // !NEW_VISUALIZATION_CODE
	void Viewable< DerivedViewableType >::SetFrameBufferCallBack( DerivedViewableType* v , const char* prompt )
#endif // NEW_VISUALIZATION_CODE
	{
#ifdef NEW_VISUALIZATION_CODE
		if( prompt.size() )
		{
			snapshotName = prompt;
			flushImage = true;
		}
#else // !NEW_VISUALIZATION_CODE
		if( prompt )
		{
			v->snapshotName = new char[ strlen(prompt)+1 ];
			strcpy( v->snapshotName , prompt );
			v->flushImage = true;
		}
#endif // NEW_VISUALIZATION_CODE
	}
}
#endif // VISUALIZATION_INCLUDED