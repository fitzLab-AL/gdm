//
//
//	Message.cpp
//
//
#include "stdafx.h"
#include <stdio.h>
#include "Message.h"

// Modified by DNL:
//#pragma warning( disable : 4996 )

//////////////////////////////////////////////////////////////////////////////////////
//
// Message functions for string arguments
//
void Message( char *message )
{
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)message, (LPCSTR)"", MB_OK );
#endif
}


void Message( char *message, char *title )
{
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)message, (LPCSTR)title, MB_OK );
#endif
}


//////////////////////////////////////////////////////////////////////////////////////
//
// Message functions for 64 bit integer arguments
//
void Message( long long value )
{
	char *p = new char [256];
	// Modified by DNL:
	//sprintf( p, "%1.0lf", (double)value );
	snprintf( p, sizeof(p)*2, "%1.0lf", (double)value );
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)"", MB_OK );
#endif
	if ( p ) delete[] p;
}


void Message( char *name, long long value )
{
	char *p = new char [256];
	// Modified by DNL:
	//sprintf( p, "%s = %1.0lf", name, (double)value );
	snprintf( p, sizeof(p)*2, "%s = %1.0lf", name, (double)value );
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)"", MB_OK );
#endif
	if ( p ) delete[] p;
}


void Message( long value, char *title )
{
	char *p = new char [256];
	// Modified by DNL:
	//sprintf( p, "%1.0lf", (double)value );
	snprintf( p, sizeof(p)*2, "%1.0lf", (double)value );
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)title, MB_OK );
#endif
	if ( p ) delete[] p;
}


void Message( char *name, long long value, char *title )
{
	char *p = new char [256];
	// Modified by DNL:
	//sprintf( p, "%s = %1.0lf", name, (double)value );
	snprintf( p, sizeof(p)*2, "%s = %1.0lf", name, (double)value );
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)title, MB_OK );
#endif
	if ( p ) delete[] p;
}


//////////////////////////////////////////////////////////////////////////////////////
//
// Message functions for integer arguments
//
void Message( int value )
{
	char *p = new char [256];
	snprintf( p, sizeof(p)*2, "%d", value );
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)"", MB_OK );
#endif
	if ( p ) delete[] p;
}


void Message( char *name, int value )
{
	char *p = new char [256];
	snprintf( p, sizeof(p)*2, "%s = %d", name, value );
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)"", MB_OK );
#endif
	if ( p ) delete[] p;
}


void Message( int value, char *title )
{
	char *p = new char [256];
	snprintf( p, sizeof(p)*2, "%d", value );
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)title, MB_OK );
#endif
	if ( p ) delete[] p;
}


void Message( char *name, int value, char *title )
{
	char *p = new char [256];
	snprintf( p, sizeof(p)*2,"%s = %d", name, value );
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)title, MB_OK );
#endif
	if ( p ) delete[] p;
}

//////////////////////////////////////////////////////////////////////////////////////
//
// Message functions for float arguments
//
void Message( float value )
{
	char *p = new char [256];
	snprintf( p, sizeof(p)*2, "%lf", value );
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)"", MB_OK );
#endif
	if ( p ) delete[] p;
}


void Message( char *name, float value )
{
	char *p = new char [256];
	snprintf( p, sizeof(p)*2, "%s = %lf", name, value );

#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)"", MB_OK );
#endif
	if ( p ) delete[] p;
}



void Message( float value, char *title )
{
	char *p = new char [256];
	snprintf( p, sizeof(p)*2, "%lf", value );
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)title, MB_OK );
#endif
	if ( p ) delete[] p;
}



void Message( char *name, float value, char *title )
{
	char *p = new char [256];
	snprintf( p, sizeof(p)*2, "%s = %lf", name, value );

#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)title, MB_OK );
#endif
	if ( p ) delete[] p;
}


//////////////////////////////////////////////////////////////////////////////////////
//
// Message functions for double arguments
//
void Message( double value )
{
	char *p = new char [256];
	// Modified by DNL:
	//sprintf( p, "%lf", value );
	snprintf( p, sizeof(p)*2, "%lf", value );
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)"", MB_OK );
#endif
	if ( p ) delete[] p;
}


void Message( char *name, double value )
{
	char *p = new char [256];
	// Modified by DNL:
	//sprintf( p, "%s = %lf", name, value );
	snprintf( p, sizeof(p)*2, "%s = %lf", name, value );

#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)"", MB_OK );
#endif
	if ( p ) delete[] p;
}



void Message( double value, char *title )
{
	char *p = new char [256];
	// Modified by DNL:
	//sprintf( p, "%lf", value );
	snprintf( p, sizeof(p)*2, "%lf", value );
#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)title, MB_OK );
#endif
	if ( p ) delete[] p;
}



void Message( char *name, double value, char *title )
{
	char *p = new char [256];
	// Modified by DNL:
	//sprintf( p, "%s = %lf", name, value );
	snprintf( p, sizeof(p)*2, "%s = %lf", name, value );

#ifdef _WIN32
	MessageBox( NULL, (LPCSTR)p, (LPCSTR)title, MB_OK );
#endif
	if ( p ) delete[] p;
}
