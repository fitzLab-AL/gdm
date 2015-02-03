//
//
//	Message.h
//
//

#ifndef	__MESSAGE_H__
#define	__MESSAGE_H__


void Message( char *message );
void Message( char *message, char *title );


void Message( long long value );
void Message( char *name, long long value );
void Message( long long value, char *title );
void Message( char *name, long long value, char *title );


void Message( int value );
void Message( char *name, int value );
void Message( int value, char *title );
void Message( char *name, int value, char *title );


void Message( float value );
void Message( char *name, float value );
void Message( float value, char *title );
void Message( char *name, float value, char *title );


void Message( double value );
void Message( char *name, double value );
void Message( double value, char *title );
void Message( char *name, double value, char *title );


#endif //__MESSAGE_H__
