#include <cdvmh_helpers.h>

#include "wtime.h"
#include <stdlib.h>

/*  Prototype  */
void wtime( double * );


/*****************************************************************/
/******         E  L  A  P  S  E  D  _  T  I  M  E          ******/
/*****************************************************************/
static double elapsed_time( void )
{
    double t;

    wtime( &t );
    return( t );
}


static double start[64], elapsed[64];

/*****************************************************************/
/******            T  I  M  E  R  _  C  L  E  A  R          ******/
/*****************************************************************/
void timer_clear( int n )
{
    elapsed[n] = 0.0;
}


/*****************************************************************/
/******            T  I  M  E  R  _  S  T  A  R  T          ******/
/*****************************************************************/
void timer_start( int n )
{
    start[n] = elapsed_time();
}


/*****************************************************************/
/******            T  I  M  E  R  _  S  T  O  P             ******/
/*****************************************************************/
void timer_stop( int n )
{
    double t, now;

    now = elapsed_time();
    t = now - start[n];
    elapsed[n] += t;

}


/*****************************************************************/
/******            T  I  M  E  R  _  R  E  A  D             ******/
/*****************************************************************/
double timer_read( int n )
{
    return( elapsed[n] );
}



void initCdvmhGlobals_c_timers_315728945() {
    dvmh_data_enter_C((const void *)start, sizeof(start));
    dvmh_data_enter_C((const void *)elapsed, sizeof(elapsed));
}

void finishCdvmhGlobals_c_timers_2141340040() {
    dvmh_data_exit_C((const void *)elapsed, 0);
    dvmh_data_exit_C((const void *)start, 0);
}

