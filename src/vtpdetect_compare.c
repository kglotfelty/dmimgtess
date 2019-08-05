/*                                                                
**  Copyright (C) 1997,1998,2007  Smithsonian Astrophysical Observatory 
*/                                                                

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */

/* H*****************************************************************
 *
 * FILE NAME:  vtpdetect_compare.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * This file contians routines used by vtpdetect to sort various
 * data structures using the "qsort" routine (from stdlib).
 * 
 *
 * REVISION HISTORY:
 *
 * Ref. No.         Date
   ----------       -----
     0.1            1997 Dec 30  Initial putback of working copy

 *
H***************************************************************** */
#include <stdlib.h>
#include "vtpdetect.h"
#include <float.h>


/*
  + -------------------------------------------------
  +
  + This subrountine is used by qsort to sort the events by Y and X
  +
  +-------------------------------------------------
  */
extern int sortYX(
		  const void *ii,      /* void pointer to an event structure */
		  const void *jj       /* void pointer to an event structure */
		  )
{
  struct Event *s1,*s2;
  int           retval = 0;   /* return value */


  s1 = (struct Event *) ii;    /* cast void pointers to event structures */
  s2 = (struct Event *) jj;
  
  
  if ( ( fabs( (s1->coord.y - s2->coord.y)/s1->coord.y) < 2*FLT_EPSILON) &&
       ( fabs( (s1->coord.x - s2->coord.x)/s1->coord.x) < 2*FLT_EPSILON)    )
    {
      retval = 0;
    }
  else if ( s1->coord.y < s2->coord.y ) 
    {
      retval = -1;    /* y1 smaller than y2            */
    }
  else if ( s1->coord.y > s2->coord.y ) 
    {
      retval = 1;     /* y1 bigger than y2             */
    }
  else if ( s1->coord.x < s2->coord.x ) 
    {
      retval = -1;    /* y1 = y2 && x1 smaller than x2 */
    }
  else if ( s1->coord.x > s2->coord.x ) 
    {
      retval = 1;     /* y1 = y2 && x1 bigger than x2  */
    }
  else
    {
      retval = 0;
    }

  return( retval );               /* =0 only if event (x,y)'s are equal */
}





/*
  +----------------------------------------
  +
  +This routine is used  to sort flux values using qsort
  +
  +-----------------------------------------
  */
extern int sortLo2Hi(
		     const void *ii,         /* void pointer to double */
		     const void *jj          /* void pointer to double */
		     )
{
  double *xx, *yy;
  int    retval=0;
  
  xx = (double *) ii;              /* cast voids to doubles */
  yy = (double *) jj;
  
  if ( *xx > *yy ) 
    {
      retval = 1;         /* xx is bigger than yy  */
    }
  else if ( *xx < *yy )
    {
      retval = -1;        /* xx is smaller than yy */
    }
  else
    {
      retval = 0;
    }


  return ( retval );        /* =0 only if equal */
}





/*
  + -------------------------------------------------
  +
  + This subrountine is used by qsort to sort the events by 
  + source number, Y, and X
  +
  +-------------------------------------------------
  */
extern int sortIdYX(
		 const void *ii,       /* void pointer to an event structure */
		 const void *jj        /* void pointer to an event structure */
		 )
{
  struct Event *s1,*s2;
  int           retval = 0;

  s1 = (struct Event *) ii;    /* cast void pointers to event structures */
  s2 = (struct Event *) jj;
  

  if ( s1->src_no < s2->src_no ) 
    {
      retval = -1;     /* sn1 is smaller than sn2                      */
    }
  else if ( s1->src_no > s2->src_no ) 
    {
      retval = 1;      /* sn1 is bigger than sn2                       */
    } 
  else if ( ( fabs( (s1->coord.y - s2->coord.y)/s1->coord.y) < 2*FLT_EPSILON) &&
	    ( fabs( (s1->coord.x - s2->coord.x)/s1->coord.x) < 2*FLT_EPSILON)    )
    {
      retval = 0;
    }
  else if ( s1->coord.y < s2->coord.y ) 
    {
      retval = -1;     /* sn1 == sn2 and y1 smaller than y2            */
    }
  else if ( s1->coord.y > s2->coord.y ) 
    {
      retval = 1;      /* sn1 == sn2 and y1 bigger than y2             */
    }
  else if ( s1->coord.x < s2->coord.x ) 
    {
      retval = -1;     /* sn1 == sn2 && y1 == y2 && x1 smaller than x2 */
    }
  else if ( s1->coord.x > s2->coord.x ) 
    {
      retval = 1;      /* sn1 == sn2 && y1 == y2 && x1 bigger than x2  */
    }

  return( retval );               /* =0 only if event (x,y)'s are equal */
}


