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
 * FILE NAME:  vtpdetect_aots.h
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * This file contains the #defines (macros), data structure
 * definitions, and subroutine prototyptes for the vtpdetect that
 * are needed to interface with the AOTS triangulation software.
 *
 * REVISION HISTORY:
 *
 * Ref. No.         Date
   ----------       -----
   0.5               23June 1998

 *
H***************************************************************** */


#ifndef VTPAOTS_H
#define VTPAOTS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* common data structures to both */


extern int   GlobalError;           /* Flag to hold error bits */

extern int     GlobalDebugAuxFile;  /* flag to hold error bits if debugging 
				       files could note be opened */
/* -------------------  DEBUGGING INFO  ------------------ */

extern int     GlobalDebugLevel;   /* Make debug available to all funs */     


#define DEBUG_L0                0  /* no debugging info                 */
#define DEBUG_L1                1  /* echo parameters                   */
#define DEBUG_L2                2  /* report during each iteration      */
#define DEBUG_L3                3  /* report during most function calls */
#define DEBUG_L4                4  /* report when searching for sources */
#define DEBUG_L5                5  /* report during percolation         */

             /* NOTE:  use levels 4 & 5 with extreme caution.  
                The output could be 100's of pages */

extern FILE *GlobalDebugFile;

/* ------------------ BITS -------------------- */

#define BIT00   0x0001
#define BIT01   0x0002        
#define BIT02   0x0004
#define BIT03   0x0008
#define BIT04   0x0010
#define BIT05   0x0020
#define BIT06   0x0040
#define BIT07   0x0080
#define BIT08   0x0100
#define BIT09   0x0200
#define BIT10   0x0400
#define BIT11   0x0800
#define BIT12   0x1000
#define BIT13   0x2000
#define BIT14   0x4000
#define BIT15   0x8000

/* ----------------- Data structures ----------------- */



struct Point	{
float x,y;
};


/* structure used both for sites and for vertices  */
struct Site	{
struct	Point	coord;
long		sitenbr;
long		refcnt;
};


struct TempTri {
  long id1, id2, id3;
};




extern void voroniInterface(
			    struct Site *sites,
			    long n,
			    struct TempTri *nodes,
			    long *ntri,
			    char *delfile
			    );



#endif
