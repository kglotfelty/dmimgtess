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
 * FILE NAME:  vtpdetect.h
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * This file contains the #defines (macros), data structure
 * definitions, and subroutine prototyptes for the vtpdetect program.
 * 
 *
 * REVISION HISTORY:
     0.5.1          1998 Jul 1   Cleaned up move of AOTS library
 *
H***************************************************************** */

#include "vtpdetect.h"


/*
  +---------------------------------------
  +
  + This interface routine simply copies the events from
  + one data format to another and runs the triangulation
  + routine.
  +
  +---------------------------------------
  */
void doTriangulation(
		     struct Event   *events,      /*i: input events       */
		     long            numEvts,     /*i: no. of events      */
		     int             excludeBack, /*i: exclude bkg events?*/
		     struct TempTri *nodes,       /*o: triangle nodes     */
		     long           *ntri,        /*o: no. triangles      */
		     char           *delfile      /*o: debug file name    */
		     )
{

  struct Site *sites; /* native triangulation data structure */
  long nsites;        /* number of sites */
  long i; /* loop variable */


  sites=(struct Site *)calloc(numEvts,sizeof(*sites));
  nsites = 0;
  for (i=0;i<numEvts;i++)
    {
      if (( !excludeBack ) || ( ( events[i].src_no != BACKGROUND   ) &&
				( events[i].src_no != TRIMMED_EVENT) ) )
	{
	  sites[nsites].coord.x = events[i].coord.x;
	  sites[nsites].coord.y = events[i].coord.y;
	  sites[nsites].sitenbr = i;
	  nsites += 1;
	}
    }
  
  sites=(struct Site *)realloc(sites, nsites*sizeof(*sites));

  /* do all the real work */
  voroniInterface( sites, nsites, nodes, ntri, delfile );

  free(sites);

  return;

}
