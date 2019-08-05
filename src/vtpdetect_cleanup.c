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
 * FILE NAME:  vtpdetect_cleanup.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * This file contains routine needed by the vtpdetect algorithm to
 * perform some memory de-allocation routine and edge removal routines.
 *
 * REVISION HISTORY:
 *
 * Ref. No.         Date
   ----------       -----
     0.1            1997 Dec 30  Initial putback of working copy

 *
H***************************************************************** */






#include <stdio.h>
#include <stdlib.h>
#include "vtpdetect.h"



/*
  +-------------------------------------
  +
  +  This routine takes in the list of sources
  +  and systematically goes thru and deletes them.
  +
  +  This routine is used between iterations.  Once the
  +  threshold of just the background events is found,
  +  the source list is purged (re-initialiazed to zero)
  +  and the percolation is repeated.
  +
  +-------------------------------------
  */
void eraseSources( 
		  struct Source *srcs     /* i/o: source list */
		  )
{

  struct Source  *atsrc, *nextsrc;          /* ptrs to sources to delete    */
  struct EvtList *atEvtNode, *nextEvtNode;  /* ptrs to event nodes to delete*/
  struct Cell    *atCell, *nextCell;
  struct Points  *atPos,  *nextPos;
  


  outputDebug( DEBUG_L3, 
	       "cleaning old source list.\n", 
	       BEFORE );                    /* output debugging statement */


  
  atsrc = srcs;        /* Start at the first one */

  while ( atsrc != NULL )     
    {

      /* go thru all srcs */

      nextsrc   = atsrc->next;
      atEvtNode = atsrc->evts;

      while ( atEvtNode != NULL )
	{

	  /* go thru all events that the source includes */

	  nextEvtNode     = atEvtNode->next;
	  atEvtNode->at   = NULL;
	  atEvtNode->next = NULL;

	  free( atEvtNode );  /* free node memory */

	  atEvtNode = nextEvtNode; /* get next node */

	} /* end while (atEvtNode... */


      atsrc->next = NULL;
      atsrc->evts = NULL;

      /* free super cells */
      atCell = atsrc->cell;
      while ( atCell != NULL )
	{
	  nextCell = atCell->next;

	  /* free pixel positions */
	  atPos = atCell->at;
	  while ( atPos != NULL )
	    {
	      nextPos     = atPos->next;
	      atPos->next = NULL;

	      free( atPos );
	      atPos = nextPos;
	    } 

	  atCell->at   = NULL;
	  atCell->next = NULL;
	  
	  free( atCell );
	  atCell = nextCell;
	}


      free( atsrc );  /* free source */

      atsrc = nextsrc; /* get next source */

    } /* end while (atsrc ... */


  outputDebug( DEBUG_L3, 
	       "cleaning old source list.\n", 
	       AFTER );                 /* output debugging statement */

  return;
}




/*
  +-------------------------------------
  +
  + This routine goes thru and identifies events
  + that are close to the edge of the field.  The
  + triangulation is typically un-reliable at/near
  + the edge.
  +
  + The idea is to use the source number to repetively
  + "shrink" the region of interest.  Events that
  + are at the very edge were identifed in the
  + calculateArea routine and given a src_no = AT_EDGE.
  +
  + This routine will build on this.  Events that are
  + part of a triangle with the AT_EDGE events will get
  + a src_no = AT_EDGE - 1.
  +
  + Events that are part of a triangle with those events
  + get a src_no = AT_ECGE - 2.  
  +
  + And so on, for "level" iterations.
  +
  +-------------------------------------
  */
void removeEdges( 
		 struct Event *events,    /* i: event array               */
		 long          numEvts,   /* i: total number of events    */
		 int           level,     /* i: no. of events to the edge */
		 long         *removed    /* o: no. of events removed     */
		 )
{
  long            ii;                /* loop variable */
  int             atlevel;           /* current edge level */
  int             setlevel;          /* setting edge level */
  int             jj,kk;             /* loops */
  struct TriList *atTriNode;         /* pointer to current triangle node */


  outputDebug( DEBUG_L3,
	       "removing events at edge.\n", 
	       BEFORE );                /* output debugging message */

  
  atlevel  = AT_EDGE;       /* 0th level is at the very edge */
  setlevel = atlevel-1;     /* 1st level is AT_EDGE -1 */


  /* for each level */

  for ( jj=0; jj<level; jj++ )                     
    {

      /* for each event */

      for ( ii=0; ii<numEvts; ii++ )      
	{

	  /* if the event is at the current inner most level
	     then we want to examine all the triangles that
	     have it as a vertex */

	  if ( events[ii].src_no == atlevel )
	    {

	      atTriNode = events[ii].tris;  /* first node */

	      while ( atTriNode != NULL )
		{

		  /* for each triangle */

		  for ( kk=0; kk<3; kk++ )
		    {

		      /* for each vertex */
		      
		      if ( atTriNode->at->verts[kk]->src_no == BACKGROUND )
			{
			  atTriNode->at->verts[kk]->src_no = setlevel;
			  atTriNode->at->verts[kk]->area   = 0.0;
			  (*removed) +=1;
			}
		      
		    } /* end for each triangle, for (kk=0;kk<3...*/
		  

		  atTriNode = atTriNode->next;  /* get next triangle */


		} /* end for each triangle, while (atTriNode .... */

	    } /* end if (events[ii].src_no ... */

	} /* end for each event, for (ii=0;ii<n */


      atlevel   = setlevel;  /* set current inner most level to level
				just found */
      setlevel -= 1;         /* set current set level to next level */
      
    } /* end for each level, for (jj=0,jj<level... */


  outputDebug( DEBUG_L3, 
	       "removing events at edge.\n", 
	       AFTER );                /* output debugging message */


  return;

}




/*
  +-----------------------------------------
  +
  + Remove events from memory
  +
  +-----------------------------------------
  */
void eraseEvents( 
		struct Event *evts,   /* i/o:  events to remove */
		long          numEvts /* i:    number of events */
		)
{
  
  struct TriList *atTriNode, *nextTriNode;  /* ptrs to triangle nodes */
  struct Event   *atEvt;                    /* ptr to current evnet   */
  long            ii;                       /* generic loop variable  */

  
  outputDebug( DEBUG_L3, 
	       "erasing events.\n", 
	       BEFORE );                    /* output debugging statement */

  
  for (ii=0; ii<numEvts; ii++ )
    {
      atEvt     = evts+ii;
      atTriNode = atEvt->tris;
      
      while ( atTriNode != NULL )
	{
	  /* go thru all triangle nodes */
	  
	  nextTriNode     = atTriNode->next;
	  atTriNode->at   = NULL;
	  atTriNode->next = NULL;
	  
	  free(atTriNode);   /* free triangle node */
	  
	  atTriNode = nextTriNode;
	  
	} /* end loop over TriNodes, while(... ) */


      /* Super cells */
      atTriNode = atEvt->superTris;
      while ( atTriNode != NULL )
	{
	  /* go thru all triangle nodes */
	  
	  nextTriNode     = atTriNode->next;
	  atTriNode->at   = NULL;
	  atTriNode->next = NULL;
	  
	  free(atTriNode);   /* free triangle node */
	  
	  atTriNode = nextTriNode;
	  
	}



    } /* end loop over events, for(ii...) */

  
  outputDebug( DEBUG_L3, 
	       "erasing events.\n", 
	       AFTER );                 /* output debugging statement */

  return;
}

