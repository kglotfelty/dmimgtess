/*                                                                
**  Copyright (C) 1997,1998,2007,2011  Smithsonian Astrophysical Observatory 
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
 * FILE NAME:  vtpdetect_voroni.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * This file contians routines used by vtpdetect to interface with
 * the adopted off the shelf code to compute the Delaunay triangulation
 * and then compute the area of the Voronoi cells.
 * 
 * 
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
#include <math.h>
#include "vtpdetect.h"
#include <float.h>


int inPoly( struct Pt p, struct Pt *points, long  npoints );
void resetId( struct Event *atEvt );
long countBack( struct Event *atEvt, struct Pt *sides, long *nSides, 
		long nCell, FILE *bkfile );




/*
  +-------------------------------------------------
  +
  + This subrountine takes the three events th make up the
  + vertices of a Delauny triangle and creates the triangle structure.
  + Essentailly this entials establishing pointers from the events
  + to the triangle and from the triangle to the events.   This data
  + is needed to a) compute the area of the Voroni cells and b)
  + perform the percolation.
  +
  +-------------------------------------------------
  */
void makeTriangle(
		  struct Event     *events,      /* i: event list             */
		  long              numEvts,     /* i: no. of events          */
		  struct TempTri   *nodes,       /* i: list of verticies      */
		  long              ntri,        /* i: no. of triangles       */
		  short             bORa,        /* i: before or after src. det*/
		  struct Triangle **allTriangles /* o: first triangle in list */
		  )
{
  
  struct Event    *evtarr[3];         /* array holding 3 events */
  double           aa,ab,ac;          /* angles to 3 events     */
  double           angarr[3];         /* angles stored in array */

  double           pi  = acos(-1.0);
  double           tpi = 2*pi;

  int              idx[3];            /* index array   */

  int              ii,jj;             /* loop variables */
  long             kk;

  struct TriList **atTriNode;         /* pointers to triangle lists */
  struct TriList  *newTriNode;

  double            xc, yc,r,xm,ym;    /* circle parameters */

  struct Triangle *newTriangle;       /* pointer to new triangle */
  
  struct Event    *triA, *triB, *triC;/* 3 events making up triangle */

  long             triid;             /* triangle number */

  struct Angles   *angStruct;         /* pointer to angles   */

  unsigned long    p1, p2;            /* offsets to triangle numbers */





  *allTriangles = (struct Triangle *) calloc( ntri, sizeof( struct Triangle ) );
  if ( allTriangles == NULL )
    {
      err_msg("ERROR: could not allocate memory for the triangles.\n");
      GlobalError = ( GlobalError | BIT00 );
    }

  angStruct = (struct Angles *) calloc( ntri, sizeof( struct Angles ) );
  if ( angStruct == NULL )
    {
      err_msg("ERROR: could not allocate memory for the triangle angles.\n");
      GlobalError = ( GlobalError | BIT00 );
    }

  if ( !GlobalError )
    {
      for ( kk=0; kk<ntri; kk++)       /* for each triangle */
	{
	  
	  triA  = &events[nodes[kk].id1];  /* assign variables to each vert */
	  triB  = &events[nodes[kk].id2];
	  triC  = &events[nodes[kk].id3];
	  triid = kk+1;
	  
	  newTriangle = (*allTriangles)+kk;
	  
	  
	  /* store three event structs in an array for easy manipulation */
	  
	  evtarr[0] = triA;   
	  evtarr[1] = triB;
	  evtarr[2] = triC;
	  
	  /* need the center(and radius) of the circles that 
	     circumscribes these points */
      
	  calculateCircle( triA, triB, triC, &xc, &yc, &r, &xm, &ym );
	  
	  /* need the angle each event makes with the center.*/
	  
	  /* This was change from the center of the circle, to the center
	     of mass of the triangle!  This is because of degenericies
	     introduced by quantitizing pixel locations.  If 4 events
	     form a perfect square, two triangles will be constructed --
	     both with the same center => will make same angle => will 
	     not be properly.  The center of mass however is unique to
	     each triangle 
	     */
	  
	  aa = atan2( triA->coord.y-ym, triA->coord.x-xm );
	  if (aa < 0.0 ) aa += tpi;
	  ab = atan2( triB->coord.y-ym, triB->coord.x-xm );
	  if (ab < 0.0 ) ab += tpi;
	  ac = atan2( triC->coord.y-ym, triC->coord.x-xm );
	  if (ac < 0.0 ) ac += tpi;
	  
	  angarr[0] = aa;
	  angarr[1] = ab;
	  angarr[2] = ac;
	  
	  /* Now we want to order the events based on the angle they make with
	     the center.  Similiarly, the triangles will be ordered by the same 
	     angle. We use the idx array to store the index of the ordered 
	     angles */
	  
	  
	  if (( aa < ab ) && ( aa < ac ) )
	    {
	      idx[0] = 0;
	      if (ab < ac) 
		{
		  idx[1] = 1;  /* arrange A, B, C */
		  idx[2] = 2;  
		}
	      else
		{
		  idx[1] = 2;  /* arrange A, C, B */
		  idx[2] = 1;
		}
	    }
	  else if ((ab < aa) && (ab < ac))
	    {
	      idx[0] = 1;
	      if (aa < ac)
		{
		  idx[1] = 0;  /* arrange B, A, C */
		  idx[2] = 2;
		}
	      else
		{
		  idx[1] = 2;  /* arrange B, C, A */
		  idx[2] = 0;
		}
	    }
	  else
	    {
	      idx[0] = 2;
	      if (aa < ab)
		{
		  idx[1] = 0;  /* arrange C, A, B */
		  idx[2] = 1;
		}
	      else
		{
		  idx[1] = 1;  /* arrange C, B, A */
		  idx[2] = 0;
		}
	    }
	  
	  /* populate the triangle structure */
	  
	  newTriangle->xc = xc;
	  newTriangle->yc = yc;
	  
	  for ( ii=0; ii<3; ii++) 
	    { 
	      newTriangle->verts[ii] = evtarr[idx[ii]];
	      
	      angStruct[kk].angs[ii] = angarr[idx[ii]];
	      

	      /* super cell update.  */
	      if ( bORa == BEFORE )
		{
		  atTriNode = &(newTriangle)->verts[ii]->tris;
		}
	      else
		{
		  atTriNode = &(newTriangle)->verts[ii]->superTris;
		}
	      
	      
	      
	      newTriNode = (struct TriList *) calloc( 1, sizeof(struct TriList) );
	      if (newTriNode == NULL )
		{
		  err_msg("ERROR: could not allocate memory for triangle node.\n");
		  GlobalError = ( GlobalError | BIT00 );
		  break;
		}
	      
	      newTriNode->at = newTriangle;
	      
	      /* Want to order the triangles in the triangle list,
		 there check to see a) if this is the first triangle
		 in the list or b) find where it fits into list and
		 insert it.
		 */
	      
	      if ( *atTriNode == NULL)    /* first triange for this event */
		{

		  if ( bORa == BEFORE )
		    {
		      newTriangle->verts[ii]->tris = newTriNode;
		    }
		  else
		    {
		      newTriangle->verts[ii]->superTris = newTriNode;
		    }

		}
	      else
		{
		  do
		    {
		      /* find which vertice index to compare to */
		      
		      for ( jj=0; jj<3; jj++)
			{
			  if ((*atTriNode)->at->verts[jj]  == 
			      newTriangle->verts[ii] )
			    {
			      break;  /* break leaving pointer where == */
			    }
			} 
		      
		      /* compare angles, if the angle is >, then
			 break out of the do ... while loop, leaving the
			 pointer to the triangle where it was found.
			 */
		      
		      
		      p1 = newTriNode->at   - *allTriangles;
		      p2 = (*atTriNode)->at - *allTriangles;
		      
		      
		      if ( angStruct[p1].angs[ii] > angStruct[p2].angs[jj] )
			{  
			  break;  /* break when angle is > */
			} 
		      else
			{ 
			  atTriNode = &(*atTriNode)->next; /* get next triangle in list*/ 
			}
		      
		    }  while ( (*atTriNode) != NULL ); /* until no more triangles */
		  
		  /* insert the triangle at the last place in the loop:
		     either as the first element, anywhere in middle, or end 
		     */
		  
		  newTriNode->next = (*atTriNode);
		  (*atTriNode)     = newTriNode;
		  
		} /* end if (*atTriNode == ... */
	      
	    } /* end for (ii=0;ii<3;... */
	  
	  if (GlobalError) break;

	} /* end for (k=0; k<ntri... */

      free(angStruct);
      
    } /* end if (!Global ... */
  
  return;

  
}




/*
  +-------------------------------------------------
  +
  + This routine computes the area of the Voroni
  + cells created from the Delauny triangles.
  + This routien also identifies events that are at
  + the field boundary and assigns them a scr_no = AT_EDGE.
  +
  + The area is computed at the sum x(i)*y(i+1) - x(i+1)*y(i)
  + for i=0 to n-2, and then adding x(0)*y(n-1) - x(n-1)*y(0).
  + We can use this formula since the events have been orgainzied
  + by angle!
  +
  +-------------------------------------------------
  */

void calculateArea (
		    struct Event  *events,    /* i: array with all the events  */
		    long           numEvts,   /* i: number of events           */
		    char          *vorfile    /* i: file name for debug output */
		    )
{
  int             jj,kk,flag;   /* loop flags */
  long            ii;

  struct TriList *atTriNode1;  /* pointers to two triangles */
  struct TriList *atTriNode2;
  double          area;

  long     n0exp =0;

  FILE *filen = NULL;  /* debug file number */



  if ( GlobalDebugLevel >= DEBUG_L3 )
    {

      if ( (strcmp( vorfile, "none" ) != 0 ) &&
	   (strcmp( vorfile, ""     ) != 0 ))
	{
	  filen = fopen( vorfile, "w" );
	  if ( filen == NULL )
	    {
	      err_msg(
		       "WARNING: could not open Voronoi regions file: %s\n", vorfile);
	      GlobalDebugAuxFile = ( GlobalDebugAuxFile | BIT01 );
	    }
	}
    }

  

  /* Go thru all events */
  
  for ( ii=0; ii<numEvts; ii++ )
    {
      area       = 0.0;

      atTriNode1 = events[ii].tris ;
      atTriNode2 = atTriNode1->next;


      /* check exposure time */
      if ( fabs( events[ii].exp ) < FLT_EPSILON )
	{
	  events[ii].area   = 0.0;
	  events[ii].src_no = AT_EDGE;
	  n0exp += 1;

	  if ( GlobalDebugLevel >= DEBUG_L5 )
	    {
	      err_msg("WARNING: pixel with zero exposure found at (%f,%f).\n.",
		       events[ii].coord.x, events[ii].coord.y );
	    }


	  continue;
	}
      
      /* There must be at least 3 triangles to compute the area.  If 2nd in
	 list is a NULL => not a triangle => this event is on the edge
	 of the field and should be tagged 
	 */

      if ( atTriNode2 == NULL )
	{
	  events[ii].area   = 0.0;
	  events[ii].src_no = AT_EDGE;   /* src_no AT_EDGE => edge */
	}
      else
	{
	  /*  Check to see if there is a 3rd triangle.  If not,
	      then this must be on the edge of the field and tagged */
	  
	  if ( atTriNode2->next == NULL )   
	    {
	      events[ii].area   = 0.0;
	      events[ii].src_no = AT_EDGE;
	    }
	  else
	    {
	      
	      /* repeat for all triangles in the triangle list */
	      
	      while ( atTriNode2 != NULL )
		{
		  
		  /* Need to check to make sure that the triangles
		     are connected at, atleast 2 vertices */
		  
		  flag = 0;
		  for ( jj=0; jj<3; jj++)
		    {
		      for ( kk=0; kk<3; kk++)
			{
			  if ( atTriNode1->at->verts[jj] == atTriNode2->at->verts[kk] )
			    {
			      flag++;

			    } /* endif */

			} /* end for (kk..) */

		    }  /* end for (jj..) */
		  
		  
		  /* If it could not find two verticies in common,
		     then the event is on the edge of the field.  This is
		     because the events were orgainzed by angle.  Going
		     counter clockwise, if an event is in the interior of
		     the field, we should be able to go from 
		     (a,b) -> (b,c) -> (c,d) ->... etc.. -> (z,a).  
		     If we go  (a,b) -> (c,d) , then the b,c edge is on the
		     edge of the field. */
		  
		  
		  if (flag < 2)
		    {
		      atTriNode2        = NULL;
		      area              = 0.0;
		      events[ii].src_no = AT_EDGE;
		    }
		  else
		    {
		      /* Compute area! */
		      
		      area  += (
			        ( ((double) atTriNode1->at->xc)*atTriNode2->at->yc ) - 
			        ( ((double) atTriNode2->at->xc)*atTriNode1->at->yc )   );
		      
		      atTriNode1 = atTriNode2;
		      atTriNode2 = atTriNode1->next;
		      
		    } /* end if (flag < 2 ) */
		  
		}  /* end while ( atTriNode2 != ....  */
	      
	      
	      /* Now go back to first triangle in list (if not at edge) */
	      
	      if (events[ii].src_no > AT_EDGE) 
		{
		  
		  atTriNode2 = events[ii].tris;
		  /* atTriNode1 is pointing to last triangle in list*/

		  /* Again, make sure that there are two vertices
		     in common between the two triangles, eg:
		            (a,b) -> (b,c) -> (c,a)
		     */

		  flag = 0;
		  for ( jj=0; jj<3; jj++ )
		    {
		      for ( kk=0; kk<3; kk++)
			{
			  if ( atTriNode1->at->verts[jj] == atTriNode2->at->verts[kk] )
			    {
			      flag++;

			    } /* endif */

			} /* end for (kk.. */

		    } /* end for (jj... */
		  
		  /* If not 2 verts in common, at edge! */

		  if ( flag < 2 )
		    {
		      atTriNode2        = NULL;
		      area              = 0.0;
		      events[ii].src_no = AT_EDGE;
		    }
		  else
		    {

		      /* Compute area and assign to event */

		      area+= (
			      ( ((double) atTriNode1->at->xc)*atTriNode2->at->yc) - 
			      ( ((double) atTriNode2->at->xc)*atTriNode1->at->yc)
			     );

		      events[ii].area = (double) fabs(area)/2.0;

		    } /* end if (flag < 2) */

		} /* end if (events[i].src_no < AT_EDGE) */

	    } /* end if (atTriNode2->next == NULL)  */
	  
	} /* end if (atTriNode == NULL) */
      

      
      if ( ( GlobalDebugLevel >= DEBUG_L3 )
	   && (events[ii].src_no > AT_EDGE)
	   && ( filen != NULL ))
	{
	  
	  fprintf( filen, "%s", "polygon(" );
	  
	  atTriNode1 = events[ii].tris;
	  while( atTriNode1 != NULL )
	    {
	      fprintf( filen,"%5.4E%s%5.4E", 
		       atTriNode1->at->xc,
		       ",",
		       atTriNode1->at->yc);
	      
              if ( atTriNode1->next != NULL)
		{
		  fprintf( filen, "%s", "," );
		}
	      
	      atTriNode1 = atTriNode1->next;
	    } /* end while( ... */

	  fprintf(filen,") # %lf\n", events[ii].area  );
	} /* end if ( (Global... */
      
      
      
    } /* end for (ii=0;ii<n;ii++) */
  
  
  if (( GlobalDebugLevel >= DEBUG_L3 ) && 
      ( filen != NULL ))
    {
      fclose( filen );
    }

}




/*
  +-------------------------------------------------
  +
  +  Calculate the center and radius of a cirlce circumscribed
  +  around three points by computing the intersection of the
  +  perpendicular bisectors.
  +
  +-------------------------------------------------
  */
void calculateCircle (
		      struct Event *e1,  /* i: Three events that make up  */
		      struct Event *e2,  /*    Delauny Triangle         */
		      struct Event *e3, 
		      double *xc,        /* o: coords of center of circumscribed */
		      double *yc,        /* o: circle around 3 points            */
		      double *r,         /* o: radius of circle                  */
		      double *xm,        /* o: x-center of mass of the triangle  */
		      double *ym         /* o: y-center of mass of the triangle  */
		     )
{
  double m12, m23;     /* slope from pt 1 to 2 and 2 to 3       */
  double x12, x23;     /* mid pt from pt 1 to 2 and 2 to 3 in X */
  double y12, y23;     /* mid pt from pt 1 to 2 and 2 to 3 in Y */
  double mp12, mp23;   /* perpincular slope from pt 1 to 2 and 2 to 3 */
  double dx,dy;        /* delta from center */
  
  double dx12, dx23, dy12, dy23;  /* delta in x and y to points */

    
  *xc = 0.0;
  *yc = 0.0;
  *r  = 0.0;

  /* compute center of mass */

  *xm = ( e1->coord.x + e2->coord.x + e3->coord.x )/3.0;
  *ym = ( e1->coord.y + e2->coord.y + e3->coord.y )/3.0;


  /* compute the delta X and delta Y */
  
  dx12 = ((double) e1->coord.x) - e2->coord.x;
  dx23 = ((double) e2->coord.x) - e3->coord.x;
  dy12 = ((double) e1->coord.y) - e2->coord.y;
  dy23 = ((double) e2->coord.y) - e3->coord.y;

  /* compute the half way points */

  x12 = (((double) e1->coord.x) + e2->coord.x)/2.0;
  x23 = (((double) e2->coord.x) + e3->coord.x)/2.0;
  y12 = (((double) e1->coord.y) + e2->coord.y)/2.0;
  y23 = (((double) e2->coord.y) + e3->coord.y)/2.0;


  /* if any of the delta's is zero => the angle from the two
     points is either 0 or 90 degress.  Handle this special case
     seperately.
     */
  if ( dx12 == 0.0 ) *yc = y12;
  if ( dx23 == 0.0 ) *yc = y23;
  if ( dy12 == 0.0 ) *xc = x12;
  if ( dy23 == 0.0 ) *xc = x23;


  /* if the points are not at right angles, then use the slopes to 
     find the intersection */

  if  ( ( dx12 != 0.0 ) && ( dx23 != 0.0 ) && (dy12 !=0.0) && (dy23 != 0.0))
    {
      /* slopes of the bisectors */

      m12 = dy12/dx12;
      m23 = dy23/dx23;

      /* slopes perpendicular to the bisectors */

      mp12 = -1./m12;
      mp23 = -1./m23;

      /* using slopes and location of mid-way point compute center */
      
      *xc = (double) (y12-y23-mp12*x12+mp23*x23)/(mp23-mp12);
      *yc = (double) mp12*(*xc-x12)+y12;
  
    }
  if (( dy12 == 0.0 ) || ( dy23 == 0.0 ))
    {
      /* Already found xc due to perpendicular or zero slope */

      if ( dy12 == 0.0 ) *yc =(double) -1.0*dx23/dy23*(*xc-x23)+y23;
      if ( dy23 == 0.0 ) *yc =(double) -1.0*dx12/dy12*(*xc-x12)+y12;
    }

  if (( dx12 == 0.0 ) || ( dx23 == 0.0 ))
    {
      /* Already found yc due to perpendicular or zero slope */

      if ( dx12 == 0.0 ) *xc = (double) -1.0*dy23/dx23*(*yc-y23)+x23;
      if ( dx23 == 0.0 ) *xc = (double) -1.0*dy12/dx12*(*yc-y12)+x12;
    }

  /* compute the radius */
  
  dx = (*xc - e2->coord.x);
  dy = (*yc - e2->coord.y);
  *r = (double) sqrt( ((double) dx*dx + dy*dy) );

}  






/*
  +----------------------------------------------------------
  +
  + Populate "cell" strucutres for Super Voronoi cells 
  +
  +----------------------------------------------------------
  */
void makeVoronoiCells(
		      struct Source   *firstSource,  /* i/o: first source in list         */
		      char            *superVfile,   /* i:  debug file for super V. cells */
		      int              debug         /* i: debug verbosity*/
		      )
{

  struct Source  *atSrc  = NULL;  /* current source   */
  struct EvtList *atEvt  = NULL;  /* current event    */
  struct Cell    *atCell = NULL;  /* current V. cell  */
  struct Points  *atPt   = NULL;  /* currnt Pixel     */
  struct TriList *atTri  = NULL;  /* current triangle */

  FILE  *fp = NULL; /* debug file */


  if ( debug >= DEBUG_L2 )
    {
      fprintf( GlobalDebugFile, "DEBUG: populating super Voronoi cells\n");
    }

  if ( debug >= DEBUG_L3 )
    {
      fp = fopen( superVfile, "w");
      if ( fp == NULL )
	{
	  err_msg( "WARNING: could not open debug file: %s\n",
		   superVfile );
	}
    }


  /* loop thru sources to populate structures */

  atSrc = firstSource;
  while ( atSrc != NULL )
    {
      if ( atSrc->evts != NULL )
      {

	atEvt = atSrc->evts;
	
	while ( atEvt != NULL )
	  {
	    
	    atCell = ( struct Cell *) calloc(1,sizeof( struct Cell ));
	    if ( atCell == NULL )
	      {
		err_msg( "ERROR: could not allocate sufficient memory\n");
		GlobalError |= BIT00;  /* CHECK THIS */
		return;
	      }
	    
	    /* put Cell into list */
	    
	    atCell->next = atSrc->cell;
	    atSrc->cell  = atCell;
	    atCell->evt  = atEvt->at;
	    
	    
	    if ( fp != NULL ) fprintf( fp, "polygon(");
	    
	    /* go thru triangles */

	    atTri = atEvt->at->superTris; 
	    while ( atTri != NULL )
	      {
		atPt = ( struct Points * )calloc( 1,sizeof( struct Points));
		if ( atPt == NULL )
		  {
		    err_msg( "ERROR: could not allocate sufficient memory\n");
		    GlobalError |= BIT00;  /* CHECK THIS */
		    return;
		  }
		
		atPt->next  = atCell->at; /* put into list */
		atCell->at  = atPt;
		
		atPt->coord.x = atTri->at->xc;
		atPt->coord.y = atTri->at->yc;
		
		if ( fp != NULL ) 
		  {
		    fprintf( fp, "%6.2f,%6.2f", atPt->coord.x, atPt->coord.y);
		    if ( atTri->next != NULL ) fprintf( fp, "," );
		  }

		/* check to see if at edge of field */
		if ( atTri->at->verts[0]->src_no < AT_EDGE ) atSrc->status |= BIT02;
		if ( atTri->at->verts[1]->src_no < AT_EDGE ) atSrc->status |= BIT02;
		if ( atTri->at->verts[2]->src_no < AT_EDGE ) atSrc->status |= BIT02;

		atTri = atTri->next;
	      } /* end while atTri != NULL */
	    
	    
	    if ( fp != NULL) fprintf( fp, ")\n");

	    atEvt = atEvt->next;
	  } /* end while atEvt != NULL */
	
      } /* end if atSrc->evts != NULL */
      
      atSrc = atSrc->next;
    } /* end while atSrc != NULL */

  if ( fp != NULL ) fclose(fp);

  
  if ( debug >= DEBUG_L2 )
    {
      fprintf( GlobalDebugFile, "DEBUG: done populating super Voronoi cells\n");
    }


  return;


}





/**************************************************************************
 *                                                                        *
 *  INPOLY.C
 *
 * Copyright (c) 1995-1996 Galacticomm, Inc.  Freeware source code.
 *
 * Please fell free to use this source code for any purpose, commerical
 * or otherwise, as long as you don't restrict anyone else's use of 
 * this source code.  Please give credit where credit is due.
 *
 * Point-in-polygon algorithm, created especially for World-Wide Web
 * servers to process image maps with mouse-clickable regions.
 *
 * Home for this file:  http: / /www.gcomm.com/develop/inpoly.c
 * [not active 
 *
 *      6/19/95 - Bob Stein & Craig Yap
 *      stei@gcomm.com
 *      craig@cse.fau.edu
 *
*************************************************************************/

int inPoly( struct Pt p,
	    struct Pt *points,
	    long       npoints
	    )
{

  double xnew, ynew;
  double xold, yold;
  double x1, y1;
  double x2, y2;

  long ii;
  int inside = 0;

  
  xold = points[npoints-1].x;
  yold = points[npoints-1].y;
  
  for ( ii=0; ii<npoints; ii++)
    {
      xnew = points[ii].x;
      ynew = points[ii].y;

      if ( xnew > xold )
	{
	  x1 = xold;
	  x2 = xnew;
	  y1 = yold;
	  y2 = ynew;
	}
      else
	{
	  x1 = xnew;
	  x2 = xold;
	  y1 = ynew;
	  y2 = yold;
	}

      if ( ( xnew < p.x ) == ( p.x <= xold ) &&
	   ( ( p.y - y1 ) * ( x2  - x1 ) <
	     ( y2  - y1 ) * ( p.x - x1 )    ) )
	{
	  inside = !inside;
	}

      xold = xnew;
      yold = ynew;
    } /* end for ii */

  return(inside);

}




/*
  +-------------------------------------------------------------
  +
  + Restore src_no after the super Cell percolation step.  This
  + is done by stepping thru triangles attached to a give enent
  + looking for ones that are ~-100
  +
  +-------------------------------------------------------------
  */
void resetId( 
	     struct Event *atEvt  /* i/o: current event to reset src_no of */
	     )
{
  int ii;    /* loop variable */
  struct TriList *atTri; /* current triangle */


  if ( atEvt->src_no < -100-(AT_EDGE) )
    {

      atEvt->src_no = -atEvt->src_no + 100; /* reset src_no */
      
      /* loop thru triangles */
      atTri = atEvt->tris;
      while ( atTri != NULL )
	{
	  for ( ii=0; ii<3; ii++)
	    {
	      resetId( atTri->at->verts[ii] );
	    } /* end for ii */

	  atTri = atTri->next;
	} /* end atTri */

    } /* end if */

  return;

}




long countBack( 
	       struct Event *atEvt, 
	       struct Pt    *sides, 
	       long         *nSides, 
	       long          nCell, 
	       FILE         *bkfile 
	       )
{

  long ii, jj;
  long count = 0;
  long offset = 0;
  struct TriList *atTri;


  if ( atEvt->src_no > AT_EDGE )
    {
      
      for ( jj=0;jj<nCell; jj++)  /* need to see if in ANY cell */
	{

	  if ( inPoly(atEvt->coord, sides+offset, nSides[jj] ) )
	    {
	      if ( ( atEvt->src_no == BACKGROUND ) ||
		   ( atEvt->src_no == TRIMMED_EVENT ) )
		{
		  count = atEvt->multi;
		  
		  if ( bkfile != NULL )
		    {
		      fprintf( bkfile, "circle(%f,%f,1.0)\n", 
			       atEvt->coord.x, atEvt->coord.y );
		    }

		} /* end if background */
	      
	      atEvt->src_no = -(100+atEvt->src_no);
	      
	      atTri = atEvt->tris;  /* looking for background */
	      
	      while ( atTri != NULL )
		{
		  for ( ii=0; ii<3; ii++ )
		    {
		      count += countBack( atTri->at->verts[ii], sides,
					  nSides, nCell, bkfile );
		    }
		  
		  atTri = atTri->next;
		}
	      /* break;  save from having to continue checks */
	    } /* end if inside */
	  
	  offset += nSides[jj];
	}
      
    } /* end if > AT_EDGE */
  
  return(count);
}







void countem( struct Source *firstSource, char *bkFileName)
{
  
  struct Source   *atSrc   = NULL;
  struct Cell     *atCell  = NULL;
  /* struct Event    *atEvt   = NULL; */
  /* struct TriList  *atTri   = NULL; */
  struct Points   *atPos   = NULL;
  struct Points   *nextPos = NULL;

  double tmpArea;

  struct Pt *sides;
  long      *nsides;
  long       ncells;
  long       npts, totpts;
  long       atcell = 0;

  char fname[100];

  FILE *bkfile = NULL;




  atSrc = firstSource;
  while ( atSrc != NULL )
    {
      if ( atSrc->evts != NULL)
	{

	  /* get no. cells in source */
	  ncells = 0;
	  atCell = atSrc->cell;
	  while ( atCell != NULL )
	    {
	      ncells += 1;
	      atCell = atCell->next;
	    }
	  nsides = ( long *)calloc( ncells, sizeof( long ));

	  /* get no. sides per cell and area */

	  totpts = 0;
	  atcell = 0;
	  atCell = atSrc->cell;
	  while ( atCell != NULL )
	    {
	      npts = 0;
	      /*  Compute total area of Super cells */
	      tmpArea = 0;
	      
	      atPos   = atCell->at;
	      nextPos = atPos->next;
	      
	      while ( nextPos != NULL )
		{
		  tmpArea += ((double )atPos->coord.x) * nextPos->coord.y -
		             ((double )nextPos->coord.x) * atPos->coord.y ;

		  atPos   = nextPos;
		  nextPos = nextPos->next;

		  npts += 1;
		  
		} /* end while nextPos != NULL */
	      nextPos = atCell->at;
	      tmpArea += ((double) atPos->coord.x) * nextPos->coord.y -
		         ((double) nextPos->coord.x) * atPos->coord.y ;
	      

	      npts += 1;
	      nsides[ atcell ] = npts;
	      totpts += npts;

	      tmpArea = fabs(tmpArea) / 2.0;

	      atCell->area = (double )tmpArea;
	      atCell->aexp = (double )tmpArea * atCell->evt->exp;

	      atSrc->superArea += atCell->area;
	      atSrc->superAexp += atCell->aexp;
	      
	      atCell  = atCell->next;
	      atcell += 1;
	    } /* while atCell != NULL */

	  sides = ( struct Pt *)calloc( totpts, sizeof( struct Pt ));

	  /* populate sides */
	  atcell = 0;

	  atCell = atSrc->cell;
	  while ( atCell != NULL )
	    {
	      atPos = atCell->at;
	      while ( atPos != NULL )
		{
		  sides[atcell].x = atPos->coord.x;
		  sides[atcell].y = atPos->coord.y;

		  atPos = atPos->next;
		  atcell += 1;
		}
	      atCell = atCell->next;
	    } /* end while atCell != NULL */

	  if ( GlobalDebugLevel >= DEBUG_L3 )
	    {
	      sprintf( fname, bkFileName, atSrc->Id );
	      bkfile = fopen( fname, "w");
	    }

	  atCell = atSrc->cell;

	  atCell->numBack = countBack( atCell->evt, sides, nsides, ncells, bkfile);
	  atSrc->superCount += atCell->numBack;
	  resetId( atCell->evt );
	  

	  if ( bkfile ) fclose( bkfile );
	  free( sides );
	  free( nsides);
	  
	  sides  = NULL;
	  nsides = NULL;

	} /* if atSrc->evts != NULL */

      atSrc = atSrc->next;

    } /* end while atSrc != NULL */

  
  return;

}




