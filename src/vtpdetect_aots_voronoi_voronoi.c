/*
 * The author of this software is Steven Fortune.  Copyright (c) 1994 by AT&T
 * Bell Laboratories.
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

#
#include "vtpdetect_aots_voronoi_defs.h"
#include <stdio.h>
#include <float.h>
#include <limits.h>

struct	Site	*bottomsite;

extern FILE *GlobalDebugFile;

/* implicit parameters: nsites, sqrt_nsites, xmin, xmax, ymin, ymax,
   deltax, deltay (can all be estimates).
   Performance suffers if they are wrong; better to make nsites,
   deltax, and deltay too big than too small.  (?) */

void voronoi(
	     int triangulate,
	     struct Site *(*nextsite)(void),
	     struct TempTri *verts,
	     long *ntri,
	     char *delfile
	     )
{
  struct Site *newsite, *bot, *top, *temp, *p;
  struct Site *v;
  struct Point newintstar;
  int pm;
  struct Halfedge *lbnd, *rbnd, *llbnd, *rrbnd, *bisector;
  struct Edge *e;
  
  newintstar.x = -FLT_MAX;
  newintstar.y = -FLT_MAX;
  
  
  FILE *filen = NULL;
  
  if ( GlobalDebugLevel >= DEBUG_L3 )
    {

      if ( ( strcmp( delfile, "none") != 0 ) &&
	   ( strcmp( delfile, ""    ) != 0 )
	   )
	{
	  filen = fopen( delfile , "w");
	  if ( filen == NULL)
	    {
	      fprintf( GlobalDebugFile, 
		       "WARNING: could not open Delaunay triangles region file: %s\n", delfile);
	      GlobalDebugAuxFile = ( GlobalDebugAuxFile | BIT01 );
	    }
	} /* end if file */
    } /* end if debug*/
  
  
  PQinitialize();
  bottomsite = (*nextsite)();
  out_site(bottomsite);
  ELinitialize();
  
  newsite = (*nextsite)();
  while(1)
    {
      if(!PQempty()) newintstar = PQ_min();
      
      if (newsite != (struct Site *)NULL 
	  && (PQempty() 
	      || newsite -> coord.y < newintstar.y
	      || (newsite->coord.y == newintstar.y 
		  && newsite->coord.x < newintstar.x)))
	{/* new site is smallest */
	  out_site(newsite);
	  lbnd = ELleftbnd(&(newsite->coord));
	  rbnd = ELright(lbnd);
	  bot = rightreg(lbnd);
	  e = bisect(bot, newsite);
	  bisector = HEcreate(e, le);
	  ELinsert(lbnd, bisector);
	  if ((p = intersect(lbnd, bisector)) != (struct Site *) NULL) 
	    {	PQdelete(lbnd);
	    PQinsert(lbnd, p, dist(p,newsite));
	    };
	  lbnd = bisector;
	  bisector = HEcreate(e, re);
	  ELinsert(lbnd, bisector);
	  if ((p = intersect(bisector, rbnd)) != (struct Site *) NULL)
	    {	PQinsert(bisector, p, dist(p,newsite));	
	    };
	  newsite = (*nextsite)();	
	}
      else if (!PQempty()) 
	/* intersection is smallest */
	{	
	  lbnd = PQextractmin();
	  llbnd = ELleft(lbnd);
	  rbnd = ELright(lbnd);
	  rrbnd = ELright(rbnd);
	  bot = leftreg(lbnd);
	  top = rightreg(rbnd);
	  out_triple(bot, top, rightreg(lbnd), verts, ntri, filen);
	  v = lbnd->vertex;
	  makevertex(v);
	  endpoint(lbnd->ELedge,lbnd->ELpm,v);
	  endpoint(rbnd->ELedge,rbnd->ELpm,v);
	  ELdelete(lbnd); 
	  PQdelete(rbnd);
	  ELdelete(rbnd); 
	  pm = le;
	  if (bot->coord.y > top->coord.y)
	    {	temp = bot; bot = top; top = temp; pm = re;}
	  e = bisect(bot, top);
	  bisector = HEcreate(e, pm);
	  ELinsert(llbnd, bisector);
	  endpoint(e, re-pm, v);
	  deref(v);
	  if((p = intersect(llbnd, bisector)) != (struct Site *) NULL)
	    {	PQdelete(llbnd);
	    PQinsert(llbnd, p, dist(p,bot));
	    };
	  if ((p = intersect(bisector, rrbnd)) != (struct Site *) NULL)
	    {	PQinsert(bisector, p, dist(p,bot));
	    };
	}
      else break;
    };
  
  for(lbnd=ELright(ELleftend); lbnd != ELrightend; lbnd=ELright(lbnd))
    {	e = lbnd -> ELedge;
    out_ep(e);
    };
  
  if (( GlobalDebugLevel >= DEBUG_L3 ) && ( filen != NULL ))
    {
      fclose(filen);
    }
  
  
  return;


}
