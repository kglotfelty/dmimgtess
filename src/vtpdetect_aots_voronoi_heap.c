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


#include "vtpdetect_aots_voronoi_defs.h"

long PQcount;
long PQhashsize;
struct	Halfedge *PQhash;
long PQmin;


long PQbucket(struct Halfedge *he);


void PQdelete( struct Halfedge *he)
{
  struct Halfedge *last;
  
  if(he ->  vertex != (struct Site *) NULL)
    {	
      last = &PQhash[PQbucket(he)];
      while (last -> PQnext != he) last = last -> PQnext;
      last -> PQnext = he -> PQnext;
      PQcount -= 1;
      deref(he -> vertex);
      he -> vertex = (struct Site *) NULL;
      he -> PQnext = NULL;
    };
  return;
}



long PQbucket( struct Halfedge *he)
{
  long bucket;
  
  bucket = (he->ystar - ymin)/deltay * PQhashsize;
  if (bucket<0) bucket = 0;
  if (bucket>=PQhashsize) bucket = PQhashsize-1 ;
  if (bucket < PQmin) PQmin = bucket;
  return(bucket);
}



int PQempty()
{
  int retval = ((PQcount==0)?1:0);
  return(retval);
}


struct Point PQ_min()
{
  struct Point answer;
  
  while(PQhash[PQmin].PQnext == (struct Halfedge *)NULL) {PQmin += 1;};
  answer.x = PQhash[PQmin].PQnext -> vertex -> coord.x;
  answer.y = PQhash[PQmin].PQnext -> ystar;
  return (answer);
}


struct Halfedge *PQextractmin()
{
  struct Halfedge *curr;
  
  curr = PQhash[PQmin].PQnext;
  PQhash[PQmin].PQnext = curr -> PQnext;
  PQcount -= 1;
  return(curr);
}



void PQinitialize()
{
  long i; 
  /*struct Point *s; */
  
  PQcount = 0;
  PQmin = 0;
  PQhashsize = 4 * sqrt_nsites;
  PQhash = (struct Halfedge *) myalloc(PQhashsize * sizeof *PQhash);
  for(i=0; i<PQhashsize; i+=1) PQhash[i].PQnext = (struct Halfedge *)NULL;
  return;
}

