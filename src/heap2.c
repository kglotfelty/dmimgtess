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

extern long PQcount;
extern long PQhashsize;
extern struct	Halfedge *PQhash;
extern  long PQmin;


extern long PQbucket(struct Halfedge *he);


void PQinsert(struct Halfedge *he, struct Site *v, float offset)
{
  struct Halfedge *last, *next;
  
  he -> vertex = v;
  ref(v);
  he -> ystar = v -> coord.y + offset;
  last = &PQhash[PQbucket(he)];
  while ((next = last -> PQnext) != (struct Halfedge *) NULL &&
	 (he -> ystar  > next -> ystar  ||
	  (he -> ystar == next -> ystar && v -> coord.x > next->vertex->coord.x)))
    {	last = next;};
  he -> PQnext = last -> PQnext; 
  last -> PQnext = he;
  PQcount += 1;
  return;
}
