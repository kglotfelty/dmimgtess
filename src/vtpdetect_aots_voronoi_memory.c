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
long n_allocs;
long nallbuf;
extern FILE *GlobalDebugFile;

void makefree(curr,fl)
     struct Freenode *curr;
     struct Freelist *fl;
{
  curr -> nextfree = fl -> head;
  fl -> head = curr;
  return;
}


void freeinit(fl, size)
     struct	Freelist *fl;
     long	size;
{
  fl -> head = (struct Freenode *) NULL;
  fl -> nodesize = size;
  return;
}



char *getfree(fl)
     struct	Freelist *fl;
{
  long i; struct Freenode *t;
  if(fl->head == (struct Freenode *) NULL)
    {	
      t =  (struct Freenode *) myalloc(sqrt_nsites * fl->nodesize);
      for(i=0; i<sqrt_nsites; i+=1) 	
	makefree((struct Freenode *)((char *)t+i*fl->nodesize), fl);
    };
  t = fl -> head;
  fl -> head = (fl -> head) -> nextfree;
  return((char *)t);
}



long total_alloc;
char *myalloc(n)
     unsigned n;
{
  
  void *t;
  
  t=calloc(1,n);
  
  if ( t == NULL )
    {    fprintf(GlobalDebugFile,"Insufficient memory processing site %ld (%ld bytes in use)\n",
		 siteidx, total_alloc);
    exit(1);
    };
  total_alloc += n;
  n_allocs++;
  
  
  
  if (( n_allocs % 4000)  == 0 )
    {
      if ( nallbuf > 0 ) 
	{
	  nallbuf++;
	  allAlloc=(char **)realloc(allAlloc, nallbuf*sizeof(char *)*4000);
	}
    }


  if ( n_allocs == 0 ) nallbuf = 1;
  
  
  (allAlloc[n_allocs])=(char *)t;
  
  
  return(t);
}





