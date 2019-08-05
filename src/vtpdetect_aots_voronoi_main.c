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
#include <stdio.h>
#include "vtpdetect_aots_voronoi_defs.h"
struct Site *readone(void);
struct Site *nextone(void);

int scomp(struct Point *s1, struct Point *s2);

void freeHash( struct Halfedge *he);
void freeNode ( struct Freelist *node );



void readsites(void);
char **allAlloc;
long		nsites;
long		siteidx;
struct	Site	*sites;
void voroniInterface(struct Site *inSites,
		     long  n, 
		     struct TempTri *verts, 
		     long *ntri, 
		     char *delfile)
{	
  long i;
  struct Site *(*next)(void);
  
  *ntri = 0; 
  
  sites = inSites;
  nsites = n;
  
  sorted = 1; triangulate = 1; plot = 0; debug = 0;
  
  n_allocs = -1;
  nallbuf=0;
  allAlloc=(char **)calloc(4000,sizeof(char *));
  
  
  freeinit(&sfl, sizeof *sites);

  readsites();
  next = nextone;
  
  siteidx = 0;
  geominit();
  if(plot) plotinit();
  
  voronoi(triangulate, next, verts, ntri, delfile); 
  
  /*   free(sites); */
  
  
  for (i=0;i<n_allocs;i++)
    {
      free(allAlloc[i] );
    }
  
  
  
  free(allAlloc);
  return;
  
}





/* sort sites on y, then x, coord */
int scomp(s1,s2)
     struct Point *s1,*s2;
{
  if(s1 -> y < s2 -> y) return(-1);
  if(s1 -> y > s2 -> y) return(1);
  if(s1 -> x < s2 -> x) return(-1);
  if(s1 -> x > s2 -> x) return(1);
  return(0);
}



/* return a single in-storage site */
struct Site *nextone(void)
{
  struct Site *s;
  if(siteidx < nsites)
    {	
      s = &sites[siteidx];
      siteidx += 1;
      return(s);
    }
  else	return( (struct Site *)NULL);
}




/* read all sites, sort, and compute xmin, xmax, ymin, ymax */
void readsites(void)
{
  long i;
  /*  Remove all this stuff...just find min and max 
      
      
      nsites=0;
      sites = (struct Site *) myalloc(4000*sizeof *sites);
      while(scanf("%f %f", &sites[nsites].coord.x, &sites[nsites].coord.y)!=EOF)
      {	
      sites[nsites].sitenbr = nsites;
      sites[nsites].refcnt = 0;
      nsites += 1;
      if (nsites % 4000 == 0) 
      sites = (struct Site *) realloc(sites,(nsites+4000)*sizeof*sites);
      };
      qsort(sites, nsites, sizeof *sites, scomp);
    */

  xmin=sites[0].coord.x; 
  xmax=sites[0].coord.x;
  for(i=1; i<nsites; i+=1)
    {	
      if(sites[i].coord.x < xmin) xmin = sites[i].coord.x;
      if(sites[i].coord.x > xmax) xmax = sites[i].coord.x;
      sites[i].refcnt=0;
    };
  ymin = sites[0].coord.y;
  ymax = sites[nsites-1].coord.y;

 return;
}



/* read one site */
struct Site *readone(void)
{
struct Site *s;

s = (struct Site *) getfree(&sfl);
s -> refcnt = 0;
s -> sitenbr = siteidx;
siteidx += 1;
if(scanf("%f %f", &(s->coord.x), &(s->coord.y)) == EOF)
	return ((struct Site *) NULL );
return(s);
}






void freeHash( struct Halfedge *he)
{

  struct Halfedge *nextL, *nextR ;
  
  if ( he != NULL )
    {

      free(he->ELedge);
      he->ELedge = NULL;


      nextL=he->ELleft;
      nextR=he->ELright;

      free(he->ELleft);
      free(he->ELright);

      he->ELleft=NULL;
      he->ELright=NULL;

      freeHash( nextL );
      freeHash( nextR );

      free(he);

    }

  return;
}




void freeNode ( struct Freelist *node )
{
  struct Freenode *at, *next;
  
  at = node->head;

  do
    {
      next=at->nextfree;
      free(at);
      at->nextfree = NULL;
      at=next;      
    } while ( at != NULL );

  return;

}

