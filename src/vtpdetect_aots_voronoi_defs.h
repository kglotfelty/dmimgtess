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


#ifndef NULL
#define NULL 0
#endif
#define DELETED -2


#include "vtpdetect_aots.h"

extern int triangulate, sorted, plot, debug;

struct	Freenode	{
struct	Freenode	*nextfree;
};
struct	Freelist	{
struct	Freenode	*head;
long			nodesize;
};

extern char *getfree(struct	Freelist *fl);
/* char *malloc(void); */
extern char *myalloc(unsigned foo);



extern float xmin, xmax, ymin, ymax, deltax, deltay;




/* Change the names to protect the innocent 

#define Site      event
#define sitenbr   Id

*/


extern struct	Site	*sites;
extern long		nsites;
extern long		siteidx;
extern long		sqrt_nsites;

extern struct 	Freelist sfl;
extern struct	Site	*bottomsite;



struct Edge	{
float		a,b,c;
struct	Site 	*ep[2];
struct	Site	*reg[2];
long		edgenbr;
};
#define le 0
#define re 1


/*typedef struct POINT Point; */


extern float dist( struct Site *s,struct Site *t);


extern struct Point PQ_min(void);
extern struct Halfedge *PQextractmin(void);
extern struct Edge *bisect( struct Site *s1, struct Site *s2);



struct Halfedge {
struct Halfedge	*ELleft, *ELright;
struct Edge	*ELedge;
long		ELrefcnt;
char		ELpm;
struct	Site	*vertex;
float		ystar;
struct	Halfedge *PQnext;
};

extern int has_endpoint(void),
  right_of(     struct Halfedge *el, struct Point *p);

extern struct	Halfedge *ELleftend, *ELrightend;
extern struct	Halfedge *HEcreate( struct Edge *e, int pm), 
  *ELleft(struct Halfedge *he), 
  *ELright(struct Halfedge *he), 
  *ELleftbnd(struct Point *p);
extern struct	Site *leftreg(struct Halfedge *he), 
  *rightreg(struct Halfedge *he);

extern struct Site *intersect(
     struct Halfedge *el1, 
     struct Halfedge *el2);
/*     struct Point *p); */


extern struct	Halfedge *PQfind(void);
extern int PQempty(void);


extern long n_allocs;
extern char **allAlloc;
extern long nallbuf;

extern void PQinitialize(void);
extern void ELinitialize(void);
extern void out_site(struct Site *s);
extern void out_ep(struct Edge *e);
extern void out_triple(struct Site *s1, struct Site *s2,  struct Site *s3,
     struct TempTri *verts, long *ntri, FILE *filen);

extern void ELdelete( struct Halfedge *he);
extern void PQdelete( struct Halfedge *he);

extern void PQinsert(struct Halfedge *he, struct Site *v, float offset);
extern void ELinsert(struct Halfedge *lb, struct Halfedge *new);

extern int makevertex(struct Site *v);
extern void endpoint(struct Edge *e, int lr, struct Site *s);

extern void freeinit( struct Freelist *fl, long size);


extern void ref(struct Site *v);
extern void deref(struct Site *v);

extern void voronoi(int triangulate, struct Site *(*)(void), struct TempTri *vts,
		    long *ntri, char *delfile);

extern void out_bisector(struct Edge *e);
extern void out_vertex(struct Site *v);

extern void geominit(void);
extern void plotinit(void);

extern void makefree(struct Freenode *, struct Freelist *);
