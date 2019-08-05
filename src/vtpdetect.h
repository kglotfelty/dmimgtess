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
 *
 * Ref. No.         Date
   ----------       -----
     0.1            1997 Dec 30  Initial putback of working copy
     0.1.1          1998 Jan 05  Misc updates.
     0.1.2          1998 Jan 07  Updated output file format
     0.2            1998 Feb 11  Added image input support and updated output
     0.3            1998 Apr 07  Update outputs
     0.4            1998 Jun 10  Working on super cell update
     0.5            1998 Jun 23  Moved AOTS lib to seperate directory
     0.5.1          1998 Jul 1   Cleaned up move of AOTS library
 *
H***************************************************************** */


#define VTP_VERSION "0.5.1"

#include "vtpdetect_aots.h"


/* ------------------- IRAF / ASCDS SETUP  --------------- */

#define SZ_FNAME     512 
#define KEY_LENGTH    80


#ifndef DETECT_VTP
#define DETECT_VTP

#include <parameter.h>
#include <string.h>

#include <stdio.h>
#include <bool.h>
#include <stdlib.h>
#include <math.h>

#ifndef STACK_H
#include "stack.h"
#endif

#ifndef PIXLIB_H
#include "pixlib.h"
#endif 

#include "dslib.h"

#endif



/* ---------------  WARNINGS AND ERROR FLAGS  ------------ */


extern int   GlobalError;           /* Flag to hold error bits */

extern long    GlobalFPEcount;      /* Number of doubleing point errors
				       (FPE) avoided */

#define FPE_DETECT           -1.0   /* Number to return if a possible doubleing
				       point under/over flow occurs */

extern int     GlobalDebugAuxFile;  /* flag to hold error bits if debugging 
				       files could note be opened */



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

#define BEFORE                  0  /* comment is at start of function */
#define AFTER                   1  /* comment is at end of function   */

extern FILE *GlobalDebugFile;

/* ------------------  PROCESSING FLAGS  ----------------- */

#define FUDGE_FACTOR_1    1.0245
#define AT_EDGE               -4
#define TRIMMED_EVENT         -3
#define LOCK_BLENDED          -2
#define BLENDED_EVENT         -1
#define BACKGROUND             0



/* ------------------ OUTPUT DEFINITIONS ----------------- */

#define UNUSED         0
#define EXTNAME   "SRCLIST"


/* ------------------- DATA STRUCTURES  ------------------ */


extern struct ParameterList *param;
extern bool wcs_present;
extern struct DS_WCS_Data wcs_data;

struct ParameterList 
{
                     /* input */
  char  evtsfile[SZ_FNAME] ;            /* input file name                */
  char  evtext[KEY_LENGTH];             /* extension (block) name         */
  char  xcol[KEY_LENGTH];               /* x (RA?) column name            */
  char  ycol[KEY_LENGTH];               /* y (DEC?) column name           */
  char  expmap[SZ_FNAME];               /* ?? exposure map file name      */
  char  ftype[SZ_FNAME];                /* input file type, image or events */
  int   imgmin;                         /* min. no. events/image pixel to consider*/

                      /* processing */
  double fitCDFbelow;                    /* Fit background CDF below this flux */
  double fitParaMin;                     /* Fit residual (real - expected) flux between */
  double fitParaMax;                     /*      this Min and Max                       */

  double possionTol;                     /* Fit tolerance of the Poission fit           */
  double possionBin;

  double falseSrcLimit;                  /* prob. false source limit                    */
  int   coarseSrcLimit;                 /* min. no. events per source (coarse cutoff)  */
  double thresholdScale;                 /* threshold scaling factor                    */

  int   inFromEdge;                     /* ignore this many events in from the edge    */
  int   maxIterations;                  /* max. iteration to allow if it doesn't converge */

  double cutoffLo;                       /* min flux cutoff threshold to use            */
  double cutoffHi;                       /* max flux cutoff threshold to use            */

  double splitLevel;                     /* NOT USED CURRENTLY                          */

                      /* output */
  char  outfile[SZ_FNAME];              /* output file name                               */

  char region_filename[DS_SZ_PATHNAME]; /* filename for ASCII regions */


  double ellipse_sigma;	        /* size of ellipse radii in sigmas */

  bool  clobber;                        /* overwrite existing output file?       */

                      /* debugging */
  int   debugLevel;                     /* debugging level                            */
  char  debugFile[SZ_FNAME];            /* debug file name:  stdout and stderr handled */

  char  sourceFile[SZ_FNAME];           /* These file names are used to dump the evnets */
  char  backFile[SZ_FNAME];             /*   belonging to a source, the background, the */
  char  trimFile[SZ_FNAME];             /*   edge, or were trimmed away to SAOImage     */
  char  edgeFile[SZ_FNAME];             /*   region files.                              */
  char  blendFile[SZ_FNAME];            /*                                              */

  char  backCDF[SZ_FNAME];              /* Cummulative Distribution Function of background */
  char  totalCDF[SZ_FNAME];             /*   fluxes, all event fluxes, and residual of     */
  char  residualCDF[SZ_FNAME];          /*   real - expected flux (filenames)              */

  char  delfile[SZ_FNAME];              /* Filename for Delauny Triangle region file      */
  char  vorfile[SZ_FNAME];              /* Filename for Voronoi Cells region file         */

  char  pixfile[SZ_FNAME];              /* pixel locations for each source                */


  double fitStart;  
  char  superDelFile[SZ_FNAME];
  char  superVorFile[SZ_FNAME];
  char  superBkgFile[SZ_FNAME];

};




struct TriList    /* Triangle Node (string together to make list)*/
{
  struct Triangle *at;            /* Points to current triangle */
  struct TriList  *next;          /* Points to next node        */
};


struct Pt        /* Point in space  */
{
  double x,y;
};


struct Event    /* Event (Voronoi cell) attributes  */
{
  struct Pt       coord;      /* event coordinates           */
  long            multi;      /* no. events at that location */
  double           area;       /* area of Voronoi cell        */
  double           exp;        /* exposure time of this point */
  int             src_no;     /* source no. flag.            */
  struct TriList *tris;       /* pointer to list of triangles*/
                              /*   that have this event as a */
                              /*   vertices                  */
  struct TriList *superTris;  /* pointer to super cell D.T.'s */

};


struct Triangle   /* A Delauny triangle            */
{
  struct Event *verts[3];    /* pointers to 3 events that make the verticies       */
  double         xc,yc;       /* x and y location of center of circumscribed circle */
};



struct Cell   /* a Voronoi cell */
{
  struct Event  *evt;
  double          area;
  double          aexp;
  long           numBack;
  struct Points *at;
  struct Cell   *next;
};


struct Points
{
  struct Pt      coord;
  struct Points *next;
};




struct Angles     /* Array of angles the events make with center of triangle */
{
  double angs[3];
};

/*
struct TempTri    /K List event indicies that make up a Delauny triangle    K/
{
  long id1, id2, id3;
};
*/


struct Source    /* Properties of each source            */
{
  long            Id;         /* source number     */

  double           aexp;       /* total (flux) of events in source   */
  double           area;       /* just total area     */
  double           exp;        /* just total exposure */


  long            multi;      /* total no. events in source         */
  long            npix;       /* total no. pixels in source      */
  double           nback;      /* estimated no. background events */
  double           dnback;     /* est. error in background events */

  double           ncorr;      /* corrected counts */
  double           dncorr;     /* error in expt counts */

  double           rate;       /* est. count rate */
  double           drate;      /* error in rate   */

  long            blendMulti; /* total no. of blended events in src */

  double           sig;        /* signifigance of source (prob fake) */

  double           x;
  double           y;          /* X and Y position in pixels         */
  double           xpsig;      
  double           ypsig;      /* position uncertanties */
  double           sigMinor;
  double           sigMajor;   /* major and minor axes of ellipse in pixels */
  double           angle;      /* angle of ellipse in degrees        */

  double           xm2;
  double           ym2;

  double           ra;
  double           dec;        /* X and Y position in degrees        */
  double           rapsig;
  double           decpsig;    /* RA & DEC positional uncertanties */
  double           sigxs;
  double           sigys;      /* major and minor axes of ellipse in arcseconds */
 
  double           tflux;
  double           bflux;      /* (TEMP) total source flux and exp. back flux */

  double           cutoff;     /* cutoff limit (different for each src) */

  int             status;     /* coded status of source: edge, blended, etc */

  long            superCount;
  double           superArea;
  double           superAexp;


  struct Source  *next;       /* points to next source                */

  struct EvtList *evts;       /* points to list of events that define this src */
                              /*   if == NULL, then this is not a real src     */
  struct Cell    *cell;
};




struct EvtList  /* Event list node... string nodes together to make event list */
{
  struct Event   *at;      /* points to current event  */
  struct EvtList *next;    /* points to next node      */
};




struct PercParams  /* Parameters used by the percolation routines */
{
  struct Event *events;          /* pointer to events                       */
  long          numEvts;         /* number of events                        */
  double         cutoff;          /* flux cutoff threshold                   */
  double         exptBack;        /* expected no. background events in field */
  double         dexptBack;       /* est. error in exptBack                  */
  double         backRate;        /* background rate                         */
  double         dbackRate;       /* est. error in backRate                  */
  double         tam2bam;         /* field to background area ratio          */
  double         tarea;           /* total area (now FLUX) of field          */
  double         tJustarea;       /* total area (only)                       */
  double         texp;            /* total exposure                          */
  long          multi;           /* total no. events in field               */
  double         falseSrcLimit;   /* false source limit                      */
  int           coarseLimit;     /* min. no. events per source limit        */
  double         fmean;           /* average flux (normalization factor)     */
  int           blendFlag;       /* flag these events as blended?           */
  double         fluxNorm;        /* normalization for total flux array      */
  double         bfluxNorm;       /* normalization for background flux array */

};



typedef struct
{
  long ngrp;
  dmDescriptor** wcs_descriptor;
  dmBlock *parentBlock;
} WCS_str;


  

/* ------------------ FUNCTION PROTOTYPES ---------------- */

/* see source code for detailed comments on parameters */



/* In  _cleanup.c  */

extern void eraseSources(struct Source *srcs);

extern void removeEdges(struct Event *events, long n, int level,
			long *removed);

extern void eraseEvents( struct Event *evts, long n );


/*  in compare.c */

extern int sortYX(const void *i, const void *j);

extern int sortIdYX(const void *i, const void *j);

extern int sortLo2Hi(const void *i, const void *j);


/* in input.c */

extern void getEvents(struct Event **events, long *n, 
		      struct ParameterList *par,
		      WCS_str *wcs);

extern void getEventsLst( struct Event **events, long *n, dmBlock *evts,
			  char *evtext, char *xcol, char *ycol,
			  WCS_str *wcs);
			  
extern void getEventsImg( struct Event **events, long *n, dmBlock *evts,
			  int imgmin,
			  WCS_str *wcs);


extern void removeMulti( struct Event *events, long n, long *multi);

extern void getExposureMap( struct Event *events, long n, char *expmap,
			    char *evtsfile, char *evtext,
			    WCS_str *wcs);


extern void getDmVersion( char **version );

/* in output.c */

extern void calculateSourceParameters(struct Source *atSrc, double a, double da,
				      double tarea, long multi, 
				      double tam2bam,
				      char *pixfile,
				      struct Event *events
				      );

extern void outputTable( struct Source *atSrc, struct ParameterList *param,
			 struct PercParams *pprms, int iter, long nedge,
			 bool superdo, WCS_str *wcs);

extern void outputDebug( int level, char *msg, int bORa);


extern void clobberOutput( char *fname, bool clobber );


/* from aots_voronoi_main.c */


/* in percolation.c */



extern void doPerc( struct PercParams *pprms, int *nsrcs, 
		    struct Source **firstSrc);


extern void makeSources( struct PercParams *pprms, long i,
			 struct Source **firstSrc, struct Source **atSrc,
			 int *nsrcs, int level);

extern void checkSourceSig ( struct PercParams *pprms, struct Source **firstSrc,
			     struct Source **atSrc, int *nsrcs, int level);


extern void recursivePercolate( struct PercParams *pprms, long at,
				int src_no, struct Source *atSrc, int level);

extern void splitBlends(  struct PercParams *pprms, int *nsrcs, 
			  struct Source **firstSrc, double split );


/* in simul.c */

extern double calcSimCDF( double a, double y);

extern double calcSimCDFdn( double a, double y);

extern double gammln( double x);

extern double fluct_max( double a, double fmin,
			double tam2bam, double fluct_n);

extern double probFake(double a, double fmin, double tam2bam, double xncorr);



/* in threshold.c */

extern void calcCutoff( struct Event *events, long n, long nold,
			struct ParameterList *param, long iter, 
			double *cutoff, double *a, double *dextpback, 
			double *backRate, double *dbackRate, double *sf,
			double *fluxn, double *bfluxn);

extern double dmaxlike(long ntot, long nfit, double *flux, double a);

extern void  fitPoissian(double *bflux, long nbpts, double below,
			 double toler, double start, double *a);

extern void calcParaFit( double *xfit, double *yfit, long nfit, double *fminf);

extern void getAscfitVersion ( char **version );

extern void newFitPoissian( double *bflux, long nbpts, double fitBelow, 
			    double binsize, double *exptBack,
			    double *dexptBack , double *backRate, double *dbackRate);


/* in voroni.c */


extern void makeTriangle( struct Event *events, long n, struct TempTri *nodes,
			  long ntri, short bORa,struct Triangle **allTriangle);

extern void calculateArea( struct Event *events, long n, 
			   char *vorfile);

extern void calculateCircle( struct Event *e1, struct Event *e2, struct Event *e3,
			     double *xc, double *yc, double *r, double *xm, double *ym);



/* in Main */

extern void getParameters( struct ParameterList *);

extern int vtpdetect(void);

extern void debugEvents( struct ParameterList *par, struct Event *evt, 
			 long n , long j);


extern void getVersion( char **ver);


extern void debugParameters( struct ParameterList *par );

extern void makeVoronoiCells( struct Source   *firstSource, char *file,
			      int debug);
extern void countem( struct Source *firstSource, char *bkFileName);

/* in interface.c */

extern void doTriangulation( struct Event *events, long numEvts,
			     int excludeBack, struct TempTri *nodes,
			     long *ntri, char *delfile );

