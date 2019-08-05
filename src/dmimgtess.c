/*                                                                
**  Copyright (C) 2004-2008  Smithsonian Astrophysical Observatory 
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


#include <vtpdetect.h>
#include <cxcregion.h>
#include <dslib.h>
#include <dsnan.h>
#include <histlib.h>
#include <dmimgio.h>
#include <float.h>

/* Okay, these are all heisted from vtpdetect, they're needed
   to make things work by grabbing it's object files. */

int   GlobalError;
FILE *GlobalDebugFile;
int   GlobalDebugAuxFile;
int   GlobalDebugLevel;
bool  wcs_present;


long *G_lAxes;
dmDescriptor *xx_desc;
dmDescriptor *yy_desc;


void outputDebug( // Needed by vtp routines
                 int  atLevel,       /* i: output if debug level is >= */
                 char *msg,          /* i: message string       */
                 int  bORa           /* i: before or after call */
                 )
{}




long load_image( char *infile, struct Event **events, double minval,
                 dmBlock **inBlock);


long load_image( char *infile, struct Event **events, double minval,
                 dmBlock **inBlock)
{
  void *data;
  regRegion *dss;
  long null;
  short has_null;
  dmDataType dt;

  long xx,yy;

  long nvals;
  short *mask = NULL;
  
  if ( NULL == ( *inBlock = dmImageOpen( infile) ) ) {
    err_msg("ERROR: Cannot open image '%s'\n", infile );
    return(-1);
  }

  if ( dmUNKNOWNTYPE == ( dt = get_image_data( *inBlock, &data,  &G_lAxes, 
                                               &dss, &null, &has_null ) ) ) {
    err_msg("ERROR: Cannot get image data or unknown image data-type for "
            "file '%s'\n", infile);
    return(-1);
  }
  if ( 0 != get_image_wcs( *inBlock, &xx_desc, &yy_desc ) ) {
    err_msg("ERROR: Cannot load WCS for file '%s'\n", infile );
    return(-1);
  }

  mask = get_image_mask( *inBlock, data, dt, G_lAxes, dss, 
    null, has_null, xx_desc, yy_desc );
  


  nvals=0;

  /* A bit sloppy to alloc this much data, but easier than 
     reallocing */
  *events = (struct Event*) calloc( G_lAxes[0]*G_lAxes[1], sizeof(struct Event));
  
  /* would normally do the yy=lAxes[0];yy--; trick, but we want to get
     the data sorted in y and then x so need to do it this way or
     else need to sort later */
  for(yy=0;yy<G_lAxes[1];yy++) {

    for (xx=0;xx<G_lAxes[0]; xx++) {

      double dat;
      
      dat = get_image_value( data, dt, xx, yy, G_lAxes, mask );
      
      //dat = get_image_value( data, dt, G_lAxes, xx, yy,
      //                       dss, xx_desc, yy_desc, null, has_null );
      
      if ( ds_dNAN( dat ) || ( dat < minval)) {
        /* do nothing */
      } else {
        /* Note:  we could convert coordinates here from logical to 
         * physical, but no need to.  We have to convert the 
         * centers of the tesslation anyway later on, so no need to do
         * 'em both ... methinks 
         */
         
        struct Event *at;

        at=(*events)+nvals;

        at->coord.x = xx; 
        at->coord.y = yy;
        at->multi = 1;
        at->area = 0;
        nvals++;
      }

    } /* end for xx */
  } /* end for yy */


  return(nvals);
}





int convert_coords( dmDescriptor *xdesc,
                    dmDescriptor *ydesc,
                    double xx,  /*Note: Using double rather than long here! */
                    double yy, /* */
                    double *xat,
                    double *yat
                    );

int convert_coords( dmDescriptor *xdesc,
                    dmDescriptor *ydesc,
                    double xx,  /*Note: Using double rather than long here! */
                    double yy, /* */
                    double *xat,
                    double *yat
                    )
{
  if ( xdesc ) {
    double lgc[2];
    double phy[2];
    lgc[0] = xx+1;
    lgc[1] = yy+1;
    dmCoordCalc_d( xdesc, lgc, phy );
    if ( ydesc ) {
      dmCoordCalc_d( ydesc, lgc+1, phy+1 );
    }
    *xat = phy[0];
    *yat = phy[1];
  } else {
    *xat = xx;
    *yat = yy;
  }
  return(0);
}



int dmimgtess()
{


  char infile[DS_SZ_FNAME];
  char outfile[DS_SZ_FNAME];
  short edge = 0;
  char meth[100];
  float minval = 1;
  short clobber;

  enum { TESS, TRI, CIRCLE } style;

  struct Event *events=NULL;
  long         numEvents=0;

  long                  ntri;         /* number of triangles        */
  struct TempTri       *nodes=NULL;        /* list of triangle nodes     */
  struct Triangle      *firstTri=NULL;     /* pointer to first triangle */
  long ii;
  long removed=0;

  regRegion *outRegion = NULL;

  dmBlock *inBlock;
  dmBlock *outBlock;


  clgetstr( "infile", infile, DS_SZ_FNAME );
  clgetstr( "outfile", outfile, DS_SZ_FNAME );
  minval = clgeti( "minpix" );
  clgetstr( "method", meth, 99 );
  edge = clgeti( "edge" );
  clobber = clgetb("clobber");

  if ( strncmp( meth, "tess", 4 ) == 0 ) {
    style = TESS;
  }  else if ( strncmp( meth, "tri", 3 ) ==0 ) {
    style = TRI;
  }  else if ( strncmp( meth, "cir", 3 ) == 0 ) {
    style = CIRCLE;
  } else {
    err_msg("ERROR: Unknown method '%s'\n", meth );
    return(-1);
  }

  if ( ds_clobber( outfile, clobber, NULL ) != 0 ) {
    return(-1);
  }




  /* Have to have at least 3 events */
  if ( 3 > (numEvents = load_image( infile, &events, minval, &inBlock ))) {
    return(-1);
  }


  ntri  = -1;
  nodes = (struct TempTri *) calloc( numEvents*4, sizeof(struct TempTri) );
  doTriangulation( events, numEvents, BEFORE, nodes, &ntri, NULL);
  makeTriangle( events, numEvents, nodes+1, ntri, BEFORE, &firstTri); 
  free( nodes );
  for (ii=0;ii<numEvents;ii++) {
    events[ii].exp = 1;
  }
  calculateArea( events, numEvents, NULL );
  removeEdges( events, numEvents, edge, &removed ); 


  /* This is the new stuff */

  outRegion = regCreateEmptyRegion();

  for (ii=0;ii<numEvents;ii++) {
    struct TriList *at;
    long npoints;
    double *xx;
    double *yy;
    
    if ( events[ii].area == 0.0 ) 
      continue;

    
    if ( TESS == style ) {
      
      at = events[ii].tris;
      npoints =0;
      while (at) {
        npoints++;
        at = at->next;
      }
      
      xx = (double*)calloc( npoints+1, sizeof(double));
      yy = (double*)calloc( npoints+1, sizeof(double));
      
      at = events[ii].tris;
      npoints=0;
      while(at) {
        convert_coords( xx_desc, yy_desc, at->at->xc, at->at->yc, xx+npoints, yy+npoints );

        /* Check for 0 length polygon sides */
        if ((npoints>0) && ( fabs(xx[npoints]-xx[npoints-1])<FLT_EPSILON ) && (fabs(yy[npoints]-yy[npoints-1])<FLT_EPSILON) ) {
            // skip it, don't increment the counter
            at = at->next;
            continue;
        }

        npoints++;
        at = at->next;
      } // end while
      
      // Close polygon, last point equals the first
      xx[npoints]=xx[0]; 
      yy[npoints]=yy[0];


      if ( npoints >= 3 ) {
        long pre, post;
        pre = regGetNoShapes( outRegion );
        regAppendShape( outRegion, "Polygon", 1, 1, xx, yy, npoints+1, NULL, NULL,
                        0, 0 );
        post = regGetNoShapes( outRegion );
        if ( pre == post ) {
            err_msg("Problem creating region");
        }
        
      }
      
      free(xx); free(yy);


    } else if ( TRI == style ) {
      double tx[4], ty[4];

      
      at = events[ii].tris;
      while (at) {
        short kk;        
        for (kk=0;kk<3;kk++) {
          convert_coords( xx_desc,yy_desc, 
                at->at->verts[kk]->coord.x, at->at->verts[kk]->coord.y,
                tx+kk, ty+kk );
        } /* end for kk*/
        
        tx[3] = tx[0];
        ty[3] = ty[0];

        regAppendShape( outRegion, "Polygon", 1, 1, tx, ty, 4, NULL, NULL,
                        0, 0 );
        
        at = at->next;
      }/* end loop over triangles */



    } else if ( CIRCLE == style ) {

      double xx[3], yy[3];
      
      at = events[ii].tris;
      while (at) {
        short kk;        
        for (kk=0;kk<3;kk++) {
          convert_coords( xx_desc,yy_desc, 
                at->at->verts[kk]->coord.x, at->at->verts[kk]->coord.y,
                xx+kk, yy+kk );
        } /* end for kk*/
        

        double x1 = xx[0]; double x2 = xx[1]; double x3 = xx[2];
        double y1 = yy[0]; double y2 = yy[1]; double y3 = yy[2];
        double x, y;
        
        if ( x1 == x2 ) {
            y = (y1+y2)/2.0;
            double m = (y3-y2)/(x3-x2);
            x = x2 + (y2 - y)*m;
        } else if ( x3 == x2 ) {
            y = (y3+y2)/2.0;
            double m = (y1-y2)/(x1-x2);
            x = x2 + (y2 - y)*m;
        } else if ( y2 == y1 ) {
            x = (x1+x2)/2.0;
            double m = (y3-y2)/(x3-x2);
            y =  y2-(x-x2)/m;
        } else {
            double mr = (y2-y1) / (x2-x1);
            double mt = (y3-y2) / (x3-x2);

            if ( mr != mt ) {
              x = (mr*mt*(y3-y1) + mr*(x2+x3) - mt*(x1+x2)) / (2*(mr-mt));
              y = (y1+y2)/2 - (x - (x1+x2)/2) / mr;
            } else { // Should never get here, but better than uninit memory
              x=0;
              y=0;
            }
        } 

        double radius = fabs(x2-x) > 0 ? hypot( (y2-y), (x2-x) ) : fabs(y2-y);

        regAppendShape( outRegion, "Circle", 1, 1, &x, &y, 1, &radius, NULL,
                          0, 0 );
        
        at = at->next;
      } /* end loop of triangles */
      
    } else {
      /* ERROR */
    }

  } /* end for ii */

  

  outBlock = dmTableWriteRegion( dmDatasetCreate( outfile ),
                           "REGION", NULL, outRegion ) ;

  
  dmBlockCopy( inBlock, outBlock, "HEADER"); 
  put_param_hist_info( outBlock, "dmimgtess", NULL, 0 );

  dmTableClose(outBlock);
  dmImageClose( inBlock );


    
  return(0);


}
