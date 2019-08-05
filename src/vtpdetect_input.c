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
 * FILE NAME:  vtpdetect_input.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * This file contains routines used by vtpdetect to read in
 * events and get the exposure map info.
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
#include "vtpdetect.h"
#include "ascdm.h"
#include <ds_wcs.h>

#include <float.h>

struct DS_WCS_Data wcs_data;

/* 
   +------------------------------------
   +
   + Routine to get evnets.  It has been generalized to
   + get events from either an event list (DEFAULT) or
   + an image
   +
   +------------------------------------
   */
void getEvents (
		struct Event  **events,      /* o: event list                */
		long   *numEvts,             /* o: number of events returned */
		struct ParameterList *param ,/* i: parameter list            */
		WCS_str  *wcs                /* o: wcs data */
		)
{

  dmBlock *tryBlock = dmBlockOpen(NULL,  param->evtsfile );


  if ( NULL == tryBlock ) {
    err_msg( "ERROR: could not open image: %s\n", param->evtsfile);
    GlobalError = ( GlobalError | BIT02 );
    return;
  }
  
  if ( dmIMAGE == dmBlockGetType( tryBlock ) )
    {
      strcpy( param->ftype, "image" );
      getEventsImg( events, numEvts,
		    tryBlock, 1, wcs);
    }
  else
    {
      strcpy( param->ftype, "events" );
      getEventsLst( events, numEvts,
		    tryBlock, param->evtext,
		    param->xcol, param->ycol, wcs);
    }

  return;
}





/* 
   +-------------------------------------------
   +
   + This routine used the data model to read in the
   + X and Y columns of the event list.
   +
   +-------------------------------------------
   */
void getEventsLst ( 
		struct Event  **events,   /* o: array of events      */
		long           *numEvts,  /* o: number of events     */
		dmBlock *evtTable,
		char           *evtext,   /* i: extension name       */
		char           *xcol,     /* i: column name for X    */
		char           *ycol,     /* i: column name for Y    */
		WCS_str  *wcs                /* o: wcs data */
		)
{

  dmDescriptor* xcolName = NULL;       /* X (RA) column                 */
  dmDescriptor* ycolName = NULL;       /* Y (DEC) column                */
  dmBlockType   tabChk = dmTABLE;


  double pos[2];
  int nCols;

  long          ii = 0;                /* generic counter */
  int           flag=0;

  outputDebug( DEBUG_L3, 
	       "getting events.\n", BEFORE ); /* output debugging msg */


  /* need to try to use native DM conventions.  Tables should be 
     entered as dataset[table][filter] .
     
     */

  flag=1;

  /* evtTable = dmDatasetTableOpen( evtsfile );*/



  if (evtTable != NULL) tabChk = dmBlockGetType(evtTable);

  nCols = dmTableGetNoCols( evtTable );


  ycolName = NULL;

  if ( nCols != 1 ) 
    {
      

      if ( nCols == 2 )
	{

	  xcolName = dmTableOpenColumnNo( evtTable, 1 );

	  ycolName = dmTableOpenColumnNo( evtTable, 2 );


	  dmDescriptorGetRange_d( xcolName, &(wcs_data.range_min[0]), 
				  &(wcs_data.range_max[0]));
	  dmDescriptorGetRange_d( ycolName, &(wcs_data.range_min[1]), 
				  &(wcs_data.range_max[1]));


	  wcs->wcs_descriptor = ( dmDescriptor **) calloc( 2,sizeof( dmDescriptor *)) ;      
	  wcs->wcs_descriptor[0] = dmDescriptorGetCoord( xcolName );
	  if ( wcs->wcs_descriptor[0] == NULL )
	    {
	      /* WARNING */
	    }
 	  
	  wcs->wcs_descriptor[1] = dmDescriptorGetCoord( ycolName );
	  if ( wcs->wcs_descriptor[1] == NULL )
	    {
	      /* WARNING */
	    }
 
	  wcs->ngrp = 2;

	} /* end nCols == 2 */
      else
	{ 

	  /* try to open column "pos" */
	  
	  xcolName = dmTableOpenColumn( evtTable, "pos" );
	  if ( xcolName == NULL )
	    {
	      /* ERROR */
	      
	      xcolName = dmTableOpenColumn( evtTable, "sky" );
	      if ( xcolName == NULL )
		{
		  /* Error */
		}
	    }
	  
	  if ( dmGetElementDim( xcolName ) != 2 )
	    {
	      /* ERROR */
	    }
      

	  /* get WCS */

	  wcs->wcs_descriptor = ( dmDescriptor **) calloc( 1,sizeof( dmDescriptor *)) ;
	  wcs->wcs_descriptor[0] = dmDescriptorGetCoord( xcolName );
	  if ( wcs->wcs_descriptor[0] == NULL )
	    {
	      /* WARNING */
	    }
	  wcs->ngrp = 1;

	  
	  dmDescriptorGetRange_d( dmGetCpt(xcolName, 1),
				  &(wcs_data.range_min[0]), 
				  &(wcs_data.range_max[0]));
	  dmDescriptorGetRange_d( dmGetCpt(xcolName, 2),
				  &(wcs_data.range_min[1]), 
				  &(wcs_data.range_max[1]));



	  
	}

    }
  else /* Open columns by number */
    {

      xcolName = dmTableOpenColumnNo( evtTable , 1 );
      if ( xcolName == NULL)
	{
	  /* ERROR */
	}


      if ( dmGetElementDim( xcolName ) != 2 )
	{
	  /* ERROR */
	}
      
      wcs->wcs_descriptor = ( dmDescriptor **) calloc( 1,sizeof( dmDescriptor *)) ;      
      wcs->wcs_descriptor[0] = dmDescriptorGetCoord( xcolName );
      if ( wcs->wcs_descriptor[0] == NULL )
	{
	  /* WARNING */
	}
      
      wcs->ngrp = 1;
      
      dmDescriptorGetRange_d( xcolName, &(wcs_data.range_min[0]), 
			      &(wcs_data.range_max[0]));
      dmDescriptorGetRange_d( ycolName, &(wcs_data.range_min[1]), 
			      &(wcs_data.range_max[1]));
      

      
    }

  *numEvts = dmTableGetNoRows( evtTable );     /* get number of events */
  
  
  *events = (struct Event *) calloc( *numEvts, sizeof(struct Event) );/* allocate memory */
  if ( *events == NULL )
    {
      err_msg("ERROR: could not allocate memory for %ld events.\n", *numEvts);
      GlobalError = ( GlobalError | BIT01 );
      return;
    }
  
  /* Get events */
  ii=0;
  
  do
    {

      if ( ycolName == NULL )
	{
	  dmGetVector_d( xcolName, pos, 2 );
	}
      else
	{
	  pos[0] = dmGetScalar_d(xcolName);
	  pos[1] = dmGetScalar_d(ycolName);
	}
      
      (*events+(ii))->coord.x = ((long)(10.0*pos[0]))/10.0;
      (*events+(ii))->coord.y = ((long)(10.0*pos[1]))/10.0;


      (*events+(ii))->multi   = 1.0;

      ii += 1;

    } while ( dmTableNextRow(evtTable ) != dmNOMOREROWS );
  

  *numEvts = ii;

  if (*numEvts < 3)
    {
      err_msg( "ERROR: not enough events to proceede: %ld\n",
	       *numEvts );
      GlobalError |= BIT00;  /* check */
      return;
    }


  wcs->parentBlock = evtTable ;


  outputDebug( DEBUG_L3, 
	       "getting events.\n", AFTER ); /* output debugging msg */

  return;

}





/*
  +--------------------------------------
  +
  + This routine is used to collapse several events that
  + happen at the same (x,y) location into a single event
  + with a multiplicity > 1.
  +
  + The input is assumed to be a sorted event list
  + (in this case by Y then X).  This is so that
  + events with the same X and Y will be arranged
  + sequentially.  
  +
  + The next event in the list is checked against the
  + current event.  If they are the same, the current
  + event's muliplicity is ++, and the next event's
  + source number is set to 1 (this will be used as
  + a flag to prune the event list later).  The next
  + event is incremented (and the procedure repeated)
  + until the (x,y) locations differ.
  +
  + This is done so that the procedure works in linear time
  + ie, O(n).
  +
  + ------------------------------------------
  */
void removeMulti(
		 struct Event *events,    /* i: event list */
		 long          numEvts,   /* i: number of events */
		 long         *multi      /* o: number of events removed
					        due to the multiplictiy
					        of n-(*multi) events
					        being increased */
		 )
{

  long ii,kk;   /* generic loop varialbes */



  outputDebug( DEBUG_L3, 
	       "collapsing multiple (x,y) events.\n", 
	       BEFORE);  /* output debugging info */
  

  *multi = 0;

  /* for all events */
  for ( ii=0; ii<numEvts-1; ii++ )
    {

      kk = ii + 1;      /* ii = current event,  kk = "next" event  */

      while ( ( fabs(events[ii].coord.x-events[kk].coord.x)/events[ii].coord.x < 2*FLT_EPSILON ) &&
	      ( fabs(events[ii].coord.y-events[kk].coord.y)/events[ii].coord.y < 2*FLT_EPSILON) &&
	      ( kk < numEvts ) )
	{

	  /* repeat while the x & y coords are the same and
	     we're not at the end of the event list */
	  
	  events[ii].multi += 1;    /* increase multiplicty of the first
				       event with this X and Y */
	  events[kk].src_no = 1;    /* set the source number of the
				       duplicate event to 1 (flag to be
				       removed later */
	  kk               += 1;    /* compare against next event in list */
	  (*multi)         += 1;    /* count number of duplicated removed */

	} /* end loop over while events (x,y) are same, while ((... */


      ii = kk - 1;  /* go back one...that's where the difference occured */
      
    } /* end loop over all events,  for (ii=1 ;i<n; .... */


  outputDebug( DEBUG_L3,
	       "collapseing multiple (x,y) events.\n", 
	       AFTER );  /* output debugging msg */

  return;

}







/*
  +------------------------------------------------
  +
  + Get the exposure time at every unique X,Y position
  +
  +
  + Temporary hack.  Use ONTIME (if available).
  +
  +
  +------------------------------------------------
  */
void getExposureMap(
		    struct Event *events,   /* i/o: array of events    */
		    long          numEvts,  /* i: number of events     */
		    char         *expfile,  /* i: name of exposure map */
		    char         *evtsfile, /* i: name of event list   */
		    char         *evtext,    /* i: extension of events  */
		    WCS_str  *wcs                /* o: wcs data */
		    )
{

  /*dmDataset*    inputDs  = NULL;*/    /* event list data set        */
  /*dmBlock*      evtTable = NULL;*/    /* event list extension       */
  /*dmDescriptor* ot       = NULL;*/    /* ontime keyword descriptor  */

  /*dmBlockType   tabChk = dmTABLE; */

  long          ii;                 /* generic loop variable      */
  double         ontime;             /* ontime to use for exp. map */
  /*int           flag;*/               /* using native DM syntax=1 or not=0*/



  dmBlock *expMap = NULL;
  dmDescriptor *expData = NULL;

  struct DS_WCS_Data exp_wcs;
  float pointing[2];
  
  double radec[2];
  double pxpy[2];
  double lxly_d[2];
  long lxly_l[2];

  long nAxis;
  long *nAxes;
  long lo_left[2] = { 1, 1};

  dmDataType exp_type;
  double *exp_f = NULL;
  double *exp_d = NULL;

  
  memset( &exp_wcs, 0, sizeof( exp_wcs ));

  outputDebug( DEBUG_L3, 
	       "getting exposure map.\n", BEFORE ); /* output debugging msg */


  if  ( ( ds_strcmp_cis( expfile , "none" ) != 0 ) &&
	( ds_strcmp_cis( expfile, "" ) != 0 ) )
    {
      expMap = dmImageOpen( expfile );
      if ( expMap == NULL )
	{
	  /* ERROR */
	  err_msg( "ERROR: Could not open exposure map file %s",
		   expfile );
	  GlobalError |= BIT00;
	  return;
	}
      
      expData = dmImageGetDataDescriptor( expMap );

      
      ds_get_wcs_data( expData, pointing, 1, 1, &exp_wcs, NULL);

      /* read in exposure map */
      
      exp_type = dmGetDataType( expData );
      
      nAxis =  dmGetArrayDimensions( expData, &nAxes);
      
      
      switch ( exp_type )
	{
	case dmFLOAT:
	  exp_f = (double *)calloc( nAxes[0]*nAxes[1] , sizeof( double ));
	  dmImageDataGetSubArray_d( expData, lo_left, nAxes, exp_f );
	  break;
	  
	case dmDOUBLE:
	  exp_d = (double *)calloc( nAxes[0]*nAxes[1], sizeof( double ));
	  dmImageDataGetSubArray_d( expData, lo_left, nAxes, exp_d );
	  break;
	  
	default:
	  err_msg( "ERROR: Unsupported datatype for exposure mape image (double or double only)\n" );
	  GlobalError |= BIT00;
	  return;
	}
      

  
      /* populate exposure time tag of events */
      
      
      for ( ii=0; ii<numEvts; ii++ )
	{
	  
	  pxpy[0] = events[ii].coord.x;
	  pxpy[1] = events[ii].coord.y;
	  
	  
	  /* take events from physical to sky */
	  
	  if ( wcs->ngrp == 1 )
	    {
	      if ( wcs->wcs_descriptor[0] != NULL )
		dmCoordCalc_d( wcs->wcs_descriptor[0], &pxpy[0], &radec[0] );
	    }
	  else
	    {
	      if ( wcs->wcs_descriptor[0] != NULL )
		dmCoordCalc_d( wcs->wcs_descriptor[0], &pxpy[0], &radec[0] );
	      if ( wcs->wcs_descriptor[1] != NULL )
		dmCoordCalc_d( wcs->wcs_descriptor[1], &pxpy[1], &radec[1] );
	    }
	  
	  /* take sky down back to logical in expmap */
	  
	  if ( exp_wcs.axis_groups == 1 )
	    {
	      dmCoordInvert_d( exp_wcs.world[0], &radec[0], &pxpy[0] );
	      dmCoordInvert_d( exp_wcs.physical[0], &pxpy[0], &lxly_d[0] );
	    }
	  else
	    {
	      dmCoordInvert_d( exp_wcs.world[0], &radec[0], &pxpy[0] );
	      dmCoordInvert_d( exp_wcs.world[1], &radec[1], &pxpy[1] );
	      
	      dmCoordInvert_d( exp_wcs.physical[0], &pxpy[0], &lxly_d[0] );
	      dmCoordInvert_d( exp_wcs.physical[1], &pxpy[1], &lxly_d[1] );
	    }
	  
	  
	  lxly_l[0] = lxly_d[0]-1;
	  lxly_l[1] = lxly_d[1]-1;
	  
	  /* check range */
	  
	  if ( ( lxly_l[0] < 1 ) || ( lxly_l[0] > nAxes[0] ) ||
	       ( lxly_l[1] < 1 ) || ( lxly_l[1] > nAxes[1] )  )
	    {
	      ontime = 0.0 ;


	      if ( GlobalDebugLevel >= DEBUG_L4 )
		{
		  err_msg( "WARNING:  Pixel (%f,%f) is outside the range of the exposure map", events[ii].coord.x, events[ii].coord.y );
		}

	    }
	  else
	    {
	      switch ( exp_type )
		{
		case dmFLOAT:
		  ontime = exp_f[ lxly_l[0] + lxly_l[1] * nAxes[0] ];
		  break;
		  
		case dmDOUBLE:
		  ontime = exp_d[ lxly_l[0] + lxly_l[1] * nAxes[0] ];
		  break;
		  
		default:
		  ontime = 0.0;
		  break;
		}
	    }
	  
	  
	  events[ii].exp = ontime;
	  
	}
      
      
      free( nAxes );
      dmImageClose( expMap );
      
    } /* end if file is not "none" */
  else
    {
      
      if ( GlobalDebugLevel >= DEBUG_L2 )
	{
	  fprintf( GlobalDebugFile, "DEBUG: No exposure map provided.  Using time=1.0.\n");
	}
      
      for ( ii=0; ii<numEvts; ii++ )
	{
	  events[ii].exp = 1.0;
	}
      
    }
  


  outputDebug( DEBUG_L3, 
	       "getting exposure map.\n", AFTER ); /* output debugging msg */
  
  return;
  
}



/* 
   +-----------------------------------
   +
   + Get datamodel version number
   +
   +-----------------------------------
   */
void getDmVersion(
		  char **version     /* datamodel version number */
		  )
{
  *version = (char * ) calloc( SZ_FNAME, sizeof(char) );
  dmGetVersion( *version, SZ_FNAME);

  return;
}














/*
  +-------------------------------
  +
  + THIS DOESN'T WORK.  STUB to try to read in images
  + as well as event lists.
  +
  +-------------------------------
  */

void getEventsImg ( 
		   struct Event  **events,   /* array of events      */
		   long           *numEvts,  /* number of events     */
		   dmBlock *image, /* input file name      */
		   int             imgmin,   /* min. no events to consider*/
		   WCS_str  *wcs                /* o: wcs data */
		   )
{

  /*  dmBlock*       image  = NULL;*/  /* image  block */
  dmDescriptor*  data   = NULL;  /* data descriptor */
  /*  dmDescriptor*  xaxis  = NULL; */
  /*dmDescriptor*  yaxis  = NULL; */

  dmDataType     dType;

  long nx, ny, pixel;
  long ii, jj;
  long val, *axes;
  long loLeft[2], upRite[2];

  long  *val_l = NULL;
  short *val_s = NULL;
  unsigned char  *val_b = NULL;



  /*double pointing[2]; */

  double lxly[2], pxpy[2], wxwy[2];

  
  memset( &wcs_data, 0, sizeof(wcs_data));
  
  /* image = dmImageOpen(evtsfile);/
  if ( image == NULL )
    {
      err_msg( "ERROR: could not open image: %s\n", evtsfile);
      GlobalError = ( GlobalError | BIT02 );
      return;
    }
  */
  
  data  = dmImageGetDataDescriptor(image);
  if ( data == NULL )
    {
      err_msg("ERROR: could not open image: \n");
      GlobalError = ( GlobalError | BIT02 );
      return;
    }

  dType = dmGetDataType( data );


  ii    = dmGetArrayDimensions( data, &axes);
  nx = axes[0];
  ny = axes[1];



  /* GET WCS STUFF */

  if (!ds_get_wcs_data( data, NULL, nx/2, ny/2, &wcs_data, NULL ))
    wcs_present=TRUE;
  else
    wcs_present=FALSE;

  wcs->ngrp = wcs_data.axis_groups;

  wcs->wcs_descriptor = ( dmDescriptor **)calloc( wcs->ngrp, 
						  sizeof( dmDescriptor*));

  
  wcs->wcs_descriptor[0] = wcs_data.world[0];


  switch ( dType )
    {
    case dmBYTE:
      val_b = (unsigned char *) calloc (nx*ny,sizeof(unsigned char));
      if ( val_b == NULL )
	{
	  err_msg( "ERROR: could not allocate memory for image\n");
	  GlobalError = ( GlobalError | BIT02 ); /* CHECK */
	  return;
	}
      break;
      
    case dmSHORT:
      val_s = (short *) calloc (nx*ny,sizeof(short));
      if ( val_s == NULL )
	{
	  err_msg( "ERROR: could not allocate memory for image\n");
	  GlobalError = ( GlobalError | BIT02 ); /* CHECK */
	  return;
	}
      break;
      
    case dmLONG:
      val_l = (long *) calloc (nx*ny,sizeof(long));
      if ( val_l == NULL )
	{
	  err_msg( "ERROR: could not allocate memory for image\n");
	  GlobalError = ( GlobalError | BIT02 ); /* CHECK */
	  return;
	}
      break;

    default:
      err_msg( "ERROR: unsupported image data type (byte, short, long, only)\n");
      GlobalError = ( GlobalError | BIT02 ); /* CHECK */
      return;
    }




  loLeft[0]=1;
  loLeft[1]=1;
  upRite[0]=nx;
  upRite[1]=ny;



  switch ( dType )
    {
    case dmBYTE:
      dmImageDataGetSubArray_ub( data, loLeft, upRite, val_b );
      break;
    case dmSHORT:
      dmImageDataGetSubArray_s( data, loLeft, upRite, val_s );
      break;
    case dmLONG:
      dmImageDataGetSubArray_l( data, loLeft, upRite, val_l );
      break;
    default:
      err_msg( 
	       "ERROR: unsupported image data type (byte, short, long, only)\n");
      GlobalError = ( GlobalError | BIT02 ); /* CHECK */
      return;
    }


  *numEvts  = 0;
  (*events) = (struct Event *) calloc(4000, sizeof(struct Event) );
  if (*events == NULL )
    {
      err_msg( "ERROR: could not allocate memory for events\n");
      GlobalError = ( GlobalError | BIT02 ); /* CHECK */
      return;
    }


  for ( ii=0; ii<nx; ii++)
    {
      for ( jj=0; jj<ny; jj++)
	{
	  
	  pixel = ii + (nx *jj);
	  
	  switch ( dType )
	    {
	    case dmBYTE:
	      val = (long )val_b[pixel];
	      break;
	    case dmSHORT:
	      val = (long )val_s[pixel];
	      break;
	    case dmLONG:
	      val = val_l[pixel];
	      break;
	    default:
	      val = -999;
	      break;
	    }

	  
	  if ( val >= imgmin )
	    {
	      
	      lxly[0]=ii+1;
	      lxly[1]=jj+1;

	      if ( wcs_present ) 
		{
		  ds_coord_calc_d( &wcs_data, &lxly[0], 
				   &pxpy[0], &wxwy[0], NULL );
		  
		  (*events+(*numEvts))->coord.x = pxpy[0];
		  (*events+(*numEvts))->coord.y = pxpy[1];
		}
	      else
		{
		  (*events+(*numEvts))->coord.x = ii+1;
		  (*events+(*numEvts))->coord.y = jj+1;
		}
	      (*events+(*numEvts))->multi   = val;
	      (*events+(*numEvts))->area    = 0.0;
	      (*events+(*numEvts))->exp     = 0.0;
	      (*events+(*numEvts))->src_no  = 0;
	      (*events+(*numEvts))->tris    = NULL;
	      (*events+(*numEvts))->superTris = NULL;

	      (*numEvts) += 1;
	      
	      if (( ((*numEvts) % 4000 ) == 0) && ( (*numEvts) > 0 ) )
		{
		  (*events) = ( struct Event *) realloc( (*events),
							 ((*numEvts)/4000 + 1)*4000*
							 sizeof(struct Event));
		  
		  if ( *events ==NULL )
		    {

		      err_msg(
			      "ERROR: could not allocate enought space for events\n");
		      GlobalError = ( GlobalError | BIT02 ); /* CHECK */
		      return;
		    }


	      
		} /* end if need to alloc more events */
	    } /* end if pixel > min */
	} /* end for jj */
    } /* end for ii */
  
  if (val_l ) free(val_l);
  if (val_s ) free(val_s);
  if (val_b ) free(val_b);
  free(axes);
  

  wcs->parentBlock = image ;

  /*
    dmImageClose(image);
    */

  return;
  
}
