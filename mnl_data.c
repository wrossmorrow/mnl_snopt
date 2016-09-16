/*
 *  mnl_data.c
 *  
 *
 *  Created by W. Ross Morrow on 10/23/13.
 *  Copyright 2013 Iowa State University. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>

// #include <vecLib/cblas.h>
// #include <vecLib/clapack.h>
#include <Accelerate/Accelerate.h>

#include "mnl_macros.h"
#include "mnl_data.h"

#include "cholmod.h"
#include "SuiteSparseQR_C.h"

#define BUFFER_SIZE 1024

void mnl_data_new( MNL_DATA * data ) 
{
	
	data->on = 0;
	
	data->M  = 0;
	data->ml = NULL;
	data->Jm = NULL;
	data->J  = 0;
	
	data->K  = 0;
	
	data->B  = 0;
	data->Bo = NULL;
	data->BN = 0;
	
	data->D  = 0;
	data->Ld = NULL;
	data->L  = 0;
	
	data->XC = NULL;
	data->XD = NULL;
	
	data->s  = NULL;
	
	data->og = 'n';
	data->so = NULL;
	
}

void mnl_data_alloc( MNL_DATA * data )
{
	if( data == NULL ) { return; }
	if( data->on == 1 ) { mnl_data_free( data ); }
	
	
}

void mnl_data_free( MNL_DATA * data )
{
	int m;
	
	if( data == NULL ) { return; }
	
	if( data->K > 0 ) { 
		for( m = 0 ; m < data->M ; m++ ) { free(data->XC[m]); } 
		free(data->XC); data->XC = NULL;
		data->K = 0;
	}
	
	if( data->D > 0 ) { 
		for( m = 0 ; m < data->M ; m++ ) { free(data->XD[m]); } 
		free(data->XD); data->XD = NULL;
		free(data->Ld); data->Ld = NULL;
		data->D = 0;
		data->L = 0;
	}
	
	if( data->B > 0 ) { 
		for( m = 0 ; m < data->M ; m++ ) { 
			free(data->Bo[m]);
			free(data->XB[m]);
		} 
		free(data->Bo);	 data->Bo  = NULL;
		free(data->BoN); data->BoN = NULL;
		free(data->XB);	 data->XB  = NULL;
		data->B = 0;
		data->BN = 0;
	} 
	
	for( m = 0 ; m < data->M ; m++ ) { free(data->s[m]); }
	free(data->s); data->s = NULL;
	
	if( data->og == 'y' ) { free(data->so); data->so = NULL; }
	
	free(data->ml); data->ml = NULL;
	free(data->Jm); data->Jm = NULL;
	
	data->M = 0;
	data->J = 0;
	
	data->on = 0;
	
}

// column codes: 'c': characteristic, 'b': binary, 'd': dummy, 'm': market index, 's': shares, 'i': ignore

int readerror() { printf("read error...\n"); return -3; }

int mnl_data_read(FILE * dfp,			// file to read from
				int headrows ,		// header rows to skip
				char coderow ,		// "y" or "n" flag denoting whether id codes are in the first row read
				int colnum ,		// number of columns given (ignored if coderow == 'y', needed if coderow === 'n')
				char * coltype ,	// colnum-vector of column types (ignored if coderow == 'y', needed if coderow === 'n')
				MNL_DATA * data)	// data file to read into
{
	int C, I; // number of rows and columns (ignored columns)
	
	int c, r; // row/column indices
	
	int k, b, d, j;
	
	int m, ml, newml;
	
	double * XCt;
	int * XDt;
	int Bot;
	int * XBt;
	double s;
	double * sharesum;
	
	int hass, hasm; // has share data, market identifiers
	
	char * ct; // internal column types
	
	char line[BUFFER_SIZE];
	char * tok;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * ARGUMENT CHECKS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	if( dfp  == NULL ) { return -1; }
	if( data == NULL ) { return -1; }
	if( coderow == 'n' && ( colnum <= 0 || coltype == NULL ) ) { return -1; }
	
	// ensure file start
	rewind( dfp );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * PRE-PROCESS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	hass = 0; hasm = 0; data->K = 0; data->B = 0; data->D = 0; I = 0;
	
	if( coderow == 'n' ) {
		
		C = colnum;
		ct = coltype;
		
		for( c = 0 ; c < C ; c++ ) {
			switch( ct[c] ) {
				case 'c': (data->K)++; break;
				case 'b': (data->B)++; break;
				case 'd': (data->D)++; break;
				case 'm': hasm = 1;	   break;
				case 's': hass = 1;    break;
				case 'i': I++;		   break;
				default: break;
			}
		}
		
		// skip header rows
		for( r = 0 ; r < headrows ; r++ ) {
			if( fgets( line , BUFFER_SIZE , dfp ) == NULL ) { return readerror(); }
		}
		
	} else { // first (read) row has column codes
		
		// we first have to get number of columns
		
		// skip header rows
		for( r = 0 ; r < headrows ; r++ ) {
			if( fgets( line , BUFFER_SIZE , dfp ) == NULL ) { return readerror(); }
		}
		
		// read through first row to identify total number of columns
		if( fgets( line , BUFFER_SIZE , dfp ) == NULL ) { return readerror(); }
		else {
			C = 0;
			tok = strtok( line, "," ); // first csv-delimited token
			while( tok != NULL ) {
				C++; // increment column counter
				tok	= strtok( NULL , "," ); // get next csv-delimited token
			}
		}
		
		// now we have to read and store the column types
		
		// allocate space for column type storage
		ct = (char *)calloc( C , sizeof(char) );
		
		// rewind file to start at the beginning
		rewind( dfp );
		
		// skip header rows again
		for( r = 0 ; r < headrows ; r++ ) {
			if( fgets( line , BUFFER_SIZE , dfp ) == NULL ) { return readerror(); }
		}
		
		// read through first row to identify column types
		if( fgets( line , BUFFER_SIZE , dfp ) == NULL ) { return readerror(); }
		else {
			c = 0; tok = strtok( line, "," );
			while( tok != NULL ) {
				while( isspace(tok[0]) ) { tok++; } // trim off any leading whitespace
				switch( tok[0] ) { // store code of first character in token
					case 'c': ct[c] = 'c'; (data->K)++;	break;
					case 'b': ct[c] = 'b'; (data->B)++;	break;
					case 'd': ct[c] = 'd'; (data->D)++;	break;
					case 'm': ct[c] = 'm'; hasm = 1;	break;
					case 's': ct[c] = 's'; hass = 1;	break;
					case 'i': ct[c] = 'i'; I++;			break;
					default: break;
				}
				tok	= strtok( NULL , "," ); c++;
			}
		}
		
	}
	
	// if we don't have shares column, don't continue
	if( hass == 0 ) { return -2; }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// continue reading to get:
	// 
	//		number of levels per (non-binary) dummy variable
	//		number of markets (if identifier specified)
	//		total number of observations
	//
	
	// first, allocate level sizes
	data->Ld = ( int * )calloc( data->D , sizeof(int) );
	data->ml = NULL;
	
	data->M = 0; data->J = 0; data->L = 0;
	
	while( fgets( line , BUFFER_SIZE , dfp ) != NULL ) { 
		
		(data->J)++; // increment total number of observations
		
		d = 0; // zero out dummy variable counter
		c = 0; // zero out column counter
		
		tok = strtok( line, "," ); // first csv-delimited token
		while( tok != NULL ) {
			
			while( isspace(tok[0]) ) { tok++; } // trim off any leading whitespace from token
			
			switch( ct[c] ) {
				case 'd': 
					(data->Ld)[d] = MAX( (int)strtol(tok,NULL,10) , (data->Ld)[d] );
					d++;
					break;
				case 'm': 
					ml = (int)strtol(tok,NULL,10); // read market label
					newml = 1; 
					for( m = 0 ; m < data->M ; m++ ) { 
						if( (data->ml)[m] == ml ) { newml = 0; break; } 
					}
					if( newml == 1 ) { // new market label
						(data->M)++;
						data->ml = (int*)realloc( data->ml , (data->M)*sizeof(int) );
						(data->ml)[m] = ml;
					}
					break;
				default: break;
			}
			
			tok	= strtok( NULL , "," ); // get next csv-delimited token
			c++; // increment column counter
			
		}
		
	}
	
	// finish computing L, total number of dummy variable levels
	for( d = 0 ; d < data->D ; d++ ) { data->L += (data->Ld)[d]; }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// now we have market labels... we might want to re-arrange to
	// preserve some sort of order to the markets based on the labeling indices
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// re-read to get: 
	//
	//		number of products per market (if identifier specified)
	//
	
	// rewind to file start
	rewind(dfp);
	
	// skip header rows
	for( r = 0 ; r < headrows ; r++ ) {
		if( fgets( line , BUFFER_SIZE , dfp ) == NULL ) { return readerror(); }
	}
	
	// skip code row, if it exists
	if( coderow == 'y' ) {
		if( fgets( line , BUFFER_SIZE , dfp ) == NULL ) { return readerror(); }
	}
	
	// allocate space ...
	data->Jm = ( int * )calloc( data->M , sizeof(int) );
	
	XCt = (double *)calloc( data->K , sizeof(double) );
	XDt = (int    *)calloc( data->D , sizeof(int   ) );
	XBt = NULL;
	
	if( data->K > 0 ) {
		data->XC = ( double ** )calloc( data->M , sizeof(double *) );
		for( m = 0 ; m < data->M ; m++ ) { (data->XC)[m] = NULL; }
	}
	
	if( data->D > 0 ) {
		data->XD = ( int ** )calloc( data->M , sizeof(int *) );
		for( m = 0 ; m < data->M ; m++ ) { (data->XD)[m] = NULL; }
	}
	
	if( data->B > 0 ) {
		data->Bo = ( int ** )calloc( data->M , sizeof(int *) );
		data->XB = ( int ** )calloc( data->M , sizeof(int *) );
		for( m = 0 ; m < data->M ; m++ ) { 
			(data->Bo)[m] = NULL; (data->XB)[m] = NULL; 
		}
		data->BoN = ( int * )calloc( data->M , sizeof(int) );
	}
	data->BN = 0;
	
	data->s = ( double ** )calloc( data->M , sizeof(double *) );
	for( m = 0 ; m < data->M ; m++ ) { (data->s)[m] = NULL; }
	
	sharesum = (double *)calloc( data->M , sizeof(double) );
	
	// read through and parse market data
	while( fgets( line , BUFFER_SIZE , dfp ) != NULL ) {
		
		// zero out ... lots of things
		Bot = 0; XBt = NULL; c = 0; k = 0; b = 0; d = 0; s = 0.0;
		
		tok = strtok( line , "," ); // first csv-delimited token
		while( tok != NULL ) {
			
			while( isspace(tok[0]) ) { tok++; } // trim off any leading whitespace from token
			switch( ct[c] ) {
					
				case 'c': XCt[k++] = strtod(tok,NULL); break;
					
				case 'b': // read "on" binary variables into a running list
					if( (int)strtol(tok,NULL,10) == 1 ) {
						Bot++; // increment this product observation's number of "on" binaries
						XBt = ( int * )realloc( XBt , Bot * sizeof(int) ); // resize list
						XBt[ Bot-1 ] = b; // add this binary variable index to the running list
					}
					b++;
					break;
					
				case 'd': XDt[d++] = (int)strtol(tok,NULL,10) - 1; break; // read in in C-style indexing
					
				case 'm':
					ml = (int)strtol(tok,NULL,10); // read market label
					m = 0; while( ml != (data->ml)[m] ) { m++; } // find ordered market index
					break;
					
				case 's': 
					s = strtod(tok,NULL); 
					if( s < 0.0 ) { 
						printf("WARNING: negative shares not allowed. Ignoring this share.\n");
						s = 0.0;
					}
					break;
					
				default: break;
			}
			
			tok	= strtok( NULL , "," ); // get next csv-delimited token
			c++; // increment column counter
			
		}
		
		// printprodvec( data->K , XCt , data->B , Bot , XBt , data->D , XDt ); printf("\n");
		
		// now, copy this data into main arrays
		
		// increment product index and number of products in this market
		(data->Jm)[m]++;
		
		// reallocate arrays
		(data->XC)[m] = (double *)realloc( (data->XC)[m] , (data->K)*((data->Jm)[m]) * sizeof(double) );
		(data->XD)[m] = (int    *)realloc( (data->XD)[m] , (data->D)*((data->Jm)[m]) * sizeof(int   ) );
		(data->Bo)[m] = (int    *)realloc( (data->Bo)[m] ,           ((data->Jm)[m]) * sizeof(int   ) );
		(data->s )[m] = (double *)realloc( (data->s )[m] ,           ((data->Jm)[m]) * sizeof(double) );
		
		// overwrite c with the right index for XC(m), XD(m) matrices
		c = (data->Jm)[m] - 1;
		
		for( k = 0 ; k < data->K ; k++ ) { ((data->XC)[m])[ (data->K)*c + k ] = XCt[k]; }
		for( d = 0 ; d < data->D ; d++ ) { ((data->XD)[m])[ (data->D)*c + d ] = XDt[d]; }
		
		((data->Bo)[m])[c] = Bot; 
		(data->BoN)[m] += Bot; data->BN += Bot;
		(data->XB)[m] = ( int * )realloc( (data->XB)[m] , (data->BoN)[m] * sizeof(int) );
		for( b = 1 ; b <= Bot ; b++ ) { ((data->XB)[m])[ (data->BoN)[m]-b ] = XBt[Bot-b]; }
		
		((data->s)[m])[ (data->Jm)[m]-1 ] = s;
		sharesum[m] += s;
		
		/*
		printprodvec( data->K , (data->XC)[m] + (data->K)*c , 
					  data->B , ((data->Bo)[m])[c] , (data->XB)[m] + (data->BoN)[m] - ((data->Bo)[m])[c] , 
					  data->D , (data->XD)[m] + (data->D)*c ); 
		printf("\n");
		 */
		
		// free XBt... as this was allocated above and will be re-allocated in next loop
		free(XBt);
		
	}
	
	// check share sums to ensure normalization
	for( m = 0 ; m < data->M ; m++ ) { 
		
		if( sharesum[m] == 0.0 ) {
			printf("WARNING: market %i has zero shares.\n",m+1);
		}
		
		if( sharesum[m] != 1.0 ) { 
			for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
				((data->s)[m])[j] = ((data->s)[m])[j] / sharesum[m];
			}
		}
		
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * FREE DATA * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	if( coderow == 'y' ){ free(ct); }
	
	free( XCt );
	free( XDt );
	// XBt freed above
	
	free(sharesum);
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	return 1;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
}

void printprodvec( int K , double * xc , int B , int Bo , int * xb , int D , int * xd )
{	
	int k, b, d;
	
	if( K > 0 ) {
		printf("( %0.2f ",xc[0]); for( k = 1 ; k < K ; k++ ) { printf(", %0.2f ",xc[k]); } printf(") ");
	} else { printf("( ) "); }
	
	if( B > 0 && Bo > 0 ) {
		printf("( %i ",xb[0]+1); for( b = 1 ; b < Bo ; b++ ) { printf(", %i ",xb[b]+1); } printf(") ");
	} else { printf("( ) "); }
	
	if( D > 0 ) {
		printf("( %i ",xd[0]+1); for( d = 1 ; d < D ; d++ ) { printf(", %i ",xd[d]+1); } printf(") ");
	} else { printf("( ) "); }
}

void mnl_data_print( MNL_DATA * data )
{
	int m, j, d, Bo, bbase;
	double * xc;
	int * xb, * xd;
	
	if( data == NULL ) { return; }
	
	printf("Basic data file info: \n");
	printf("  %i characteristics\n",data->K);
	printf("  %i binary variables\n",data->B);
	printf("    %i total %con%c variables\n",data->BN,'"','"');
	printf("  %i dummy variables:\n",data->D);
	for( d = 0 ; d < data->D ; d++ ) {
		printf("    dummy %i has %i levels\n",d+1,(data->Ld)[d]);
	}
	printf("  %i markets:\n",data->M);
	for( m = 0 ; m < data->M ; m++ ) {
		printf("    market %i (%c%i%c) has %i products\n",m+1,'"',(data->ml)[m],'"',(data->Jm)[m]);
		bbase = 0;
		for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
			if( data->K > 0 ) { xc = (data->XC)[m] + (data->K)*j; }
			if( data->B > 0 ) { 
				Bo = ((data->Bo)[m])[j]; 
				if( Bo > 0 ) { xb = (data->XB)[m] + bbase; bbase += Bo; } 
			}
			if( data->D > 0 ) { xd = (data->XD)[m] + (data->D)*j; }
			printprodvec( data->K , xc , data->B , Bo , xb , data->D , xd ); 
			printf("\n");
		}
	}
	printf("  %i total product observations\n",data->J);
	/*
	for( m = 0 ; m < data->M ; m++ ) {
		printf("    market %i product shares: ( %0.3f ",m+1,100.0*((data->s)[m])[0]);
		for( j = 1 ; j < (data->Jm)[m] ; j++ ) {
			printf(", %0.3f ",100.0*((data->s)[m])[j]);
		}
		printf(")\n");
	}*/
	
	
	
	printf("\n");
}

// This function transforms a data structure into another data structure
// by mapping each product definition into another product definition. 
// The basic mapping is given by a function "map" having the following prototype: 
//
//		map( int K_s , double xc_s , int B_s , int   Bo_s , int * xb_s , int D_s , int * xd_s , 
//			 int K_d , double xc_d , int B_d , int * Bo_d , int * xb_d , int D_d , int * xd_d , 
//			 void * params )
//
// The arguments define the structure of the mapping between characteristics, binaries, 
// and dummies in the source ("_s") and destination ("_d") data structures. The "prm"
// pointer is a pointer to user data that will be passed to all calls of "map". Note that
// the number of "on" binary variables in the destination is not necessarily known; thus
// the argument Bo_d is a pointer, into which the number of on variables should be written, 
// and the array xb_d should be of length at least B_d. Only the first Bo_d[0] (on return)
// values will be written into. 
//
// The information about destination sizes is presumed to already be in place in the
// "dst" argument. That is, K_d, B_d, D_d as well as Ld_d and L are assumed to be defined in dest. This
// will result in array sizing decisions, so specification of these terms consistent
// with the mapping function provided is important to avoid memory errors. 
// 
int mnl_data_map( MNL_DATA * src , MNL_DATA * dst , MNL_DATA_MAP map , void * prm )
{
	int m, j, k, b, d;
	
	double * xc_d;
	int * xb_d;
	int * xd_d;
	int Bo;
	
	int sbbase;
	
	if( src == NULL ) { return 0; }
	if( dst == NULL ) { return 0; }
	if( map == NULL ) { return 0; }
	if( dst->K <= 0 && dst->B <= 0 && dst->D <= 0 ) { return 0; }
	
	dst->M = src->M; // copy number of markets
	dst->J = src->J; // copy total number of product observations
	
	dst->ml = (int *)calloc( dst->M , sizeof(int) ); // allocate market labels
	dst->Jm = (int *)calloc( dst->M , sizeof(int) ); // allocate number of products per market
	for( m = 0 ; m < src->M ; m++ ) {
		(dst->ml)[m] = (src->ml)[m]; // copy market labels
		(dst->Jm)[m] = (src->Jm)[m]; // copy number of products per market
	}
	
	dst->s = (double **)calloc( dst->M , sizeof(double *) ); // allocate market share pointers
	for( m = 0 ; m < dst->M ; m++ ) {
		(dst->s)[m] = (double *)calloc( (dst->Jm)[m] , sizeof(double) ); // allocate market shares
		cblas_dcopy( (dst->Jm)[m] , (src->s)[m] , 1 , (dst->s)[m] , 1 ); // copy market shares
	}
	
	if( src->og == 'y' ) {
		dst->og = 'y'; // set outside good flag
		dst->so = (double *)calloc( dst->M , sizeof(double) ); // allocate outside good shares
		cblas_dcopy( dst->M , src->so , 1 , dst->so , 1 ); // copy outside good shares
	} else { dst->og = 'n'; dst->so = NULL; }
	
	// now, transform characteristics, binaries and dummies for each product
	
	// allocate space for characteristics
	if( dst->K > 0 ) {
		dst->XC = (double **)calloc( dst->M , sizeof(double *) );
		for( m = 0 ; m < dst->M ; m++ ) {
			(dst->XC)[m] = (double *)calloc( (dst->K)*((dst->Jm)[m]) , sizeof(double) );
		}
	} else { dst->XC = NULL; }
	
	// allocate space for binaries
	if( dst->B > 0 ) {
		dst->Bo = (int **)calloc( dst->M , sizeof(int *) );
		dst->XB = (int **)calloc( dst->M , sizeof(int *) );
		for( m = 0 ; m < dst->M ; m++ ) {
			(dst->Bo)[m] = (int *)calloc( ((dst->Jm)[m]) , sizeof(int) );
			(dst->XB)[m] = NULL; // sized and allocated below
		}
		dst->BoN = (int *)calloc( dst->M , sizeof(int) );
	} else { dst->Bo = NULL; dst->BoN = NULL; dst->XB = NULL; xb_d = NULL; }
	dst->BN = 0;
	
	// allocate space for dummies
	dst->XD = (int **)calloc( dst->M , sizeof(int *) );
	if( dst->D > 0 ) { 
		for( m = 0 ; m < dst->M ; m++ ) {
			(dst->XD)[m] = (int *)calloc( (dst->D)*((dst->Jm)[m]) , sizeof(int) );
		}
	} else {  dst->XD = NULL; }
		
	// actual transformations
	
	xc_d = (double *)calloc( dst->K , sizeof(double) );
	xb_d = (int    *)calloc( dst->B , sizeof(int   ) );
	xd_d = (int    *)calloc( dst->D , sizeof(int   ) );
	
	for( m = 0 ; m < dst->M ; m++ ) {
		
		sbbase = 0;
		for( j = 0 ; j < (dst->Jm)[m] ; j++ ) {
			
			// use map to define destination characteristics, binaries, and dummies
			// from the source destination characteristics, binaries, and dummies
			
			map(src->K ,  (src->XC)[m] + (src->K)*j , 
				src->B , ((src->Bo)[m])[j]			, (src->XB)[m] + sbbase , 
				src->D ,  (src->XD)[m] + (src->D)*j , 
				dst->K ,  xc_d , 
				dst->B ,  &Bo , xb_d , 
				dst->D ,  xd_d , 
				prm );
			
			/*
			printprodvec(src->K ,  (src->XC)[m] + (src->K)*j , 
						 src->B , ((src->Bo)[m])[j]			, (src->XB)[m] + sbbase , 
						 src->D ,  (src->XD)[m] + (src->D)*j );
			printf(" -> ");
			printprodvec(dst->K ,  xc_d , 
						 dst->B ,  Bo , xb_d , 
						 dst->D ,  xd_d );
			printf("\n");
			/**/
			
			if( dst->K > 0 ) { 
				for( k = 0 ; k < dst->K ; k++ ) { ((dst->XC)[m])[(dst->K)*j+k] = xc_d[k]; }
			}
			
			if( dst->B > 0 ) {
				
				((dst->Bo)[m])[j] = Bo; // define number of "on" binary variables
				
				(dst->BoN)[m] += ((dst->Bo)[m])[j]; // increment number of "on" binary variables in this market
				dst->BN += ((dst->Bo)[m])[j]; // increment total number of "on" binary variables
				
				// set "on" binary variables in the destination data structure
				(dst->XB)[m] = (int *)realloc( (dst->XB)[m] , ((dst->BoN)[m])*sizeof(int) );
				for( b = 1 ; b <= ((dst->Bo)[m])[j] ; b++ ) {
					((dst->XB)[m])[ (dst->BoN)[m] - b ] = xb_d[ Bo - b ];
				}
				
			}
			
			if( dst->D > 0 ) { 
				for( d = 0 ; d < dst->D ; d++ ) { ((dst->XD)[m])[(dst->D)*j+d] = xd_d[d]; }
			}
			
			sbbase += ((src->Bo)[m])[j]; // increment source binary base by number of "on" binary variables in the source
			
		}
		
	}
	
	if( xc_d != NULL ) { free(xc_d); }
	if( xb_d != NULL ) { free(xb_d); }
	if( xd_d != NULL ) { free(xd_d); }
	
	return 1;
	
}

// this routine should check (1) the rank of the utility constraints and (2) how close
// the image of the constraints gets to a vector of ones
int mnl_data_check( MNL_DATA * data , int * rank , double * norm )
{
	
	int i, m, j, k, b, d, l, ibase, jbase, bbase, Abase;
	int b_v_start, p_v_start, g_v_start;
	
	int minrow, maxrow, mincol, maxcol;
	
	long int N , M, Annz;
	
	cholmod_common c;
	
	cholmod_triplet * At;
	double * Adata;
	long int * Arows , * Acols;
	
	cholmod_sparse * AT; 
    cholmod_dense  * y , * w;
	double * yx;
	
	SuiteSparseQR_C_factorization * QR;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	if( data == NULL ) { return 0; }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// define matrix sizes
	M = (long int)( data->D + data->J );			// number of rows
	N = (long int)( data->K + data->B + data->L );	// number of columns
	
	b_v_start = 0;
	p_v_start = b_v_start + data->K;
	g_v_start = p_v_start + data->B;
	
	Annz  = data->L;				// total of L terms in D dummy variable constraints
	Annz += (data->K)*(data->J);	// K characteristic terms in J utility equations
	Annz += data->BN;				// BN binary variable terms (one for each "on" variable)
	Annz += (data->D)*(data->J);	// D dummy terms in J utility equations (one nonzero per dummy, per equation)
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// start cholmod common object
	cholmod_l_start( &c );
	
	// allocate space for a "triplet" object to store matrix in
	At = cholmod_l_allocate_triplet( M , N , Annz , 0 , CHOLMOD_REAL , &c );
	if( At == NULL ) { 
		printf("Triplet form sparse matrix could not be constructed.\n"); 
		cholmod_l_finish( &c ); 
		return 0; 
	}
	
	// set number of non-zeros
	At->nnz = Annz;
	
	// 
	Adata = ( double   * )( At->x ); 
	Arows = ( long int * )( At->i ); 
	Acols = ( long int * )( At->j );
	
	Abase = 0;
	
	// dummy variable constraints (sum-to-one within variables)
	
	ibase = 0;
	for( d = 0 ; d < data->D ; d++ ) {
		for( l = 0 ; l < (data->Ld)[d] ; l++ ) {
			Arows[Abase] = d;
			Acols[Abase] = g_v_start + ibase + l;
			Adata[Abase] = 1.0;
			Abase++;
		}
		ibase += (data->Ld)[d];
	}
	
	// utility equation constraints: characteristic terms
	
	if( data->K > 0 ) { // don't do any of this if K == 0
		
		jbase = 0;
		for( m = 0 ; m < data->M ; m++ ) {
			
			// create linear characteristics part of utility equations with transformed coordinates
			for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
				for( k = 0 ; k < data->K ; k++ ) {
					Arows[Abase] = data->D + jbase;
					Acols[Abase] = b_v_start + k;
					Adata[Abase] = ((data->XC)[m])[ (data->K)*j + k ]; 
					// Adata[Abase] = ( ((data->XC)[m])[ (data->K)*j + k ] - XCL[k] ) / XCM[k] - 1.0;
					Abase++;
				}
				jbase++;
			}
		}
		
	}
	
	// utility equation constraints: "on" binary variable terms
	
	if( data->B > 0 ) { // don't do any of this if B == 0
		jbase = 0;
		for( m = 0 ; m < data->M ; m++ ) {
			bbase = 0;
			for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
				for( b = 0 ; b < ((data->Bo)[m])[j] ; b++ ) {
					Arows[Abase] = data->D + jbase;
					Acols[Abase] = p_v_start + ((data->XB)[m])[ bbase ];
					Adata[Abase] = 1.0;
					bbase++;
					Abase++;
				}
				jbase++;
			}
		}
	}
	
	// utility equation constraints: dummy variable terms
	
	if( data->D > 0 ) { // don't do any of this if D == 0
		jbase = 0;
		for( m = 0 ; m < data->M ; m++ ) {
			for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
				ibase = 0;
				// printf("  ");
				for( d = 0 ; d < data->D ; d++ ) {
					Arows[Abase] = data->D + jbase;
					Acols[Abase] = g_v_start + ibase + ((data->XD)[m])[ (data->D)*j + d ];
					// printf("%i (%i,%i), ",((data->XD)[m])[ (data->D)*j + d ],Arows[Abase]+1,Acols[Abase]+1);
					Adata[Abase] = 1.0;
					Abase++;
					ibase += (data->Ld)[d];
				}
				// printf("\n");
				jbase++;
			}
		}
	}
	 
	/*
	printf("A: \n");
	minrow = Arows[0]; maxrow = Arows[0];
	mincol = Acols[0]; maxcol = Acols[0];
	for( j = 0 ; j < Annz ; j++ ) {
		// printf("  A(%i,%i) = %0.6f\n",Arows[j]+1,Acols[j]+1,Adata[j]);
		minrow = MIN( minrow , Arows[j] ); maxrow = MAX( maxrow , Arows[j] );
		mincol = MIN( mincol , Acols[j] ); maxcol = MAX( maxcol , Acols[j] );
	}
	
	printf("  min row: %i, max row: %i\n",minrow,maxrow);
	printf("  min col: %i, max col: %i\n",mincol,maxcol);
	*/
	
	// convert At to a cholmod_sparse object
	AT = cholmod_l_triplet_to_sparse( At , Annz , &c );
	
	// free the triplet object
	cholmod_l_free_triplet( &At , &c );
	
	// make sure we were successful in creating AT
	if( AT == NULL ) { 
		
		printf("Cholmod form sparse matrix could not be constructed.\n"); 
		
		
		printf("A: \n");
		minrow = Arows[0]; maxrow = Arows[0];
		mincol = Acols[0]; maxcol = Acols[0];
		for( j = 0 ; j < Annz ; j++ ) {
			printf("  A(%li,%li) = %0.6f\n",Arows[j]+1,Acols[j]+1,Adata[j]);
			minrow = MIN( minrow , Arows[j] ); maxrow = MAX( maxrow , Arows[j] );
			mincol = MIN( mincol , Acols[j] ); maxcol = MAX( maxcol , Acols[j] );
		}
		printf("  min row: %i, max row: %i\n",minrow,maxrow);
		printf("  min col: %i, max col: %i\n",mincol,maxcol);
		
		
		cholmod_l_finish( &c ); 
		
		return 0; 
		
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// form QR factorization of AT
	QR = SuiteSparseQR_C_factorize( SPQR_ORDERING_DEFAULT , 1e-20 , AT , &c );
	
	// check for successful factorization
	if( QR == NULL ){ 
		printf("QR factorization failed.\n"); 
		cholmod_l_free_sparse( &AT , &c );
		cholmod_l_finish( &c );
		return 0; 
	}
	
	// store rank estimate from factorization
	if( rank != NULL ) { rank[0] = (int)( (c.SPQR_istat)[4] ); }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// construct a (cholmod) vector RHS y = [ 0 ; 1 ]
	y = cholmod_l_ones( M , 1 , CHOLMOD_REAL , &c ); yx = (double *)(y->x);
	for( d = 0 ; d < data->D ; d++ ) { yx[d] = 0.0; }
	
	// form w <- Q' y (space allocated inside routine?)
	w = SuiteSparseQR_C_qmult( SPQR_QTX , QR , y , &c );
	if( w != NULL ) { 
		
		// compute || w(N+1:M) ||_2. This is an estimate of the minimum distance between
		// the image of the characteristics-binaries-dummies matrix and the vector of all ones
		if( norm != NULL ) { norm[0] = cblas_dnrm2( M-N , (double *)(w->x) + N , 1 ); }
		
		// clean up memory
		
		SuiteSparseQR_C_free ( &QR , &c );
		cholmod_l_free_sparse( &AT , &c );
		cholmod_l_free_dense ( &y  , &c );
		cholmod_l_free_dense ( &w  , &c );
		
		cholmod_l_finish( &c );
		
		// return success code
		
		return 1;
		
	} else { 
		
		printf("Warning - w = Q'1 could not be formed.\n"); 
		
		// clean up memory
		
		SuiteSparseQR_C_free ( &QR , &c );
		cholmod_l_free_sparse( &AT , &c );
		cholmod_l_free_dense ( &y  , &c );
		cholmod_l_free_dense ( &w  , &c );
		
		cholmod_l_finish( &c );
		
		// return failure code
		
		return 0;
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
}

void mnl_model_new( MNL_MODEL * model )
{
	if( model == NULL ) { return; }
	
	model->on = 'n';
	
	model->K = 0;
	
	model->B = 0;
	
	model->D = 0;
	model->Ld = NULL;
	model->L = 0;
	
	model->cc = NULL;
	model->bc = NULL;
	model->dc = NULL;
	
	model->og = 'n';
	model->uc = 0.0;
	
}

void mnl_model_alloc( MNL_MODEL * model , int K , int B , int D , int * Ld )
{
	int d;
	
	if( model == NULL ) { return; }
	if( K == 0 && B == 0 && D == 0 ) { return; }
	if( D > 0 && Ld == NULL ) { return; }
	
	model->K = K;
	model->cc = ( double * )calloc( K , sizeof( double ) );
	
	model->B = B;
	model->bc = ( double * )calloc( B , sizeof(double) );
	
	model->D = D; 
	model->L = 0;
	model->Ld = ( int * )calloc( D , sizeof(int) );
	for( d = 0 ; d < D ; d++ ) { (model->Ld)[d] = Ld[d]; model->L += Ld[d]; }
	model->dc = ( double * )calloc( model->L , sizeof(double) );
	
	model->on = 'y';
}

void mnl_model_free( MNL_MODEL * model )
{
	if( model == NULL ) { return; }
	
	if( model->Ld != NULL ) { free(model->Ld); }
	if( model->cc != NULL ) { free(model->cc); }
	if( model->bc != NULL ) { free(model->bc); }
	if( model->dc != NULL ) { free(model->dc); }
	
	model->K = 0;
	model->B = 0;
	model->D = 0;
	model->L = 0;
	model->og = 'n';
	model->uc = 0.0;
	
	model->on = 'n';
}

int mnl_model_probabilities( MNL_MODEL * model , MNL_DATA * data , double ** PL , double * LL )
{
	int k, b, d, m, j, base, dbase, bbase;
	
	double * U, * Us, * E, * Es;
	
	if( model == NULL || data  == NULL ) { return 0; }
	
	if( model->K != data->K ) { return -1; }
	if( model->B != data->B ) { return -2; }
	if( model->D != data->D ) { return -3; }
	for( d = 0 ; d < model->D ; d++ ) {
		if( (model->Ld)[d] > (data->Ld)[d] ) {  return -4; }
	}
	if( model->og != data->og ) { return -5; }
	
	// allocate internal space
	
	U  = (double *)calloc( data->J , sizeof(double) );
	E  = (double *)calloc( data->J , sizeof(double) );
	Us = (double *)calloc( data->M , sizeof(double) );
	Es = (double *)calloc( data->M , sizeof(double) );
	
	// compute utilities
	
	base = 0; 
	for( m = 0 ; m < data->M ; m++ ) {
		bbase = 0; // reset binary variable indexer in each market
		for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
			
			// characteristic part
			U[base] = cblas_ddot( data->K , (data->XC)[m] + (data->K)*j , 1 , model->cc , 1 );
			
			// binary variable part
			// number of "on" binary variables for this product: (Bo[m])[j]
			// indices of "on" variables: (XB[m])(bbase,bbase+1,...)
			if( data->B > 0 ) {
				for( b = 0 ; b < ((data->Bo)[m])[j] ; b++ ) {
					U[base] += (model->bc)[ ((data->XB)[m])[ bbase++ ] ];
				}
			}
			
			// dummy variable part
			dbase = 0;
			for( d = 0 ; d < data->D ; d++ ) {
				U[base] += (model->dc)[ ((data->XD)[m])[ (data->D)*j + d] + dbase ];
				dbase += (data->Ld)[d];
			}
			
			Us[m] = MAX( U[base] , Us[m] );
			base++;
			
		}
	}
	
	/*
	printf("utilities: \n");
	base = 0;
	for( m = 0 ; m < data->M ; m++ ) {
		for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
			printf("  u(%i,%i) = %0.16f\n",m+1,j+1,U[base++]);
		}
	}
	 */
	
	// exponentiate shifted utilities
	
	base = 0;
	for( m = 0 ; m < data->M ; m++ ) {
		Es[m] = 0.0;
		if( Us[m] > 0.0 ) {
			for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
				E[base] = exp( U[base] - Us[m] ); Es[m] += E[base]; base++;
			}
		} else {
			for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
				E[base] = exp( U[base] ); Es[m] += E[base]; base++;
			}
		}
	}
	
	if( LL != NULL ) { // compute log likelihood
		
		base = 0; LL[0] = 0.0;
		for( m = 0 ; m < data->M ; m++ ) {
			LL[0] += cblas_ddot( (data->Jm)[m] , (data->s)[m] , 1 , U + base , 1 );
			LL[0] -= ( Us[m] + log( Es[m] ) );
			base += (data->Jm)[m];
		}
		
	}
	
	if( PL != NULL ) { // compute probabilities
		
		base = 0;
		for( m = 0 ; m < data->M ; m++ ) {
			if( PL[m] != NULL ) { 
				if( (data->Jm)[m] == 1 ) { (PL[m])[0] = 1.0; }
				else{ 
					for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
						(PL[m])[j] = E[base] / Es[m]; base++;
					}
				}
			}
		}
		
	}
	
	free(U);
	free(Us);
	free(E);
	free(Es);
	
	return 1;
	
}

// compute fischer information matrix
int mnl_model_fisinfmat( MNL_MODEL * model , MNL_DATA * data , double ** PL , double * LL )
{
	
	// 0. consistency checks?
	
	// 1. compute probabilities
	
	// 2. form negative hessian matrix of the log likelihood function
	
	base = 0;
	for( m = 0 ; m < data->M ; m++ ) {
		
		
		if( PL[m] != NULL ) { 
			if( (data->Jm)[m] == 1 ) { ; }
			else{ 
				for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
					
					
					
				}
			}
		}
		
		
	}
	
}

void mnl_model_print( MNL_MODEL * model )
{
	int k, b, d, l, dbase;
	
	if( model == NULL ) { return; }
	
	printf("K: %i\n",model->K);
	if( model->K > 0 ) { 
		printf( "  %0.16f " , (model->cc)[0] );  
		for( k = 1 ; k < model->K ; k++ ) { 
			printf( ", %0.16f " , (model->cc)[k] );  
		}
		printf("\n");
	}
	
	printf("B: %i\n",model->B);
	if( model->B > 0 ) { 
		printf( "  %0.16f " , (model->bc)[0] );  
		for( b = 1 ; b < model->B ; b++ ) { 
			printf( ", %0.16f " , (model->bc)[b] );  
		}
		printf("\n");
	}
	
	printf("D: %i\n",model->D);
	if( model->D > 0 ) { 
		dbase = 0;
		for( d = 0 ; d < model->D ; d++ ) { 
			printf( "  %0.16f " , (model->dc)[dbase++] ); 
			for( l = 1 ; l < (model->Ld)[d] ; l++ ) {
				printf( ", %0.16f " , (model->dc)[dbase++] );  
			}
			printf("\n");
		}
	}
	
}

void mnl_model_fprint( MNL_MODEL * model , FILE * fp )
{
	int k, b, d;
	
	if( model == NULL || fp == NULL ) { return; }
	
	fprintf( fp , "%i , %i , %i " , model->K , model->B , model->D ); 
	if( model->D > 0 ) { 
		fprintf( fp , ", %i " , (model->Ld)[0] );  
		for( d = 1 ; d < model->D ; d++ ) { 
			fprintf( fp , ", %i " , (model->Ld)[d] );  
		}
	}
	if( model->K > 0 ) { 
		fprintf( fp , ", %0.16f " , (model->cc)[0] );  
		for( k = 1 ; k < model->K ; k++ ) { 
			fprintf( fp , ", %0.16f " , (model->cc)[k] );  
		}
	}
	if( model->B > 0 ) { 
		fprintf( fp , ", %0.16f " , (model->bc)[0] );  
		for( b = 1 ; b < model->B ; b++ ) { 
			fprintf( fp , ", %0.16f " , (model->bc)[b] );  
		}
	}
	if( model->D > 0 ) { 
		fprintf( fp , ", %0.16f " , (model->dc)[0] );  
		for( d = 1 ; d < model->L ; d++ ) { 
			fprintf( fp , ", %0.16f " , (model->dc)[d] );  
		}
	}
}
