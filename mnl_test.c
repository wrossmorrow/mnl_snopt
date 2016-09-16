/*
 *  mnl_test.c
 *  
 *
 *  Created by W. Ross Morrow on 10/23/13.
 *  Copyright 2013 Iowa State University. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "mnl_macros.h"
#include "mnl_data.h"
#include "mnl_snopt.h"
	
#ifdef __cplusplus
	}
#endif

#define BIGGEST_INT 276447232

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int getnext( int N , int * L , int * U , int * x ) 
{
	int i , j;
	
	for( i = 0 ; i < N ; i++ ) { 
		if( x[i] < U[i] ) { 
			x[i]++; for( j = 0 ; j <= i-1 ; j++ ) { x[j] = L[j]; } return 1;
		}
	}

	return 0;
}

// Calling codes:
//
// Characteristics:
//
//	price:			1, 2			(2)
//	op. cost:		0, 1, 2, 3, 4	(5)
//	acceleration:	0, 1, 2, 3, 4	(5)
//	size:			0, 1, 2, 3, 4	(5)
//	style:			0, 1			(2)
//	previous share: 0, 1, 2			(3)
//
// Binaries:
// 
//	luxury:			0, 1			(2)
//	transmission:	0, 1			(2)
//
// Dummies:
// 
//	brand:			0, 1			(2)
//	class:			0, 1			(2)
//
static int codeN = 10;
static int codeL[10] = { 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
static int codeU[10] = { 2 , 4 , 4 , 4 , 1 , 2 , 1 , 1 , 1 , 1 };
static int codes[10] = { 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }; // default starting code (== to codeL)
// static int codes[10] = { 2 , 4 , 3 , 4 , 1 , 0 , 0 , 1 , 1 , 0 }; // failed without scaling
// static int codes[10] = { 2 , 1 , 3 , 4 , 1 , 0 , 1 , 1 , 1 , 0 }; // 2134101110


static int K = 6; // price, op. cost, acc, size, style, prev. shr
static int B = 2; // luxury, transmission
static int D = 2; // brand, class
static int Ld[2] = { 37 , 10 }; // limits on brand and class

void definesizes( MNL_DATA * data , int * codes )
{
	int k, b, d, base;
	
	base = 0;
	
	data->K = 0;
	for( k = 0 ; k < K ; k++ ) {
		switch( codes[base++] ) {
			case  0: break;
			default: data->K++; break;
		}
	}
	
	data->B = 0;
	for( b = 0 ; b < B ; b++ ) {
		switch( codes[base++] ) {
			case  0: break;
			default: data->B++; break;
		}
	}
	
	data->D = 0; data->Ld = NULL; data->L = 0;
	for( d = 0 ; d < D ; d++ ) {
		switch( codes[base++] ) {
			case  0: break;
			default:
				data->D++;
				data->Ld = (int *)realloc( data->Ld , (data->D)*sizeof(int) );
				data->Ld[ data->D - 1 ] = Ld[d];
				data->L += Ld[d];
				break;
		}
	}
	
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int map( int K_s , double * xc_s , int B_s , int   Bo_s , int * xb_s , int D_s , int * xd_s , 
		 int K_d , double * xc_d , int B_d , int * Bo_d , int * xb_d , int D_d , int * xd_d , 
		 void * params )
{
	int i, k, b, d, base, bp, bbase;
	
	int * codes;
	
	double mpg, hp, wgt, len, wid, hgt, gpr, psh, pri;
	
	// presume we have been passed codes as params
	codes = (int *)params;
	
	mpg = xc_s[0];
	hp  = xc_s[1];
	wgt = xc_s[2];
	len = xc_s[3];
	wid = xc_s[4];
	hgt = xc_s[5];
	gpr = xc_s[6];
	psh = xc_s[7];
	pri = xc_s[8];
	
	base = 0;
	
	// characteristics
	
	k = 0;
	
	// price
	switch( codes[base++] ) {
		case  1: xc_d[k++] = pri;		break; // 
		case  2: xc_d[k++] = log(pri);	break;
		default: break;
	}
	
	// operating cost
	switch( codes[base++] ) {
		case  1: xc_d[k++] = gpr / mpg; break;
		case  2: xc_d[k++] = mpg / gpr; break;
		case  3: xc_d[k++] = mpg;		break; // 
		case  4: xc_d[k++] = 1.0 / mpg; break;
		default: break;
	}
	
	// acceleration
	switch( codes[base++] ) {
		case  1: xc_d[k++] = hp / wgt;	break;
		case  2: xc_d[k++] = wgt / hp;	break;
		case  3: xc_d[k++] = exp ( - 0.00275 * pow( hp/wgt , -0.776 ) ); break; // 
		case  4: xc_d[k++] = hp;		break;
		default: break;
	}
	
	// size
	switch( codes[base++] ) {
		case  1: xc_d[k++] = len;		break;
		case  2: xc_d[k++] = wid;		break;
		case  3: xc_d[k++] = len - wid; break;
		case  4: xc_d[k++] = len * wid; break; // 
		default: break;
	}
	
	// style
	switch( codes[base++] ) {
		case  1: xc_d[k++] = len * wid / hgt; break;
		default: break;
	}
	
	// previous share
	switch( codes[base++] ) {
		case  1: xc_d[k++] = psh;		break; // 
		case  2: xc_d[k++] = log(psh);	break;
		default: break;
	}
	
	// binaries... this has to be dealt with differently. 
	
	Bo_d[0] = 0; 
	
	// luxury and transmission
	
	bbase = 0;
	
	for( b = 0 ; b < B ; b++ ) {
		
		if( codes[base++] == 1 ) { // if this variable is included 
			
			// search through "on" binary variables in source for this 
			for( bp = 0 ; bp < Bo_s ; bp++ ) { 
				if( xb_s[bp] == b ) { // if this component is "on" in source (and coded)
					xb_d[Bo_d[0]] = bbase; Bo_d[0]++;
				}
			}
			
			// increment binary variable counter
			bbase++;
			
		}
		
	}
	
	// dummies
	
	d = 0;
	
	// brand
	switch( codes[base++] ) {
		case  1: xd_d[d++] = xd_s[0]; break;
		default: break;
	}
	
	// class
	switch( codes[base++] ) {
		case  1: xd_d[d++] = xd_s[1]; break;
		default: break;
	}
	
	/**
	printf("map: ");
	printprodvec( K_s , xc_s , B_s , Bo_s , xb_s , D_s , xd_s );
	printf(" -> ");
	printprodvec( K_d , xc_d , B_d , Bo_d[0] , xb_d , D_d , xd_d );
	printf("\n");
	**/
	 
	return 1;
}

int computemetrics(MNL_MODEL * model , 
				   MNL_DATA * data , 
				   double * metrics , 
				   double ** minae , 
				   double ** maxae , 
				   double ** bestp , 
				   double ** beste , 
				   double * bestm )
{
	int m, j, jp, info, tmpj;
	double ** PL, ** AE, ** RE , ** tmpE, tmpmin; 
	double loglik;
	
	if( model == NULL || data == NULL || metrics == NULL ) { return 0; }
	
	// allocate space for Logit choice probabilities
	PL = (double **)calloc( data->M , sizeof(double *) );
	for( m = 0 ; m < data->M ; m++ ) {
		PL[m] = (double *)calloc( (data->Jm)[m] , sizeof(double) );
	}
	
	// compute Logit choice probabilities
	info = mnl_model_probabilities( model , data , PL , &loglik );
	
	// act based on outcome
	switch( info ) {
			
		case 1: // success in computing choice probabilities; compute metric(s)
			
			// allocate memory for errors
			AE = (double **)calloc( data->M , sizeof(double *) );
			RE = (double **)calloc( data->M , sizeof(double *) );
			for( m = 0 ; m < data->M ; m++ ) {
				AE[m] = (double *)calloc( (data->Jm)[m] , sizeof(double) );
				RE[m] = (double *)calloc( (data->Jm)[m] , sizeof(double) );
			}
			
			// compute share errors
			metrics[0] = 0.0;
			metrics[1] = 0.0;
			metrics[2] = 0.0;
			for( m = 0 ; m < data->M ; m++ ) {
				for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
					
					// absolute share error
					(AE[m])[j] = ABS( ((data->s)[m])[j] - (PL[m])[j] );
					
					// relative share error
					(RE[m])[j] = (AE[m])[j] / ((data->s)[m])[j];
					
					// first metric: KLD divergence (one of several likelihood metrics)
					// metrics[0] += ((data->s)[m])[j] * log( ((data->s)[m])[j] / (PL[m])[j] );
					metrics[0] += ((data->s)[m])[j] * log( ((data->s)[m])[j] );
					
					// second metric: average absolute share error
					metrics[1] += (AE[m])[j];
					
					// third metric: weighted average absolute share error
					metrics[2] += ((data->s)[m])[j] * (AE[m])[j];
					
				}
			}
			metrics[0] -= loglik;
			metrics[1] /= ( 1.0 * (data->J) ); // finish average absolute share error
			metrics[2] /= ( 1.0 * (data->M) ); // finish weighted average absolute share error
			
			if( bestp != NULL && beste != NULL && bestm != NULL ) { 
				
				if( metrics[0] < bestm[0] ) { // replace best metric value and copy best shares
					bestm[0] = metrics[0];
					for( m = 0 ; m < data->M ; m++ ) {
						if( bestp[m] != NULL ) {
							for( j = 0 ; j < (data->Jm)[m] ; j++ ) { 
								(bestp[m])[j] = (PL[m])[j]; 
							}
						}
						if( beste[m] != NULL ) {
							for( j = 0 ; j < (data->Jm)[m] ; j++ ) { 
								(beste[m])[j] = (AE[m])[j]; 
							}
						}
					}
				}
				
			}
			
			// min/max envelopes
			if( minae != NULL || maxae != NULL ) {
				
				tmpE = (double **)calloc( data->M , sizeof(double *) );
				for( m = 0 ; m < data->M ; m++ ) {
					tmpE[m] = (double *)calloc( (data->Jm)[m] , sizeof(double) );
				}
				
				// sort errors to be in ascending order (inefficient method; quicksort is better)
				for( m = 0 ; m < data->M ; m++ ) {
					for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
						tmpmin = 1.0e20; tmpj = 0; // initialize
						for( jp = 0 ; jp < (data->Jm)[m] ; jp++ ) {
							if( (AE[m])[jp] < tmpmin ) { tmpmin = (AE[m])[jp]; tmpj = jp; }
						}
						(tmpE[m])[j]  = tmpmin; // store the next largest entry in temp storage
						(AE[m])[tmpj] = 1.0e20; // make the element just stored large enough to ignore
					}
				}
				for( m = 0 ; m < data->M ; m++ ) { // replace
					for( j = 0 ; j < (data->Jm)[m] ; j++ ) { (AE[m])[j] = (tmpE[m])[j]; }
				}
				
				// update min/max envelopes
				if( minae != NULL ) { 
					for( m = 0 ; m < data->M ; m++ ) {
						if( minae[m] != NULL ) {
							for( j = 0 ; j < (data->Jm)[m] ; j++ ) { 
								(minae[m])[j] = MIN( (minae[m])[j] , (AE[m])[j] );
							}
						}
					}
				}
				if( maxae != NULL ) { 
					for( m = 0 ; m < data->M ; m++ ) {
						if( maxae[m] != NULL ) {
							for( j = 0 ; j < (data->Jm)[m] ; j++ ) { 
								(maxae[m])[j] = MAX( (maxae[m])[j] , (AE[m])[j] );
							}
						}
					}
				}
				
			}
			
			// release memory allocated for errors
			for( m = 0 ; m < data->M ; m++ ) { free(RE[m]); } free(RE);
			for( m = 0 ; m < data->M ; m++ ) { free(AE[m]); } free(AE);
			
			break;
			
		default: break;
			
	}
	
	// free Logit choice probability memory 
	for( m = 0 ; m < data->M ; m++ ) { free(PL[m]); } free(PL);
	
	return info;
	
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int run_mle(char * efpn , 
			char * pfpn , 
			int T , 
			int maxloops , 
			int checkders , 
			int scale , 
			char * lfpn , 
			char * rfpn , 
			char * cfpn ,
			char * bfpn )
{
	int i, n, m, j, k, b, d, l;
	
	int info, rank, loopcounter;
	
	FILE * est_dfp , * pred_dfp , * rfp , * cfp;
	
	MNL_DATA  * edata, * pdata , * emodeldata , * pmodeldata;
	MNL_MODEL * model , * pmodel;
	
	double ** minae , ** maxae , ** bestp , ** beste , bestm;
	
	double norm, loglik, OptTol, FeasTol;
	
	double rates[4] = { 0.0 };
	
	double metrics[3] = { 0 };
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	est_dfp  = fopen( efpn , "r" ); // open estimation data file
	pred_dfp = fopen( pfpn , "r" ); // open prediction data file
	if( est_dfp == NULL || pred_dfp == NULL ) { return -1; }
	
	// allocate and instantiate new mnl_data structures
	edata = (MNL_DATA *)malloc( sizeof(MNL_DATA) ); mnl_data_new( edata );
	pdata = (MNL_DATA *)malloc( sizeof(MNL_DATA) ); mnl_data_new( pdata );
	
	// attempt to read estimation data from csv file
	switch( mnl_data_read( est_dfp , 3 , 'y' , 0 , NULL , edata ) ) {
		case -3: printf("estimation data read error (-3)\n"); break;
		case -2: printf("estimation data read error (-2)\n"); break;
		case -1: printf("estimation data read error (-1)\n"); break;
		default: printf("estimation data read correctly.\n"); break;
	}
	fclose(est_dfp);
	
	// mnl_data_print( edata );
	
	// attempt to read prediction data from csv file
	switch( mnl_data_read( pred_dfp , 3 , 'y' , 0 , NULL , pdata ) ) {
		case -3: printf("prediction data read error (-3)\n"); break;
		case -2: printf("prediction data read error (-2)\n"); break;
		case -1: printf("prediction data read error (-1)\n"); break;
		default: printf("prediction data read correctly.\n"); break;
	}
	fclose(pred_dfp);
	
	// we can do this here, and need to do this here, because markets and products
	// are copied when mapping data
	minae = (double **)calloc( pdata->M , sizeof(double *) );
	maxae = (double **)calloc( pdata->M , sizeof(double *) );
	bestp = (double **)calloc( pdata->M , sizeof(double *) );
	beste = (double **)calloc( pdata->M , sizeof(double *) );
	for( m = 0 ; m < pdata->M ; m++ ) {
		minae[m] = (double *)calloc( (pdata->Jm)[m] , sizeof(double) );
		maxae[m] = (double *)calloc( (pdata->Jm)[m] , sizeof(double) );
		bestp[m] = (double *)calloc( (pdata->Jm)[m] , sizeof(double) );
		beste[m] = (double *)calloc( (pdata->Jm)[m] , sizeof(double) );
		for( j = 0 ; j < (pdata->Jm)[m] ; j++ ) { (minae[m])[j] = 1.0e20; }
	}
	
	// mnl_data_print( pdata );
	
	// allocate and instantiate an empty model structure for estimation
	model  = (MNL_MODEL*)malloc( sizeof(MNL_MODEL) );
	pmodel = (MNL_MODEL*)malloc( sizeof(MNL_MODEL) );
	
	// allocated estimation and prediction data set structures (modified from data via mapping)
	emodeldata = ( MNL_DATA * )calloc( 1 , sizeof( MNL_DATA ) ); 
	pmodeldata = ( MNL_DATA * )calloc( 1 , sizeof( MNL_DATA ) ); 
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * RUN * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	rfp = fopen( rfpn , "w" );
	if( rfp == NULL ) { printf("error initializing results file. \n"); exit(1); }
	
	cfp = fopen( cfpn , "w" );
	if( cfp == NULL ) { printf("error initializing coefficients file. \n"); exit(1); }
	
	loopcounter = 0;
	
	bestm = 1.0e20; // initialized high enough to guarantee overwrite
	
	do { // while( getnext( codeN , codeL , codes ) ); see below
		
		// open results file for this code
		rfp = fopen( rfpn , "a" );
		if( rfp == NULL ) { printf("error reopening results file. \n"); exit(1); }
		
		cfp = fopen( cfpn , "a" );
		if( cfp == NULL ) { printf("error reopening coefficients file. \n"); exit(1); }
		
		// print out current codes to both results and coefficients files
		
		fprintf( rfp , "%i" , codes[0] );
		for( i = 1 ; i < 10 ; i++ ) { fprintf( rfp , "%i" , codes[i] ); }
		
		fprintf( cfp , "%i" , codes[0] );
		for( i = 1 ; i < 10 ; i++ ) { fprintf( cfp , "%i" , codes[i] ); }
		
		// modify data for model estimation
		
		// re-inialize estimation data structure (make sure this is released via mnl_model_free below)
		mnl_data_new( emodeldata );
		
		// define sizes based on the current codes
		definesizes( emodeldata , codes );
		
		// map all data into estimation data
		info = mnl_data_map( edata , emodeldata , (MNL_DATA_MAP)map , codes );
		
		// mnl_data_print( emodeldata );
		
		// print out map code
		fprintf( rfp , ", %i " , info ); 
		fprintf( cfp , ", %i " , info ); 
		
		switch( info ) {
				
			default: break;
				
			case 1: // mapping successful
				
				// print out mapped sizes
				fprintf( rfp , ", %i , %i , %i , %i " , emodeldata->K , emodeldata->B , emodeldata->D , emodeldata->L ); 
				fprintf( cfp , ", %i , %i , %i , %i " , emodeldata->K , emodeldata->B , emodeldata->D , emodeldata->L ); 
				
				// check model
				switch( mnl_data_check( emodeldata , &rank , &norm ) ) {
						
					case  1: 
						printf("check succeeded: ");
						printf("rank = %i (%i x %i), ",rank,emodeldata->D+emodeldata->J,emodeldata->K+emodeldata->B+emodeldata->L);
						printf("norm = %0.4e\n",norm); 
						break;
						
					default: 
						rank = -1; 
						norm = -1.0;
						printf("check failed\n"); 
						exit(1);
						break;
						
				}
				
				// print out rank and norm
				fprintf( rfp , ", %i , %0.16f " , rank , norm ); 
				fprintf( cfp , ", %i , %0.16f " , rank , norm ); 
				
				// exit(1);
				
				// initialize model structure (and free below)
				mnl_model_new( model );
				
				// attempt to solve estimation problem (store info to use in file print)
				info = mnl_mle_snopt(emodeldata , 
									 model , 
									 T , 
									 &loglik , 
									 rates , 
									 scale , 
									 checkders , 
									 'r' ,			// random initial conditions
									 &OptTol , 
									 &FeasTol , 
									 lfpn );
				
				// print out info and estimation success/failure rates to results file
				fprintf( rfp , ", %i " , info ); 
				fprintf( rfp , ", %0.4e " , OptTol  );
				fprintf( rfp , ", %0.4e " , FeasTol );
				fprintf( rfp , ", %0.2f " , 100.0 * rates[0] );
				fprintf( rfp , ", %0.2f " , 100.0 * rates[1] );
				fprintf( rfp , ", %0.2f " , 100.0 * rates[2] );
				fprintf( rfp , ", %0.2f " , 100.0 * rates[3] ); 
				
				// print out info and estimation success/failure rates to coefficients file
				fprintf( cfp , ", %i " , info ); 
				fprintf( cfp , ", %0.4e " , OptTol  );
				fprintf( cfp , ", %0.4e " , FeasTol );
				fprintf( cfp , ", %0.2f " , 100.0 * rates[0] );
				fprintf( cfp , ", %0.2f " , 100.0 * rates[1] );
				fprintf( cfp , ", %0.2f " , 100.0 * rates[2] );
				fprintf( cfp , ", %0.2f " , 100.0 * rates[3] ); 
				
				switch( info ) {
						
					default: break;
					case 1: 
						
						// print out (negative) log likelihood to results file
						fprintf( rfp , ", %0.16f " , -loglik ); 
						
						// print computed model coefficients to coefficients file
						fprintf( cfp , ", " ); mnl_model_fprint( model , cfp );
						
						// modify data for model prediction
						
						// re-inialize prediction data structure (make sure freed below)
						mnl_data_new( pmodeldata );
						
						// define sizes based on this code
						definesizes( pmodeldata , codes );
						
						// map all data into prediction data
						info = mnl_data_map( pdata , pmodeldata , (MNL_DATA_MAP)map , codes );
						fprintf( rfp , ", %i "  , info ); 
						switch( info ) {
								
							default: break;
								
							case 1:
								
								info = computemetrics( model , pmodeldata , metrics , minae , maxae , bestp , beste , &bestm );
								fprintf( rfp , ", %i "  , info ); 
								
								switch( info ) {
										
									default: break;
										
									case 1: // success
										
										// print out metrics to results file
										fprintf( rfp , ", %0.16f " , metrics[0] ); 
										fprintf( rfp , ", %0.16f " , metrics[1] ); 
										fprintf( rfp , ", %0.16f " , metrics[2] ); 
										
										// check for best model fit
										
										
										// estimate model on * prediction * data and compute metrics to describe overfitting
										
										// we don't free model so that we can use it as an initial condition
										
										// attempt to solve prediction data estimation problem
										info = mnl_mle_snopt(pmodeldata , 
															 model , 
															 T , 
															 &loglik , 
															 rates , 
															 scale , 
															 checkders , 
															 'm' , 
															 &OptTol , 
															 &FeasTol , 
															 lfpn );
										fprintf( rfp , ", %i "  , info ); 
										fprintf( rfp , ", %0.4e " , OptTol  );
										fprintf( rfp , ", %0.4e " , FeasTol );
										switch( info ) {
											default: break;
											case 1: 
												// print out (negative) log likelihood to results file
												fprintf( rfp , ", %0.16f " , -loglik );
												// compute metrics for prediction model fit on prediction data
												info = computemetrics( model , pmodeldata , metrics , NULL , NULL , NULL , NULL , NULL );
												fprintf( rfp , ", %i "  , info ); 
												switch( info ) {
													case 1: // success: print out metrics to results file
														fprintf( rfp , ", %0.16f " , metrics[0] ); 
														fprintf( rfp , ", %0.16f " , metrics[1] );  
														fprintf( rfp , ", %0.16f " , metrics[2] ); 
														break;
													default: break;
												}
												break;
										}
										
										break;
										
								}
								
								break;
								
						} // end prediction mapping switch
						
						
						break;
						
				} // end estimation switch
				
				// free model structure so that it can be used again
				mnl_model_free( model );
				
				break;
				
		} // end estimation data map switch
		
		// free estimation data structure for next call
		mnl_data_free( emodeldata );
		
		// end line and close results file
		fprintf( rfp , "\n" ); fclose( rfp );
		
		// end line and close coefficients file
		fprintf( cfp , "\n" ); fclose( cfp );
		
		// 
		loopcounter++;
		if( loopcounter >= maxloops ) { break; }
		
	} while( getnext( codeN , codeL , codeU , codes ) );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * CREATE BEST PREDICTION FILE * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	rfp = fopen( bfpn , "w" );
	if( rfp != NULL ) { 
		
		for( m = 0 ; m < pdata->M ; m++ ) {
			for( j = 0 ; j < (pdata->Jm)[m] ; j++ ) { 
				fprintf( rfp , "%i , %i "  , m+1 , j+1 );
				fprintf( rfp , ", %0.16f " , ((pdata->s)[m])[j] );
				fprintf( rfp , ", %0.16f " , (bestp[m])[j] );
				fprintf( rfp , ", %0.16f " , (beste[m])[j] );
				fprintf( rfp , ", %0.16f " , (minae[m])[j] );
				fprintf( rfp , ", %0.16f\n", (maxae[m])[j] ); 
			}
		}
		
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * CLEAN UP  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// free data storage
	
	for( m = 0 ; m < pdata->M ; m++ ) { 
		free(minae[m]); free(maxae[m]); free(bestp[m]); free(beste[m]); 
	}
	free(minae); free(maxae); free(bestp); free(beste);
	
	mnl_data_free( edata ); free(edata); 
	mnl_data_free( pdata ); free(pdata);
	
	free(model); free(pmodel);
	
	return 1;
	
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int run_gmm(char * efpn , 
			char * pfpn , 
			int T , 
			int maxloops , 
			int checkders , 
			int scale , 
			char * lfpn , 
			char * rfpn , 
			char * cfpn ,
			char * bfpn )
{
	int i, n, m, j, k, b, d, l;
	
	int info, rank, loopcounter;
	
	FILE * est_dfp , * pred_dfp , * rfp , * cfp;
	
	MNL_DATA  * edata, * pdata , * emodeldata , * pmodeldata;
	MNL_MODEL * model , * pmodel;
	
	double ** minae , ** maxae , ** bestp , ** beste , bestm;
	
	double norm, gmmres, OptTol, FeasTol;
	
	double rates[4] = { 0.0 };
	
	double metrics[3] = { 0 };
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	est_dfp  = fopen( efpn , "r" ); // open estimation data file
	pred_dfp = fopen( pfpn , "r" ); // open prediction data file
	if( est_dfp == NULL || pred_dfp == NULL ) { return -1; }
	
	// allocate and instantiate new mnl_data structures
	edata = (MNL_DATA *)malloc( sizeof(MNL_DATA) ); mnl_data_new( edata );
	pdata = (MNL_DATA *)malloc( sizeof(MNL_DATA) ); mnl_data_new( pdata );
	
	// attempt to read estimation data from csv file
	switch( mnl_data_read( est_dfp , 3 , 'y' , 0 , NULL , edata ) ) {
		case -3: printf("estimation data read error (-3)\n"); break;
		case -2: printf("estimation data read error (-2)\n"); break;
		case -1: printf("estimation data read error (-1)\n"); break;
		default: printf("estimation data read correctly.\n"); break;
	}
	fclose(est_dfp);
	
	// mnl_data_print( edata );
	
	// attempt to read prediction data from csv file
	switch( mnl_data_read( pred_dfp , 3 , 'y' , 0 , NULL , pdata ) ) {
		case -3: printf("prediction data read error (-3)\n"); break;
		case -2: printf("prediction data read error (-2)\n"); break;
		case -1: printf("prediction data read error (-1)\n"); break;
		default: printf("prediction data read correctly.\n"); break;
	}
	fclose(pred_dfp);
	
	// we can do this here, and need to do this here, because markets and products
	// are copied when mapping data
	minae = (double **)calloc( pdata->M , sizeof(double *) );
	maxae = (double **)calloc( pdata->M , sizeof(double *) );
	bestp = (double **)calloc( pdata->M , sizeof(double *) );
	beste = (double **)calloc( pdata->M , sizeof(double *) );
	for( m = 0 ; m < pdata->M ; m++ ) {
		minae[m] = (double *)calloc( (pdata->Jm)[m] , sizeof(double) );
		maxae[m] = (double *)calloc( (pdata->Jm)[m] , sizeof(double) );
		bestp[m] = (double *)calloc( (pdata->Jm)[m] , sizeof(double) );
		beste[m] = (double *)calloc( (pdata->Jm)[m] , sizeof(double) );
		for( j = 0 ; j < (pdata->Jm)[m] ; j++ ) { (minae[m])[j] = 1.0e20; }
	}
	
	// mnl_data_print( pdata );
	
	// allocate and instantiate an empty model structure for estimation
	model = (MNL_MODEL*)malloc( sizeof(MNL_MODEL) );
	pmodel = (MNL_MODEL*)malloc( sizeof(MNL_MODEL) );
	
	// allocated estimation and prediction data set structures (modified from data via mapping)
	emodeldata = ( MNL_DATA * )calloc( 1 , sizeof( MNL_DATA ) ); 
	pmodeldata = ( MNL_DATA * )calloc( 1 , sizeof( MNL_DATA ) ); 
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * RUN * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	rfp = fopen( rfpn , "w" );
	if( rfp == NULL ) { printf("error initializing results file. \n"); exit(1); }
	
	cfp = fopen( cfpn , "w" );
	if( cfp == NULL ) { printf("error initializing coefficients file. \n"); exit(1); }
	
	loopcounter = 0;
	
	bestm = 1.0e20; // initialized high enough to guarantee overwrite
	
	do { // while( getnext( codeN , codeL , codes ) ); see below
		
		// open results file for this code
		rfp = fopen( rfpn , "a" );
		if( rfp == NULL ) { printf("error reopening results file. \n"); exit(1); }
		
		cfp = fopen( cfpn , "a" );
		if( cfp == NULL ) { printf("error reopening coefficients file. \n"); exit(1); }
		
		// print out current codes to both results and coefficients files
		
		fprintf( rfp , "%i" , codes[0] );
		for( i = 1 ; i < 10 ; i++ ) { fprintf( rfp , "%i" , codes[i] ); }
		
		fprintf( cfp , "%i" , codes[0] );
		for( i = 1 ; i < 10 ; i++ ) { fprintf( cfp , "%i" , codes[i] ); }
		
		// modify data for model estimation
		
		// re-inialize estimation data structure (make sure this is released via mnl_model_free below)
		mnl_data_new( emodeldata );
		
		// define sizes based on the current codes
		definesizes( emodeldata , codes );
		
		// map all data into estimation data
		info = mnl_data_map( edata , emodeldata , (MNL_DATA_MAP)map , codes );
		
		// mnl_data_print( emodeldata );
		
		// print out map code
		fprintf( rfp , ", %i " , info ); 
		fprintf( cfp , ", %i " , info ); 
		
		switch( info ) {
				
			default: break;
				
			case 1: // mapping successful
				
				// print out mapped sizes
				fprintf( rfp , ", %i , %i , %i , %i " , emodeldata->K , emodeldata->B , emodeldata->D , emodeldata->L ); 
				fprintf( cfp , ", %i , %i , %i , %i " , emodeldata->K , emodeldata->B , emodeldata->D , emodeldata->L ); 
				
				// check model
				switch( mnl_data_check( emodeldata , &rank , &norm ) ) {
						
					case  1: 
						printf("check succeeded: ");
						printf("rank = %i (%i x %i), ",rank,emodeldata->D+emodeldata->J,emodeldata->K+emodeldata->B+emodeldata->L);
						printf("norm = %0.4e\n",norm); 
						break;
						
					default: 
						rank = -1; 
						norm = -1.0;
						printf("check failed\n"); 
						exit(1);
						break;
						
				}
				
				// print out rank and norm
				fprintf( rfp , ", %i , %0.16f " , rank , norm ); 
				fprintf( cfp , ", %i , %0.16f " , rank , norm ); 
				
				// exit(1);
				
				// initialize model structure (and free below)
				mnl_model_new( model );
				
				// attempt to solve estimation problem (store info to use in file print)
				info = mnl_gmm_snopt(emodeldata , 
									 model , 
									 T , 
									 &gmmres , 
									 rates , 
									 scale , 
									 checkders , 
									 'r' ,			// random initial conditions
									 &OptTol , 
									 &FeasTol , 
									 lfpn );
				
				// print out info and estimation success/failure rates to results file
				fprintf( rfp , ", %i " , info ); 
				fprintf( rfp , ", %0.4e " , OptTol  );
				fprintf( rfp , ", %0.4e " , FeasTol );
				fprintf( rfp , ", %0.2f " , 100.0 * rates[0] );
				fprintf( rfp , ", %0.2f " , 100.0 * rates[1] );
				fprintf( rfp , ", %0.2f " , 100.0 * rates[2] );
				fprintf( rfp , ", %0.2f " , 100.0 * rates[3] ); 
				
				// print out info and estimation success/failure rates to coefficients file
				fprintf( cfp , ", %i " , info ); 
				fprintf( cfp , ", %0.4e " , OptTol  );
				fprintf( cfp , ", %0.4e " , FeasTol );
				fprintf( cfp , ", %0.2f " , 100.0 * rates[0] );
				fprintf( cfp , ", %0.2f " , 100.0 * rates[1] );
				fprintf( cfp , ", %0.2f " , 100.0 * rates[2] );
				fprintf( cfp , ", %0.2f " , 100.0 * rates[3] ); 
				
				switch( info ) {
						
					default: break;
					case 1: 
						
						// print out (negative) log likelihood to results file
						fprintf( rfp , ", %0.16f " , gmmres ); 
						
						// print computed model coefficients to coefficients file
						fprintf( cfp , ", " ); mnl_model_fprint( model , cfp );
						
						// modify data for model prediction
						
						// re-inialize prediction data structure (make sure freed below)
						mnl_data_new( pmodeldata );
						
						// define sizes based on this code
						definesizes( pmodeldata , codes );
						
						// map all data into prediction data
						info = mnl_data_map( pdata , pmodeldata , (MNL_DATA_MAP)map , codes );
						fprintf( rfp , ", %i "  , info ); 
						switch( info ) {
								
							default: break;
								
							case 1:
								
								info = computemetrics( model , pmodeldata , metrics , minae , maxae , bestp , beste , &bestm );
								fprintf( rfp , ", %i "  , info ); 
								
								switch( info ) {
										
									default: break;
										
									case 1: // success
										
										// print out metrics to results file
										fprintf( rfp , ", %0.16f " , metrics[0] ); 
										fprintf( rfp , ", %0.16f " , metrics[1] ); 
										fprintf( rfp , ", %0.16f " , metrics[2] ); 
										
										// check for best model fit
										
										
										// estimate model on * prediction * data and compute metrics to describe overfitting
										
										// we don't free model so that we can use it as an initial condition
										
										// attempt to solve prediction data estimation problem
										info = mnl_gmm_snopt(pmodeldata , 
															 model , 
															 T , 
															 &gmmres , 
															 rates , 
															 scale , 
															 checkders , 
															 'm' , 
															 &OptTol , 
															 &FeasTol , 
															 lfpn );
										fprintf( rfp , ", %i "  , info ); 
										fprintf( rfp , ", %0.4e " , OptTol  );
										fprintf( rfp , ", %0.4e " , FeasTol );
										switch( info ) {
											default: break;
											case 1: 
												// print out (negative) log likelihood to results file
												fprintf( rfp , ", %0.16f " , gmmres );
												// compute metrics for prediction model fit on prediction data
												info = computemetrics( model , pmodeldata , metrics , NULL , NULL , NULL , NULL , NULL );
												fprintf( rfp , ", %i "  , info ); 
												switch( info ) {
													case 1: // success: print out metrics to results file
														fprintf( rfp , ", %0.16f " , metrics[0] ); 
														fprintf( rfp , ", %0.16f " , metrics[1] );  
														fprintf( rfp , ", %0.16f " , metrics[2] ); 
														break;
													default: break;
												}
												break;
										}
										
										break;
										
								}
								
								break;
								
						} // end prediction mapping switch
						
						
						break;
						
				} // end estimation switch
				
				// free model structure so that it can be used again
				mnl_model_free( model );
				
				break;
				
		} // end estimation data map switch
		
		// free estimation data structure for next call
		mnl_data_free( emodeldata );
		
		// end line and close results file
		fprintf( rfp , "\n" ); fclose( rfp );
		
		// end line and close coefficients file
		fprintf( cfp , "\n" ); fclose( cfp );
		
		// 
		loopcounter++;
		if( loopcounter >= maxloops ) { break; }
		
	} while( getnext( codeN , codeL , codeU , codes ) );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * CREATE BEST PREDICTION FILE * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	rfp = fopen( bfpn , "w" );
	if( rfp != NULL ) { 
		
		for( m = 0 ; m < pdata->M ; m++ ) {
			for( j = 0 ; j < (pdata->Jm)[m] ; j++ ) { 
				fprintf( rfp , "%i , %i "  , m+1 , j+1 );
				fprintf( rfp , ", %0.16f " , ((pdata->s)[m])[j] );
				fprintf( rfp , ", %0.16f " , (bestp[m])[j] );
				fprintf( rfp , ", %0.16f " , (beste[m])[j] );
				fprintf( rfp , ", %0.16f " , (minae[m])[j] );
				fprintf( rfp , ", %0.16f\n", (maxae[m])[j] ); 
			}
		}
		
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * CLEAN UP  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// free data storage
	
	for( m = 0 ; m < pdata->M ; m++ ) { 
		free(minae[m]); free(maxae[m]); free(bestp[m]); free(beste[m]); 
	}
	free(minae); free(maxae); free(bestp); free(beste);
	
	mnl_data_free( edata ); free(edata); 
	mnl_data_free( pdata ); free(pdata);
	
	free(model); free(pmodel);
	
	return 1;
	
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int main( int argc , char * argv[] )
{
	int i, n;
	
	int T, scale, checkders, maxloops, method;
	
	char runnm[256], lfpn[256], rfpn[256], cfpn[256], bfpn[256], efpn[256], pfpn[256];
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * DEFAULTS  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	T = 1;						// number of solve attempts per code
	maxloops  = BIGGEST_INT;	// maximum number of loops
	scale = 0;					// do not scale
	checkders = 0;				// do not check derivatives
	method = 'l';				// use MLE
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	printf("\nReading command-line options...\n");
	
	// read arguments
	for( i = 1 ; i < argc ; i++ ) {
		
		// parse next command line option (only if starting with the '-' character)
		if( (argv[i])[0] == '-' ) {
			
			switch( (argv[i])[1] ) {
					
				default: 
					printf("  Option%c%c%c not recognized.\n",'"',(argv[i])[1],'"');
					break;
					
				case 'T': // reading number of trials
					T = (int)strtol( argv[i] + 2 , NULL , 10 );
					if( T < 1 ) {
						printf("  Invalid number of trials (%i): T must be > 0.\n",T);
						T = 1;
					}
					break;
					
				case 'C': // reading code
					for( n = 9 ; n >= 0 ; n-- ) {
						codes[n] = MIN( codeU[n] , MAX( codeL[n] , (int)strtol( argv[i] + 2 + n , NULL , 10 ) ) );
						(argv[i])[2+n] = '\0';
					}
					printf("  Read the following code: %i",codes[0]);
					for( n = 1 ; n < 10 ; n++ ) { printf(", %i",codes[n]); }
					printf("\n");
					break;
					
				case 'M': // reading method
					switch( (argv[i])[2] ) {
						case 'm':
						case 'l': method = (argv[i])[2]; break;
						default:  method = 'l'; break;
					}
					break;
					
				case 'L': // reading number of loops allowed
					maxloops = (int)strtol( argv[i] + 2 , NULL , 10 );
					if( maxloops < 1 ) {
						printf("  Invalid number of max loops (%i): must be > 0.\n",maxloops);
						maxloops = BIGGEST_INT;
					}
					break;
					
				case 'S': 
					printf("  Scaling characteristics\n");
					scale = 1; 
					break;
					
				case 'D': 
					printf("  Checking derivatives\n");
					checkders = 1; 
					break;
					
					
			}
			
		}
		
	}
	
	printf("\n");
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * SETUP * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	sprintf(runnm,"MNL_SNOPT");
	switch( method ) {
		case 'l': strcat(runnm,"_MLE"); break;
		case 'm': strcat(runnm,"_GMM"); break;
		default: break;
	}
	
	strcpy(lfpn,runnm);
	strcat(lfpn,"_log.log");
	
	strcpy(rfpn,runnm);
	strcat(rfpn,"_results.log");
	
	strcpy(cfpn,runnm);
	strcat(cfpn,"_coeffs.log");
	
	strcpy(bfpn,runnm);
	strcat(bfpn,"_best.log");
	
	strcpy(efpn,"GraceCase_2005-2006.csv");
	strcpy(pfpn,"GraceCase_2007.csv");
	
	switch( method ) {
		case 'l': run_mle(efpn,pfpn,T,maxloops,checkders,scale,lfpn,rfpn,cfpn,bfpn); break;
		case 'm': run_gmm(efpn,pfpn,T,maxloops,checkders,scale,lfpn,rfpn,cfpn,bfpn); break;
		default: break;
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	return 0;
	
}