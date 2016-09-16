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

#include <vecLib/cblas.h>

#include "mnl_macros.h"
#include "mnl_data.h"
#include "rcl_data.h"

#define BUFFER_SIZE 1024

static int rngon = 0;
static UNUR_DISTR * distr;
static UNUR_PAR   * param;
static UNUR_URNG  * unrng_un;
static UNUR_GEN   * unrng;
static UNUR_URNG  * unrng_sn;
static UNUR_GEN   * snrng;
static unsigned long seed[] = {111u, 222u, 333u, 444u, 555u, 666u};

int rcl_imn_sampler(int I , RCL_MODEL * model , 
					double * cc , double * bc , double * dc , 
					void * params )
{
	int i, k, b, d, l;
	double * VC , * VB , * VD;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	if( model == NULL ) { return -1; }
	if( model->K > 0 && cc == NULL ) { return -1; }
	if( model->B > 0 && bc == NULL ) { return -1; }
	if( model->D > 0 && dc == NULL ) { return -1; }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// allocate sample arrays
	VC = (double *)calloc( I*(model->K) , sizeof(double) ); // continuous coefficients
	VB = (double *)calloc( I*(model->B) , sizeof(double) ); // binary coefficients
	VD = (double *)calloc( I*(model->L) , sizeof(double) ); // dummy coefficients
	
	// make sure sampler is started
	if( rngon ) {
		
		// define package seed
		RngStream_SetPackageSeed( seed );
		
		// Make a uniform generator stream object for standard normal random number generator.
		unrng_sn = unur_urng_rngstream_new( "unrng_sn" );
		if( unrng_sn == NULL ) { 
			printf("ERROR - Uniform generator could not be constructed.\n");
			exit(EXIT_FAILURE); 
		}
		
		// use predefined distribution - standard normal
		distr = unur_distr_normal( NULL , 0 );
		
		// use "TDR" method (why, not sure)
		param = unur_tdr_new( distr );
		
		// Set uniform generator in parameter object
		unur_set_urng( param , snrng );
		
		// Create the uniform random number generator object.
		snrngen = unur_init( param ); // param is "destroyed" here
		if( snrngen == NULL ) { 
			printf("ERROR - standard normal generator could not be constructed.\n");
		}
		
		// free distribution object
		unur_distr_free( distr );
		
	}
	
	// sample 
	for( k = 0 ; k < data->K ; k++ ) { 
		for( i = 0 ; i < I ; i++ ) {
			VC[I*k+i] = unur_sample_cont(snrng);
		}
	}
	
	for( b = 0 ; b < data->B ; b++ ) { 
		for( i = 0 ; i < I ; i++ ) {
			VB[I*b+i] = unur_sample_cont(snrng);
		}
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// compute coefficients
	
	if( data->K > 0 ) {
		for( i = 0 ; i < I ; i++ ) {
			cblas_dcopy( data->K , model->cc_m , 1 , cc + (data->K)*i , 1 ); // cc(i,:) <- cc_m
		}
		for( k = 0 ; k < data->K ; k++ ) { // cc(:,k) <- cc(:,k) + cc_v(k) VC(:,k)
			cblas_daxpy( I , (model->cc_v)[k] , VC + I*k , 1 , cc + k , data->K );
		}
	}
	
	if( data->B > 0 ) {
		for( i = 0 ; i < I ; i++ ) {
			cblas_dcopy( data->B , model->bc_m , 1 , bc + (data->B)*i , 1 ); // bc(i,:) <- bc_m
		}
		for( b = 0 ; b < data->B ; b++ ) { // bc(:,b) <- bc(:,k) + bc_v(k) VB(:,b)
			cblas_daxpy( I , (model->bc_v)[b] , VB + I*b , 1 , bc + b , data->B );
		}
	}
	
	if( data->D > 0 ) {
		
		
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	free(VC); 
	free(VB); 
	free(VD);
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
}

void rcl_model_new( MNL_MODEL * model )
{
	if( model == NULL ) { return; }
	
	model->on = 'n';
	model->og = 'n';
	
	model->K = 0;
	
	model->B = 0;
	
	model->D = 0;
	model->Ld = NULL;
	model->L = 0;
	
	model->cc_m = NULL;
	model->bc_m = NULL;
	model->dc_m = NULL;
	model->uc_m = 0.0;
	
	model->cc_v = NULL;
	model->bc_v = NULL;
	model->dc_v = NULL;
	model->uc_v = 0.0;
	
}

void rcl_model_alloc( MNL_MODEL * model , int K , int B , int D , int * Ld )
{
	int d;
	
	if( model == NULL ) { return; }
	if( K == 0 && B == 0 && D == 0 ) { return; }
	if( D > 0 && Ld == NULL ) { return; }
	
	model->K = K;
	model->cc_m = ( double * )calloc( K , sizeof( double ) );
	model->cc_v = ( double * )calloc( K , sizeof( double ) );
	
	model->B = B;
	model->bc_m = ( double * )calloc( B , sizeof(double) );
	model->bc_v = ( double * )calloc( B , sizeof(double) );
	
	model->D = D; 
	model->L = 0;
	model->Ld = ( int * )calloc( D , sizeof(int) );
	for( d = 0 ; d < D ; d++ ) { (model->Ld)[d] = Ld[d]; model->L += Ld[d]; }
	model->dc_m = ( double * )calloc( model->L , sizeof(double) );
	model->dc_v = ( double * )calloc( model->L , sizeof(double) );
	
	model->on = 'y';
}

void rcl_model_free( MNL_MODEL * model )
{
	if( model == NULL ) { return; }
	
	if( model->Ld   != NULL ) { free(model->Ld);   }
	if( model->cc_m != NULL ) { free(model->cc_m); }
	if( model->cc_v != NULL ) { free(model->cc_v); }
	if( model->bc_m != NULL ) { free(model->bc_m); }
	if( model->bc_v != NULL ) { free(model->bc_v); }
	if( model->dc_m != NULL ) { free(model->dc_m); }
	if( model->dc_v != NULL ) { free(model->dc_v); }
	
	model->K = 0;
	model->B = 0;
	model->D = 0;
	model->L = 0;
	model->og = 'n';
	model->uc = 0.0;
	
	model->on = 'n';
}

int rcl_model_probabilities( RCL_MODEL * model , MNL_DATA * data , double ** P , double * LL )
{
	int k, b, d, m, j, i, base, dbase, bbase;
	
	int I;
	
	double ** Us, ** Es, ** PL , * cc , * bc , * dc ;
	
	if( model == NULL || data  == NULL ) { return 0; }
	
	if( model->K != data->K ) { return -1; }
	if( model->B != data->B ) { return -2; }
	if( model->D != data->D ) { return -3; }
	for( d = 0 ; d < model->D ; d++ ) {
		if( (model->Ld)[d] > (data->Ld)[d] ) {  return -4; }
	}
	if( model->og != data->og ) { return -5; }
	
	// allocate internal space
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	Us = (double **)calloc( data->M , sizeof(double *) );
	Es = (double **)calloc( data->M , sizeof(double *) );
	PL = (double **)calloc( data->M , sizeof(double *) );
	for( m = 0 ; m < data->M ; m++ ) {
		Us[m] = (double *)calloc( I , sizeof(double) );
		Es[m] = (double *)calloc( I , sizeof(double) );
		PL[m] = (double *)calloc( I*((data->Jm)[m]) , sizeof(double) );
	}
	cc = (double *)calloc( I*(data->K) , sizeof(double) ); // continuous coefficients
	bc = (double *)calloc( I*(data->B) , sizeof(double) ); // binary coefficients
	dc = (double *)calloc( I*(data->L) , sizeof(double) ); // dummy coefficients
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	I = 10;
	
	// compute coefficients using sampler function
	(model->smpl)( I , model , cc , bc , dc , NULL );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// compute utilities
	
	base = 0; 
	for( m = 0 ; m < data->M ; m++ ) {
		bbase = 0; // reset binary variable indexer in each market
		for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
			
			// characteristic part
			(PL[m])[I*j+i] = cblas_ddot( data->K , (data->XC)[m] + (data->K)*j , 1 , model->cc , 1 );
			
			// binary variable part
			// number of "on" binary variables for this product: (Bo[m])[j]
			// indices of "on" variables: (XB[m])(bbase,bbase+1,...)
			if( data->B > 0 ) {
				for( b = 0 ; b < ((data->Bo)[m])[j] ; b++ ) {
					(PL[m])[I*j+i] += (model->bc)[ ((data->XB)[m])[ bbase++ ] ];
				}
			}
			
			// dummy variable part
			dbase = 0;
			for( d = 0 ; d < data->D ; d++ ) {
				(PL[m])[I*j+i] += (model->dc)[ ((data->XD)[m])[ (data->D)*j + d] + dbase ];
				dbase += (data->Ld)[d];
			}
			
			(Us[m])[i] = MAX( (PL[m])[I*j+i] , (Us[m])[i] );
			
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
	
	// we can free the coefficients (cc, bc, dc) now, but not sure why
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	if( LL != NULL ) { // compute likelihood, if given
		
		base = 0; LL[0] = 0.0;
		for( m = 0 ; m < data->M ; m++ ) {
			LL[0] += cblas_ddot( (data->Jm)[m] , (data->s)[m] , 1 , U + base , 1 );
			LL[0] -= ( Us[m] + log( Es[m] ) );
			base += (data->Jm)[m];
		}
		
	}
	
	if( PL != NULL ) { // compute mixture probabilities
		
		for( m = 0 ; m < data->M ; m++ ) {
			if( P[m] != NULL ) { 
				if( (data->Jm)[m] == 1 ) { (P[m])[0] = 1.0; }
				else{ 
					for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
						(P[m])[j] = 0.0;
						for( i = 0 ; i < I ; i++ ) { (P[m])[j] += (PL[m])[I*j+i]; }
						(P[m])[j] /= (double)I;
					}
				}
			}
		}
		
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	for( m = 0 ; m < data->M ; m++ ) { free(Us[m]); } free(Us);
	for( m = 0 ; m < data->M ; m++ ) { free(Es[m]); } free(Es);
	for( m = 0 ; m < data->M ; m++ ) { free(PL[m]); } free(PL);
	free(cc);
	free(bc);
	free(dc);
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	return 1;
	
}

void rcl_model_print( MNL_MODEL * model )
{
	int k, b, d, l, dbase;
	
	if( model == NULL ) { return; }
	
	printf("K: %i\n",model->K);
	if( model->K > 0 ) { 
		printf( "  m: %0.16f " , (model->cc_m)[0] );  
		for( k = 1 ; k < model->K ; k++ ) { 
			printf( ", %0.16f " , (model->cc_m)[k] );  
		}
		printf("\n");
		printf( "  v: %0.16f " , (model->cc_v)[0] );  
		for( k = 1 ; k < model->K ; k++ ) { 
			printf( ", %0.16f " , (model->cc_v)[k] );  
		}
		printf("\n");
	}
	
	printf("B: %i\n",model->B);
	if( model->B > 0 ) { 
		printf( "  m: %0.16f " , (model->bc_m)[0] );  
		for( b = 1 ; b < model->B ; b++ ) { 
			printf( ", %0.16f " , (model->bc_m)[b] );  
		}
		printf("\n");
		printf( "  v: %0.16f " , (model->bc_v)[0] );  
		for( b = 1 ; b < model->B ; b++ ) { 
			printf( ", %0.16f " , (model->bc_v)[b] );  
		}
		printf("\n");
	}
	
	printf("D: %i\n",model->D);
	if( model->D > 0 ) { 
		dbase = 0;
		for( d = 0 ; d < model->D ; d++ ) { 
			printf( "  m(%i): %0.16f " , d+1, (model->dc_m)[dbase++] ); 
			for( l = 1 ; l < (model->Ld)[d] ; l++ ) {
				printf( ", %0.16f " , (model->dc_m)[dbase++] );  
			}
			printf("\n");
			printf( "  v(%i): %0.16f " , d+1, (model->dc_m)[dbase++] ); 
			for( l = 1 ; l < (model->Ld)[d] ; l++ ) {
				printf( ", %0.16f " , (model->dc_m)[dbase++] );  
			}
			printf("\n");
		}
	}
	
}

void rcl_model_fprint( MNL_MODEL * model , FILE * fp )
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
		fprintf( fp , ", %0.16f " , (model->cc_m)[0] );  
		for( k = 1 ; k < model->K ; k++ ) { 
			fprintf( fp , ", %0.16f " , (model->cc_m)[k] );  
		}
		fprintf( fp , ", %0.16f " , (model->cc_v)[0] );  
		for( k = 1 ; k < model->K ; k++ ) { 
			fprintf( fp , ", %0.16f " , (model->cc_v)[k] );  
		}
	}
	if( model->B > 0 ) { 
		fprintf( fp , ", %0.16f " , (model->bc_m)[0] );  
		for( b = 1 ; b < model->B ; b++ ) { 
			fprintf( fp , ", %0.16f " , (model->bc_m)[b] );  
		}
		fprintf( fp , ", %0.16f " , (model->bc_v)[0] );  
		for( b = 1 ; b < model->B ; b++ ) { 
			fprintf( fp , ", %0.16f " , (model->bc_v)[b] );  
		}
	}
	if( model->D > 0 ) { 
		fprintf( fp , ", %0.16f " , (model->dc_m)[0] );  
		for( d = 1 ; d < model->L ; d++ ) { 
			fprintf( fp , ", %0.16f " , (model->dc_m)[d] );  
		}
		fprintf( fp , ", %0.16f " , (model->dc_v)[0] );  
		for( d = 1 ; d < model->L ; d++ ) { 
			fprintf( fp , ", %0.16f " , (model->dc_v)[d] );  
		}
	}
}
