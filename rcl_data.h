/*
 *  rcl_data.h
 *  
 *
 *  Created by W. Ross Morrow on 11/3/13.
 *  Copyright 2013 Iowa State University. All rights reserved.
 *
 */

#ifndef RCLDATA_H
#define RCLDATA_H

typedef int (* RCL_SAMPLER)(int I , RCL_MODEL * model , 
							double * cc , double * bc , double * dc , 
							void * params );

int rcl_imn_sampler(int I , RCL_MODEL * model , 
					double * cc , double * bc , double * dc , 
					void * params );
				 
// RCL model data
typedef struct { 
	
	char		on;		// is this model structure "on"? ('y' or 'n')
	char		og;		// "y" or "n" flag for outside good
	
	int			K;		// number of characteristics
	
	int			B;		// number of binary variables
	
	int			D;		// number of dummy variables
	int *		Ld;		// D-vector, number of "levels" per dummy variable
	int			L;		// total number of dummy variables (== sum_d L(d))
	
	// mean coefficients
	double *	cc_m;	// K-vector of continuous coefficients
	double *	bc_m;	// B-vector of binary variable "penalties"
	double *	dc_m;	// L-vector of dummy part-worths
	double		uc_m;	// utility constant, if there is an outside good
	
	// variance coefficients
	double *	cc_v;	// K-vector of continuous coefficients
	double *	bc_v;	// B-vector of binary variable "penalties"
	double *	dc_v;	// L-vector of dummy part-worths
	double		uc_v;	// utility constant, if there is an outside good
	
	// sampler
	RCL_SAMPLER smpl;
	
} RCL_MODEL;

void rcl_model_fprint( RCL_MODEL * model , FILE * fp );

void rcl_model_new( RCL_MODEL * data );
void rcl_model_alloc( RCL_MODEL * model , int K , int B , int D , int * Ld );
void rcl_model_free( RCL_MODEL * data );

// compute probabilities a given model prescribes to a given data set
// model and dataset must be consistent
int rcl_model_probabilities( RCL_MODEL * model , MNL_DATA * data , double ** P , double * loglik );

#endif