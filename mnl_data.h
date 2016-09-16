/*
 *  mnl_data.h
 *  
 *
 *  Created by W. Ross Morrow on 10/23/13.
 *  Copyright 2013 Iowa State University. All rights reserved.
 *
 */

#ifndef MNLDATA_H
#define MNLDATA_H

typedef struct { 
	
	int			on;		// "on" (or "off") flag. denotes whether structure has been allocated
	
	int			M;		// number of markets
	int *		ml;		// market "labels"
	int *		Jm;		// M-vector, number of products in each market
	int			J;		// number of product-observations (== sum_m J(m))
	
	int			K;		// number of (continous) characteristics
	
	int			B;		// number of binary variables
	int **		Bo;		// M-vector of J(m) vectors, numbers of binary variables that are "on" (= 1) per product-observation
	int	*		BoN;	// M-vector, total number of "on" binary variables for each market
	int			BN;		// total number of "on" binary variables for each product in each market
	
	int			D;		// number of dummy variables
	int	*		Ld;		// D-vector, number of "levels" per dummy variable
	int			L;		// total number of dummy variables (== sum_d L(d))
	
	double **	XC;		// M-vector of K x J(m) matrices (column-major) of continuous characteristics
	int **		XB;		// M-vector of BoN(m)-vectors of indices (from {1,...,B}) of "on" binary variables
	int	**		XD;		// M-vector of D x J(m) matrices (column-major) of dummy level values for each product-observation
	
	double **	s;		// M-vector of J(m)-vectors of product market shares. Sum to one within markets if there is no outside good, and less than one otherwise. 
	
	char		og;		// "y" or "n" flag for outside good
	double *	so;		// M-vector, outside good shares for each market. Shares within markets plus these terms sum to one. 
	
} MNL_DATA;

void mnl_data_new( MNL_DATA * data );
void mnl_data_alloc( MNL_DATA * data );
void mnl_data_free( MNL_DATA * data );
void mnl_data_print( MNL_DATA * data );
int  mnl_data_read(FILE * dfp,			// file to read from
				int headrows ,		// header rows to skip
				char coderow ,		// "y" or "n" flag denoting whether id codes are in the first row read
				int colnum ,		// number of columns
				char * coltype ,	// colnum-vector of column types: 
				MNL_DATA * data);	// data file to read into


void printprodvec( int K , double * xc , int B , int Bo , int * xb , int D , int * xd );

typedef int (* MNL_DATA_MAP)( int K_s , double * xc_s , int B_s , int   Bo_s , int * xb_s , int D_s , int * xd_s , 
							  int K_d , double * xc_d , int B_d , int * Bo_d , int * xb_d , int D_d , int * xd_d , 
							  void * params );

int mnl_data_map( MNL_DATA * source , MNL_DATA * dest , MNL_DATA_MAP map , void * params );

int mnl_data_check( MNL_DATA * data , int * rank , double * norm );
				 
// MNL model data
typedef struct { 
	
	char		on;		// is this model structure "on"? 
	
	int			K;		// number of characteristics
	
	int			B;		// number of binary variables
	
	int			D;		// number of dummy variables
	int *		Ld;		// D-vector, number of "levels" per dummy variable
	int			L;		// tota number of dummy variables (== sum_d L(d))
	
	double *	cc;		// K-vector of continuous coefficients
	double *	bc;		// B-vector of binary variable "penalties"
	double *	dc;		// L-vector of dummy part-worths
	
	char		og;		// "y" or "n" flag for outside good
	double		uc;		// utility constant, if there is an outside good
	
} MNL_MODEL;

void mnl_model_fprint( MNL_MODEL * model , FILE * fp );

void mnl_model_new( MNL_MODEL * data );
void mnl_model_alloc( MNL_MODEL * model , int K , int B , int D , int * Ld );
void mnl_model_free( MNL_MODEL * data );

// compute probabilities a given model prescribes to a given data set
// model and dataset must be consistent
int mnl_model_probabilities( MNL_MODEL * model , MNL_DATA * data , double ** PL , double * loglik );
int mnl_model_probabilities( MNL_MODEL * model , MNL_DATA * data , double ** PL , double * loglik );

#endif