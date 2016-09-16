/*
 *  mnl_snopt.c
 *  
 *
 *  Created by W. Ross Morrow on 10/23/13.
 *  Copyright 2013 Iowa State University. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Use these includes on wrmorrow-01,wrmorrow-02
//
// #include <vecLib/cblas.h>
// #include <vecLib/clapack.h>

// Use this include on wrmorrow-MBP
//
#include <Accelerate/Accelerate.h>

#include <unuran.h>
#include <unuran_urng_rngstreams.h>

#include "f2c.h"
#include "snfilewrapper.h"
#include "snopt.h"
#include "snoextras.h"

#include "mnl_macros.h"
#include "mnl_data.h"
#include "mnl_snopt.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static MNL_DATA * data;
static MNL_MODEL * model;

static int b_v_start;
static int p_v_start;
static int g_v_start;
static int u_v_start;
static int x_v_start;

static int ob_r_start;
static int dv_r_start;
static int ue_r_start;
static int re_r_start;

static double * XCL;	// K-vector: lower bounds on characteristics
static double * XCU;	// K-vector: upper bounds on characteristics
static double * XCM;	// K-vector: midpoints of characteristic ranges

static double * Us;	// M-vector, maximum in-market utilities
static double * E;	// J-vector, exponentiated utilities
static double * Es;	// M-vector, exponentiated utility sums

static char initcond;
static int scalechar;

static UNUR_URNG *  unrng;
static UNUR_GEN *   urngen;

static clock_t tmp_ticks;
static clock_t func_eval_ticks;

void printl( char * filename , char * message )
{
	FILE * fp  = fopen( filename , "a" );
	if( fp != NULL ) { fprintf( fp , message ); fclose( fp ); }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * MULTINOMIAL LOGIT MAXIMUM LIKELIHOOD ESTIMATION * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This function attempts to fit a Multinomial Logit model to share data from
 * multiple markets using maximum likelihood. 
 * 
 * Each attempt solves the following problem: 
 * 
 *	max		s' u - sum_{m=1}^M log( sum_{j=1}^{J(m)} exp( u(m,j) ) )
 *	
 *	wrt		b(1),...,b(K)										K characteristic coefficients
 *			p(1),...,p(B)										B binary variable coefficients
 *			g(1,1),...,g(1,L(1)),...,g(D,1),...,g(D,L(D))		L = sum_d L(d) dummy variable coefficients
 *			u(1,1),...,u(1,J(1)),...,u(M,1),...,u(M,J(M))		J = sum_m J(m) variables
 * 
 *	sto		g(d,1) + ... + g(d,L(d)) = 0	for all d			D linear constraints
 *			XC' b + p(XB) + g(XD) - u = 0						J linear constraints
 * 
 * "p(XB)" and "g(XD)" above denote linear "selection" operators wherein the right  
 * binary "penalties" (p) and dummy "part-worths" (g) are chosen based on the relevant
 * index data for the products. 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * ARGUMENTS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * RETURNS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * NOTES * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void mnl_mle_size( integer * N , integer * M , integer * Annz , integer * Gnnz )
{
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// number of variables (see formulation above for details) and C-style indexers
	
	N[0]  = 0;
	
	N[0] += data->K;	b_v_start = 0;						// beta (characteristic) coefficients
	N[0] += data->B;	p_v_start = b_v_start + data->K;	// rho (binary "penalty") coefficients
	N[0] += data->L;	g_v_start = p_v_start + data->B;	// gamma (dummy variable) coefficients
	N[0] += data->J;	u_v_start = g_v_start + data->L;	// utilities
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// number of constraints (see formulation above for details) 
	
	M[0]  = 0;			ob_r_start = 0;						// (always have an objective)
	M[0] += data->D;	dv_r_start = 1;						// D dummy variable constraints
	M[0] += data->J;	ue_r_start = dv_r_start + data->D;	// J utility equation constraints
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// number of nonzeros in linear part:
	
	Annz[0]  = 0;					// objective is nonlinear
	Annz[0] += data->L;				// total of L terms in D dummy variable constraints
	Annz[0] += (data->K)*(data->J); // K characteristic terms in J utility equations
	Annz[0] += data->BN;			// BN binary variable terms (one for each "on" variable)
	Annz[0] += (data->D)*(data->J); // D dummy terms in J utility equations (one nonzero per dummy, per equation)
	Annz[0] += data->J;				// J utility variables in J utility equations
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// number of objective gradient and constraint Jacobian non-zeros
	
	Gnnz[0]  = data->J;	// objective (utility derivatives)
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
}

void mnl_mle_bounds( integer N , integer M , double * xLoBnds , double * xUpBnds , double * FLoBnds , double * FUpBnds )
{
	int n;
	
	// no bounds
	
	for( n = 0 ; n < N ; n++ ) { xLoBnds[n] = -1.0e20; xUpBnds[n] = 1.0e20; }
	
	// objective has no bounds (but will be negative)
	
	FLoBnds[ob_r_start] = -1.0e20;	FUpBnds[ob_r_start] = 1.0e20;
	
	// other constraints are all negative-null form
	
	for( n = 1 ; n <= M ; n++ ) { FLoBnds[n] = 0.0;	FUpBnds[n] = 0.0; }
	
}

void mnl_mle_initcond( integer N , double * x0 ) 
{
	int n, m, j, k, b, d, base, bbase, dbase;
	double * u0;
	
	u0 = x0 + u_v_start;
	
	switch( initcond ) {
			
		case 'm': 
			
			// take coefficients from (given) model structure; consistency check done elsewhere
			for( n = 0 ; n < model->K ; n++ ) { x0[b_v_start+n] = (model->cc)[n]; }
			for( n = 0 ; n < model->B ; n++ ) { x0[p_v_start+n] = (model->bc)[n]; }
			for( n = 0 ; n < model->L ; n++ ) { x0[g_v_start+n] = (model->dc)[n]; }
			
			// define utilities consistently
			base = 0; 
			for( m = 0 ; m < data->M ; m++ ) {
				bbase = 0; // reset binary variable indexer in each market
				for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
					
					// characteristic part
					u0[base] = cblas_ddot( data->K , (data->XC)[m] + (data->K)*j , 1 , model->cc , 1 );
					
					// binary variable part
					if( data->B > 0 ) {
						for( b = 0 ; b < ((data->Bo)[m])[j] ; b++ ) {
							u0[base] += (model->bc)[ ((data->XB)[m])[ bbase++ ] ];
						}
					}
					
					// dummy variable part
					dbase = 0;
					for( d = 0 ; d < data->D ; d++ ) {
						u0[base] += (model->dc)[ ((data->XD)[m])[ (data->D)*j + d] + dbase ];
						dbase += (data->Ld)[d];
					}
					
					// increment base
					base++;
				}
			}
			break;
			
		default: 
			for( n = 0 ; n < N ; n++ ) { x0[n] = 20.0 * unur_sample_cont(urngen) - 10.0; }
			break;
			
	}
	
}

void mnl_mle_linmap( integer N , integer M , integer Annz , integer * Arows , integer * Acols , double * Adata )
{
	int d, l, k, b, m, j, Abase, jbase, bbase, ibase;
	
	// linear part nonzero pattern and elements
	
	Abase = 0;
	
	// dummy variable constraints (sum-to-one within variables)
	
	ibase = 0;
	for( d = 0 ; d < data->D ; d++ ) {
		for( l = 0 ; l < (data->Ld)[d] ; l++ ) {
			Arows[Abase] = dv_r_start + d;
			Acols[Abase] = g_v_start + ibase + l;
			Adata[Abase] = 1.0;
			Abase++;
		}
		ibase += (data->Ld)[d];
	}
	
	// utility equation constraints: characteristic terms
	
	if( data->K > 0 ) { // don't do any of this if K == 0
		
		if( scalechar ) {
			
			// find characteristics range and midpoint to modify characteristics in problem
			for( k = 0 ; k < data->K ; k++ ) {
				XCL[k] = 1.0e20; XCU[k] = -1.0e20;
				for( m = 0 ; m < data->M ; m++ ) {
					for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
						XCL[k] = MIN( ((data->XC)[m])[ (data->K)*j + k ] , XCL[k] );
						XCU[k] = MAX( ((data->XC)[m])[ (data->K)*j + k ] , XCU[k] );
					}
				}
				XCM[k] = 0.5 * ( XCU[k] - XCL[k] );
			}
			
			jbase = 0;
			for( m = 0 ; m < data->M ; m++ ) {
				
				// create linear characteristics part of utility equations with transformed coordinates
				for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
					for( k = 0 ; k < data->K ; k++ ) {
						Arows[Abase] = ue_r_start + jbase;
						Acols[Abase] = b_v_start  + k;
						Adata[Abase] = ( ((data->XC)[m])[ (data->K)*j + k ] - XCL[k] ) / XCM[k] - 1.0;
						Abase++;
					}
					jbase++;
				}
			}
			
		} else {
			
			jbase = 0;
			for( m = 0 ; m < data->M ; m++ ) {
				
				// create linear characteristics part of utility equations with transformed coordinates
				for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
					for( k = 0 ; k < data->K ; k++ ) {
						Arows[Abase] = ue_r_start + jbase;
						Acols[Abase] = b_v_start  + k;
						Adata[Abase] = ((data->XC)[m])[ (data->K)*j + k ]; 
						Abase++;
					}
					jbase++;
				}
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
					Arows[Abase] = ue_r_start + jbase;
					Acols[Abase] = p_v_start  + ((data->XB)[m])[ bbase ];
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
				for( d = 0 ; d < data->D ; d++ ) {
					Arows[Abase] = ue_r_start + jbase;
					Acols[Abase] = g_v_start  + ibase + ((data->XD)[m])[ (data->D)*j + d ];
					Adata[Abase] = 1.0;
					Abase++;
					ibase += (data->Ld)[d];
				}
				jbase++;
			}
		}
	}
	
	// utility equation constraints: utility variable terms
	
	for( j = 0 ; j < data->J ; j++ ) {
		Arows[Abase] = ue_r_start + j;
		Acols[Abase] = u_v_start + j;
		Adata[Abase] = - 1.0;
		Abase++;
	}
	
	// enforce FORTRAN-style indexing
	for( d = 0 ; d < Annz ; d++ ) { Arows[d] += 1; Acols[d] += 1; }
}

void mnl_mle_structure( integer N , integer M , integer Gnnz , integer * Grows , integer * Gcols )
{
	int m, j, base;
	
	// nonlinear part Jacobian sparsity pattern
	
	base = 0;
	
	for( m = 0 ; m < data->M ; m++ ) {
		for( j = 0 ; j < data->Jm[m] ; j++ ) {
			Grows[base] = ob_r_start;
			Gcols[base] = u_v_start + base;
			base++;
		}
	}
	
	// enforce FORTRAN-style indexing
	for( m = 0 ; m < Gnnz ; m++ ) { Grows[m] += 1; Gcols[m] += 1; }
	
}

int mnl_mle_callback_(integer		*Status,	// SNOPT status code
				  integer		*N,			// number of variables
				  doublereal	*x,			// current variable values
				  integer		*needF,		// 0 if f(x) is not needed, > 0 if it is
				  integer		*FN,		// length of the vector of objective and constraint values
				  doublereal	*Fvals,		// values (to be calculated) for objective and constraint values
				  integer		*needG,     // 0 if G(x) not needed, > 0 if it is
				  integer		*Gnnz,		// length of arrays iGvar and jGfun
				  doublereal	*Gvals,		// derivative values (MMF format)
				  char			*cu,		// character workspace
				  integer		*lencu,		// length of character workspace
				  integer		*iu,		// integer workspace
				  integer		*leniu,		// length of integer workspace
				  doublereal	*ru,		// double workspace
				  integer		*lenru )	// length of double workspace
{
	int m, j, base;
	
	double * u = x + u_v_start;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	base = 0;
	for( m = 0 ; m < data->M ; m++ ) {
		Us[m] = 0.0;
		for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
			Us[m] = MAX( u[base] , Us[m] ); 
			base++; // can't use ++ in the macro it seems
		}
	}
	
	base = 0;
	for( m = 0 ; m < data->M ; m++ ) {
		Es[m] = 0.0;
		if( Us[m] > 0.0 ) {
			for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
				E[base] = exp( u[base] - Us[m] ); Es[m] += E[base]; base++;
			}
		} else {
			for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
				E[base] = exp( u[base] ); Es[m] += E[base]; base++;
			}
		}
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	if( needF[0] > 0 ) {
		
		Fvals[0] = 0.0; base = 0;
		for( m = 0 ; m < data->M ; m++ ) {
			Fvals[0] += cblas_ddot( (data->Jm)[m] , (data->s)[m] , 1 , u + base , 1 );
			Fvals[0] -= ( Us[m] + log( Es[m] ) ); // Us[m] + log( sum_j exp( u[j] - Us[m] ) ) = log( sum_j exp( u[j] ) )
			base += (data->Jm)[m];
		}
		
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	if( needG[0] > 0 ) {
		
		base = 0;
		for( m = 0 ; m < data->M ; m++ ) {
			cblas_dcopy( (data->Jm)[m] , (data->s)[m] , 1 , Gvals + base , 1 );
			base += (data->Jm)[m];
		}
		
		base = 0;
		for( m = 0 ; m < data->M ; m++ ) {
			for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
				Gvals[base] -= E[base] / Es[m]; 
				base++;
			}
		}
		
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
    return 0;
	
}

// set-up, solve, and post-process SNOPT formulation
int mnl_mle_snopt(MNL_DATA * estdata , 
				  MNL_MODEL * estmodel , 
				  int trials , 
				  double * loglik , 
				  double * rates , 
				  int scale , 
				  int checkders , 
				  char ic , 
				  double * OptTol , 
				  double * FeasTol , 
				  char * logfn )
{
	integer			Start = 0; // Cold = 0, Basis = 1, Warm = 2;
	
	integer			l;
	
    integer          N, M, FN;
	
	integer			 Annz;
	double			*Adata;
	integer			*Arows, *Acols;
	
	integer			 Gnnz;
    double			*Gdata;
	integer			*Grows, *Gcols;
	
	integer			*xState;
	double			*x, *xLoBnds, *xUpBnds, *lambdaB;
	
	integer			*FState;
	double			*Fvals, *FLoBnds, *FUpBnds, *lambdaC;
	
	double			ObjAdd; // constant to add to objective
	integer			ObjRow; // row of objective in combined function
	integer			INFO;	// 
	
	integer			minrw, miniw, mincw;
	
	// USER workspace
	
	// real (double) workspace
	integer			lenru = 500;
	double			ru[8*500];
	
	// integer workspace
	integer			leniu = 500;
	integer			iu[8*500]; 
	
	// char workspace
	integer			lencu = 500;
	char			cu[8*500];
	
	// SNOPT workspace
	// we initialize to 500, and re-allocate below
	
	// real (double) workspace
	integer			lenrw = 500;
	double			*rw;
	
	// integer workspace
	integer			leniw = 500;
	integer			*iw;
	
	// char workspace
	integer			lencw = 500;
	char			*cw;
	
	integer			nxName = 1; // do not use variable names
	char			xNames[1*8];
	
	integer			nFName = 1; // do not use constraint names
	char			FNames[1*8];
	
	integer			iPrint = 9; // "unit number for the print file"
	integer			iSumm  = 6; // "unit number for the Summary file"
	
	integer			prnt_len; // 
	integer			iSpecs = 4,  spec_len;
	
	char			Prob[200];
	char			printname[200];
	char			specname[200];
	
	integer			nS, nInf, npname = 1;
	double			sInf;
	
	integer			iSum, iPrt, strOpt_len;
	char			strOpt[200];
	
	int besti;
	double bestll;
	double * bestx;
	
	clock_t ticks;
	
	int n, k, b, d, t, tr, base;
	
	FILE * resfp = NULL;
	FILE * logfp = NULL;
	
	UNUR_DISTR * distr;
	UNUR_PAR *   param;
	
	// The RNGSTREAMS library sets a package seed. 
	unsigned long seed[] = {111u, 222u, 333u, 444u, 555u, 666u};
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	if( estdata == NULL || estmodel == NULL ) { return -1; }
	
	data = estdata; model = estmodel;
	
	scalechar = scale; initcond = ic;
	
	if( ic == 'm' ) { 
		
		// verify model-data consistency to use model initial conditions
		if( model->K != data->K ) { return -1; }
		if( model->B != data->B ) { return -2; }
		if( model->D != data->D ) { return -3; }
		for( d = 0 ; d < model->D ; d++ ) {
			if( (model->Ld)[d] > (data->Ld)[d] ) {  return -4; }
		}
		if( model->og != data->og ) { return -5; }
		
	}
	
    Us = ( double * ) calloc ( data->M , sizeof( double ) );
    E  = ( double * ) calloc ( data->J , sizeof( double ) );
    Es = ( double * ) calloc ( data->M , sizeof( double ) );
	
	if( data->K > 0 && scalechar ) {
		XCL = ( double * ) calloc ( data->K , sizeof( double ) );
		XCU = ( double * ) calloc ( data->K , sizeof( double ) );
		XCM = ( double * ) calloc ( data->K , sizeof( double ) );
	} else { XCL = NULL; XCU = NULL; XCM = NULL; }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * CREATE RANDOM NUMBER GENERATORS * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	RngStream_SetPackageSeed(seed);
	
	// Make an object for uniform random number generator.
	unrng = unur_urng_rngstream_new("unrng");
	if( unrng == NULL ) { 
		printf("ERROR - Uniform generator could not be constructed.\n");
		exit(EXIT_FAILURE); 
	}
	
	// use predefined distribution - uniform
	distr = unur_distr_uniform(NULL, 0);
	
	// use "auto" method (why, not sure)
	param = unur_auto_new(distr);
	
	// Set uniform generator in parameter object
	unur_set_urng(param, unrng);
	
	// Create the uniform random number generator object.
	urngen = unur_init(param); // param is "destroyed" here
	if( urngen == NULL ) { 
		printf("ERROR - Uniform generator could not be constructed.\n");
	}
	
	// free distribution object
	unur_distr_free(distr);
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	mnl_mle_size( &N , &M , &Annz , &Gnnz ); FN = M + 1;
	
	logfp = fopen( logfn , "a" );
	if( logfp != NULL ) {
		fprintf(logfp,"Problem Sizes:\n");
		fprintf(logfp,"  %i variables\n",(int)N);
		fprintf(logfp,"  %i constraints\n",(int)M);
		fprintf(logfp,"  %i linear part nonzeros (%0.1f%% dense)\n",(int)Annz,100.0*(double)Annz/((double)(N*M)));
		fprintf(logfp,"  %i Jacobian nonzeros (%0.1f%% dense)\n",(int)Gnnz,100.0*(double)Gnnz/((double)(N*M)));
		fclose( logfp );
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// problem variables
    x			 = (double *) calloc (N, sizeof(double));
    bestx    	 = (double *) calloc (N, sizeof(double));
    xState		 = (integer *)calloc (N, sizeof(integer)); // initialized to zero for no information
    xLoBnds      = (double *) calloc (N, sizeof(double));
    xUpBnds      = (double *) calloc (N, sizeof(double));
    lambdaB      = (double *) calloc (N, sizeof(double)); // bounds multipliers
	
	// combined function (objective and constraints)
    Fvals		 = (double *) calloc (FN, sizeof(double)); // initializes to zero
    FState       = (integer *)calloc (FN, sizeof(integer));
    FLoBnds      = (double *) calloc (FN, sizeof(double)); // initializes to zero
    FUpBnds      = (double *) calloc (FN, sizeof(double)); // initializes to zero
    lambdaC      = (double *) calloc (FN, sizeof(double)); // constraint multipliers
	
	// linear part of the objective and constraints
    Adata		 = (double *) calloc (Annz, sizeof(double));
    Arows		 = (integer *)calloc (Annz, sizeof(integer));
    Acols		 = (integer *)calloc (Annz, sizeof(integer));
	
	// Jacobian of the nonlinear part of objective and constraints
    Gdata		 = (double *) calloc (Gnnz, sizeof(double));
    Grows		 = (integer *)calloc (Gnnz, sizeof(integer));
    Gcols		 = (integer *)calloc (Gnnz, sizeof(integer));
	
	// initial SNOPT workspace; resized below
	cw			 = (char*)   calloc(8*lencw,sizeof(char   ));
	iw			 = (integer*)calloc(  leniw,sizeof(integer));
	rw			 = (double*) calloc(  lenrw,sizeof(double ));
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * INITIALIZE FORTRAN-STYLE FILE REFERENCES  * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// label spec (options) file, using SNOPT's FORTRAN utilities
	sprintf(specname ,   "%s", "mle_snopt.spc"); spec_len = strlen(specname);
	
	// open Print file, using SNOPT's FORTRAN utilities
	sprintf(printname,   "%s", "mle_snopt.out"); prnt_len = strlen(printname);
	snopenappend_( &iPrint, printname,   &INFO, prnt_len );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * INITIALIZE SNOPT  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	printl( logfn , "Starting SNOPT...\n" );
	
	sninit_( &iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500 );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * SNOPT MEMORY ALLOCATION * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	snmema_(&INFO,
			&FN,
			&N,
			&nxName, &nFName,
			&Annz, &Gnnz, 
			&mincw, &miniw, &minrw, 
			cw, &lencw, 
			iw, &leniw, 
			rw, &lenrw, 
			8*500);
	
	// if memory was NOT sized successfully, 
	if( INFO != 104 ) {
		
		printl( logfn , "WARNING:: SNOPT could not estimate memory requirements correctly.\n" );
		
	} else {
		
		printl( logfn , "SNOPT estimated memory requirements.\n" );
		// fprintf(logfp,"SNOPT estimated memory requirements: %i, %i, %i.\n",(int)mincw,(int)miniw,(int)minrw);
		
		// re-initializing SNOPT workspace, if needed
		
		if( lencw < mincw ) { 
			lencw = mincw; 
			cw = (char*)realloc(cw, 8*lencw*sizeof(char));
		}
		
		if( leniw < miniw ) {
			leniw = miniw;
			iw = (integer*)realloc(iw, leniw*sizeof(integer));
		}
		
		if( lenrw < minrw ) {
			lenrw = minrw;
			rw = (double*) realloc(rw, lenrw*sizeof(double));
		}
		
		printl( logfn , "Re-initializating SNOPT.\n" );
		// fprintf(logfp,"Re-initializating SNOPT with sizes (%li,%li,%li)\n",lencw,leniw,lenrw);
			
		sninit_( &iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500 );
		
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * INITIALIZE SNOPT OPTIONS  * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// options
	snfilewrapper_(specname, 
				   &iSpecs, 
				   &INFO, 
				   cw, &lencw,
				   iw, &leniw, 
				   rw, &lenrw, 
				   spec_len,
				   8*lencw);
	
	if( INFO != 101 ) {
		printl( logfn , "WARNING: trouble reading specs file. Using default options.\n" );
    }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * INITIALIZE SNOPT PROBLEM DATA * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	Start = 0;													// cold start
	
	strcpy(Prob,"MNLMLE");										// Problem name
	
	ObjRow = 1; ObjAdd = 0.0;									// objective information
	
	mnl_mle_bounds( N , M , xLoBnds , xUpBnds , FLoBnds , FUpBnds ); // define bounds
	printl( logfn , "defined bounds...\n" );
	
	mnl_mle_linmap( N , M , Annz , Arows , Acols , Adata );		// define linear part
	printl( logfn , "defined linear map...\n" );
	
	mnl_mle_structure( N , M , Gnnz , Grows , Gcols );			// define nonlinear Jacobian structure
	printl( logfn , "defined Jacobian structure...\n" );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	
	printf("Variables (%li):\n",N);
	for( n = 0 ; n < N ; n++ ) {
		printf("  %0.4f <= x_%i <= %0.4f\n",xLoBnds[n],n+1,xUpBnds[n]);
	}
	
	printf("Combined Map: \n");
	for( n = 0 ; n < FN ; n++ ) {
		printf("  %0.8f <= F(%i) <= %0.8f\n",FLoBnds[n],n+1,FUpBnds[n]);
	}
	
	printf("Linear Part: \n");
	for( n = 0 ; n < Annz ; n++ ) {
		printf("  A(%i,%i) = %0.8f\n",(int)Arows[n],(int)Acols[n],Adata[n]);
	}
	
	printf("Nonlinear Part Jacobian Structure: \n");
	for( n = 0 ; n < Gnnz ; n++ ) {
		printf("  Jdata(%i) ~ J(%i,%i)\n",n+1,(int)Grows[n],(int)Gcols[n]);
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	printl( logfn , "\n\nMNL:: Starting...\n" );
	
	if( rates != NULL ) { 
		rates[0] = 0.0; rates[1] = 0.0; rates[2] = 0.0; rates[3] = 0.0; 
	}
	
	besti = 0; bestll = - 1.0e20; // guaranteed to replace bestll if successful
	
	for( t = 0 ; t < trials ; t++ ) {
		
		// get initial condition
		mnl_mle_initcond( N , x );
		
		// optional derivative check
		if( checkders ) {
			snopta_eval_G_check(N, FN,
								mnl_mle_callback_,
								Gnnz, Grows, Gcols, 
								x, xLoBnds, xUpBnds,
								cu, &lencu, 
								iu, &leniu,
								ru, &lenru);
			exit(1);
		}
		
		// (timed) solve
		ticks = clock();
		
		snopta_(&Start,								// 0: Cold, 1: Basis, 2: Warm
				&FN,								// FN = M + 1
				&N,									// number of variables
				&nxName,							// 1 if no names used, otherwise == N
				&nFName,							// 1 if no names used, otherwise == FN
				&ObjAdd,							// scalar to add to objective (for printing)
				&ObjRow,							// row of the objective in F (FORTRAN-Style)
				Prob,								// problem name
				mnl_mle_callback_,						// combined function F(x) needed by snopta_()
				Arows, Acols, &Annz, &Annz, Adata,	// sparse "A" matrix such that F(x) = userfun_(x) + Ax
				Grows, Gcols, &Gnnz, &Gnnz,			// Jacobian structure for G(x) = DF(x)
				xLoBnds, xUpBnds, xNames,			// lower bounds, upper bounds, and names for x variables
				FLoBnds, FUpBnds, FNames,			// lower bounds, upper bounds, and names for F values
				x, xState, lambdaB,					// x values, "states" (see pg. 18), and associated dual variables (multipliers)
				Fvals, FState, lambdaC,				// F values, "states" (see pg. 18?), and associated dual variables (multipliers)
				&INFO,								// result of call to snopta_(). See docs, pg. 19 for details
				&mincw, &miniw, &minrw,				// minimum values of SNOPT workspace sizes for snopta_()
				&nS, &nInf, &sInf,					// see docs, pg. 18 & 20
				cu, &lencu,							// character user workspace
				iu, &leniu,							// integer user workspace
				ru, &lenru,							// real user workspace
				cw, &lencw,							// character SNOPT workspace (at leat 500 + N + NF if names used)
				iw, &leniw,							// integer SNOPT workspace; minimum values given by snmema_()
				rw, &lenrw,							// real SNOPT workspace; minimum values given by snmema_()
				npname, 8*nxName, 8*nFName,
				8*500, 
				8*500);
		
		ticks = clock() - ticks;
		
		// print out rates
		switch( INFO ) { 
				
			case 1:	 
				if( rates != NULL ) { rates[0]++; } 
				
				besti = 1;
				if( Fvals[0] > bestll ) {
					bestll = Fvals[0]; 
					cblas_dcopy( N , x , 1 , bestx , 1 );
				}
				
				break;
				
			case 3:  if( rates != NULL ) { rates[1]++; } break;
			case 41: if( rates != NULL ) { rates[2]++; } break;
			default: if( rates != NULL ) { rates[3]++; } break;
		}
		
		// write out solve status
		logfp = fopen( logfn , "a" );
		if( logfp != NULL ) {
			if( INFO != 1 ) {
				fprintf(logfp,"MNL:: SNOPT failed to solve the problem, final status = %d\n", (int)INFO);
			} else {
				
				fprintf(logfp,"MNL:: SNOPT successful!\n");
				fprintf(logfp,"      Log likelihood: %0.6f\n", Fvals[0]);
				fprintf(logfp,"      Took %0.6f seconds\n", ticks/(double)CLOCKS_PER_SEC);
				
				// print out solution to log file
				if( data->K > 0 ) {
					fprintf(logfp,"      beta coefficients: %0.6e ",x[b_v_start]);
					for( k = 1 ; k < data->K ; k++ ) {
						fprintf(logfp,", %0.6e ",x[b_v_start+k]);
					}
					fprintf(logfp,"\n");
				} else { fprintf(logfp,"      no beta coefficients\n"); }
				
				if( data->B > 0 ) {
					fprintf(logfp,"      rho coefficients: %0.6e ",x[p_v_start]);
					for( b = 1 ; b < data->B ; b++ ) {
						fprintf(logfp,", %0.6e ",x[p_v_start+b]);
					}
					fprintf(logfp,"\n");
				} else { fprintf(logfp,"      no rho coefficients\n"); }
				
				if( data->D > 0 ) {
					fprintf(logfp,"      gamma coefficients: ");
					b = 0;
					for( d = 0 ; d < data->D ; d++ ) {
						fprintf(logfp," ( %0.6e ",x[g_v_start+b++]);
						for( l = 1 ; l < (data->Ld)[d] ; l++ ) {
							fprintf(logfp,", %0.6e ",x[g_v_start+b++]);
						}
						fprintf(logfp,") ");
					}
					fprintf(logfp,"\n");
				} else { fprintf(logfp,"      no gamma coefficients\n"); }
				
			}
			
			fprintf(logfp,"\n");
			
			fclose( logfp );
		}
		
	}
	
	if( rates != NULL ) {
		rates[0] /= 1.0*trials; 
		rates[1] /= 1.0*trials; 
		rates[2] /= 1.0*trials; 
		rates[3] /= 1.0*trials;
	}
	
	// parse multistart results
	
	switch( besti ) {
			
		case 1:
			
			// refine best solution? Only if we are passed both tolerance variables
			if( OptTol != NULL && FeasTol != NULL ) {
				
				strcpy( strOpt , "Major optimality tolerance" ); 
				strOpt_len = strlen( strOpt );
				sngetr_(strOpt, 
						OptTol,
						&INFO,
						cw, &lencw,
						iw, &leniw, 
						rw, &lenrw,
						strOpt_len, 
						8*strOpt_len);
				printf("maj opt tol: %2e\n",OptTol[0]);
				
				strcpy( strOpt , "Major feasibility tolerance" ); 
				strOpt_len = strlen( strOpt );
				sngetr_(strOpt, 
						FeasTol,
						&INFO,
						cw, &lencw,
						iw, &leniw, 
						rw, &lenrw,
						strOpt_len, 
						8*strOpt_len);
				printf("maj fea tol: %2e\n",FeasTol[0]);
				
				cblas_dcopy( N , bestx , 1 , x , 1 );
				INFO = 1;
				while( INFO == 1 ) {
					
					OptTol[0] *= 0.1; FeasTol[0] *= 0.1;
					
					strcpy( strOpt , "Major optimality tolerance" ); 
					strOpt_len = strlen( strOpt );
					snsetr_(strOpt, 
							OptTol,
							&iPrint, &iSumm, &INFO , 
							cw, &lencw,
							iw, &leniw, 
							rw, &lenrw,
							strOpt_len, 
							8*strOpt_len);
					
					strcpy( strOpt , "Major feasibility tolerance" ); 
					strOpt_len = strlen( strOpt );
					snsetr_(strOpt, 
							FeasTol,
							&iPrint, &iSumm, &INFO , 
							cw, &lencw,
							iw, &leniw, 
							rw, &lenrw,
							strOpt_len, 
							8*strOpt_len);
					
					// attempt a solve
					snopta_(&Start,								// 0: Cold, 1: Basis, 2: Warm
							&FN,								// FN = M + 1
							&N,									// number of variables
							&nxName,							// 1 if no names used, otherwise == N
							&nFName,							// 1 if no names used, otherwise == FN
							&ObjAdd,							// scalar to add to objective (for printing)
							&ObjRow,							// row of the objective in F (FORTRAN-Style)
							Prob,								// problem name
							mnl_mle_callback_,						// combined function F(x) needed by snopta_()
							Arows, Acols, &Annz, &Annz, Adata,	// sparse "A" matrix such that F(x) = userfun_(x) + Ax
							Grows, Gcols, &Gnnz, &Gnnz,			// Jacobian structure for G(x) = DF(x)
							xLoBnds, xUpBnds, xNames,			// lower bounds, upper bounds, and names for x variables
							FLoBnds, FUpBnds, FNames,			// lower bounds, upper bounds, and names for F values
							x, xState, lambdaB,					// x values, "states" (see pg. 18), and associated dual variables (multipliers)
							Fvals, FState, lambdaC,				// F values, "states" (see pg. 18?), and associated dual variables (multipliers)
							&INFO,								// result of call to snopta_(). See docs, pg. 19 for details
							&mincw, &miniw, &minrw,				// minimum values of SNOPT workspace sizes for snopta_()
							&nS, &nInf, &sInf,					// see docs, pg. 18 & 20
							cu, &lencu,							// character user workspace
							iu, &leniu,							// integer user workspace
							ru, &lenru,							// real user workspace
							cw, &lencw,							// character SNOPT workspace (at leat 500 + N + NF if names used)
							iw, &leniw,							// integer SNOPT workspace; minimum values given by snmema_()
							rw, &lenrw,							// real SNOPT workspace; minimum values given by snmema_()
							npname, 8*nxName, 8*nFName,
							8*500, 
							8*500);
					
					if( INFO == 1 ) { 
						cblas_dcopy( N , x , 1 , bestx , 1 );
						bestll = Fvals[0];
					} else {
						OptTol[0] *= 10.0; FeasTol[0] *= 10.0;
					}
					
				}
				
			}
			
			// success: store solve data in model structure. We presume the 
			// model structure is formatted correctly, relative to the given
			// data. 
			
			// return likelihood value
			if( loglik != NULL ) { loglik[0] = bestll; }
			
			// free data in model structure, if it is "on" (might have been used for initial conditions)
			if( model->on == 'y' ) { mnl_model_free( model ); }
			mnl_model_alloc( model , data->K , data->B , data->D , data->Ld );
			
			// store characteristic coefficients
			if( data->K > 0 ) {
				if( scalechar ) { // transform coefficients out of unit ([-1,1]) space
					for( k = 0 ; k < data->K ; k++ ) { bestx[b_v_start+k] /= XCM[k]; }
				}
				cblas_dcopy( model->K , bestx + b_v_start , 1 , model->cc , 1 );
			}
			
			// store binary variable "penalties"
			if( data->B > 0 ) {
				cblas_dcopy( model->B , bestx + p_v_start , 1 , model->bc , 1 );
			}
			
			// store dummy variable coefficients
			if( data->D > 0 ) {
				cblas_dcopy( model->L , bestx + g_v_start , 1 , model->dc , 1 );
			}
			
			break;
			
		default: break;
			
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	free(x);
    free(xState);
    free(xLoBnds);
    free(xUpBnds);
	free(lambdaB);
	
	free(Fvals);
    free(FState);
    free(FLoBnds);
    free(FUpBnds);
	free(lambdaC);
	
	free(Adata);
    free(Arows);
    free(Acols);
	
	free(Gdata);
    free(Grows);
    free(Gcols);
	
	// SNOPT workspace
	free(cw);
	free(iw);
	free(rw);
	
	// close print and specs files
	snclose_( &iPrint );
	snclose_( &iSpecs );
	
	free(bestx);
	
	if( scalechar ) {
		if( data->K > 0 ) { free( XCL ); free( XCU ); free( XCM ); }
	}
	
	free(Us);
	free(E);
	free(Es);
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	return besti;
	
}




/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * MULTINOMIAL LOGIT "GMM" ESTIMATION (LOGISTIC REGRESSION)  * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This function attempts to fit a Multinomial Logit model to share data from
 * multiple markets using maximum likelihood. 
 * 
 * Each attempt solves the following problem: 
 * 
 *	min		r' r / 2
 *	
 *	wrt		b(1),...,b(K)										K characteristic coefficients
 *			p(1),...,p(B)										B binary variable coefficients
 *			g(1,1),...,g(1,L(1)),...,g(D,1),...,g(D,L(D))		L = sum_d L(d) dummy variable coefficients
 *			u(1,1),...,u(1,J(1)),...,u(M,1),...,u(M,J(M))		J = sum_m J(m) variables
 *			r(1,1),...,r(1,J(1)),...,r(M,1),...,r(M,J(M))		J = sum_m J(m) variables
 * 
 *	sto		 g(d,1) + ... + g(d,L(d)) = 0	for all d			D linear constraints
 *			XC' b + p(XB) + g(XD) - u = 0						J linear constraints
 *			                    u + r = log( s )				J linear constraints
 * 
 * "p(XB)" and "g(XD)" above denote linear "selection" operators wherein the right  
 * binary "penalties" (p) and dummy "part-worths" (g) are chosen based on the relevant
 * index data for the products. 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * ARGUMENTS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * RETURNS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * NOTES * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void mnl_gmm_size( integer * N , integer * M , integer * Annz , integer * Gnnz )
{
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// number of variables (see formulation above for details) and C-style indexers
	
	N[0]  = 0;
	
	N[0] += data->K;	b_v_start = 0;						// beta (characteristic) coefficients
	N[0] += data->B;	p_v_start = b_v_start + data->K;	// rho (binary "penalty") coefficients
	N[0] += data->L;	g_v_start = p_v_start + data->B;	// gamma (dummy variable) coefficients
	N[0] += data->J;	u_v_start = g_v_start + data->L;	// utilities
	N[0] += data->J;	x_v_start = u_v_start + data->J;	// utility residuals
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// number of constraints (see formulation above for details) 
	
	M[0]  = 0;			ob_r_start = 0;						// (always have an objective)
	M[0] += data->D;	dv_r_start = 1;						// D dummy variable constraints
	M[0] += data->J;	ue_r_start = dv_r_start + data->D;	// J utility equation constraints
	M[0] += data->J;	re_r_start = ue_r_start + data->J;	// J utility residual equations
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// number of nonzeros in linear part:
	
	Annz[0]  = 0;					// objective is nonlinear
	Annz[0] += data->L;				// total of L terms in D dummy variable constraints
	Annz[0] += (data->K)*(data->J); // K characteristic terms in J utility equations
	Annz[0] += data->BN;			// BN binary variable terms (one for each "on" variable)
	Annz[0] += (data->D)*(data->J); // D dummy terms in J utility equations (one nonzero per dummy, per equation)
	Annz[0] += data->J;				// J utility variables in J utility equations
	Annz[0] += data->J;				// J utility variables in J utility residual equations
	Annz[0] += data->J;				// J utility residuals in J utility residual equations
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// number of objective gradient and constraint Jacobian non-zeros
	
	Gnnz[0]  = data->J;	// objective (utility residual derivatives)
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
}

void mnl_gmm_bounds( integer N , integer M , double * xLoBnds , double * xUpBnds , double * FLoBnds , double * FUpBnds )
{
	int n, m, j, d, base;
	double tmp;
	
	// no bounds
	
	for( n = 0 ; n < N ; n++ ) { xLoBnds[n] = -1.0e20; xUpBnds[n] = 1.0e20; }
	
	// objective has no bounds (but will be negative)
	
	FLoBnds[ob_r_start] = -1.0e20;	FUpBnds[ob_r_start] = 1.0e20;
	
	// constraints have a different form
	
	// sum-to-zero constraints have null form
	for( d = 0 ; d < data->D ; d++ ) { 
		FLoBnds[dv_r_start+d] = 0.0; FUpBnds[dv_r_start+d] = 0.0; 
	}
	
	// utility equations are null form, residual equations have log-share RHSs
	base = 0;
	for( m = 0 ; m < data->M ; m++ ) { 
		for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
			
			FLoBnds[ue_r_start+base] = 0.0; 
			FUpBnds[ue_r_start+base] = 0.0;
			
			tmp = log( ((data->s)[m])[j] );
			FLoBnds[re_r_start+base] = tmp; 
			FUpBnds[re_r_start+base] = tmp; 
			
			base++;
			
		}
	}
	
}

void mnl_gmm_initcond( integer N , double * x0 ) 
{
	int n, m, j, k, b, d, base, bbase, dbase;
	double * u0;
	
	u0 = x0 + u_v_start;
	
	switch( initcond ) {
			
		case 'm': 
			
			// take coefficients from (given) model structure; consistency check done elsewhere
			for( n = 0 ; n < model->K ; n++ ) { x0[b_v_start+n] = (model->cc)[n]; }
			for( n = 0 ; n < model->B ; n++ ) { x0[p_v_start+n] = (model->bc)[n]; }
			for( n = 0 ; n < model->L ; n++ ) { x0[g_v_start+n] = (model->dc)[n]; }
			
			// define utilities and utility residuals consistently
			base = 0; 
			for( m = 0 ; m < data->M ; m++ ) {
				bbase = 0; // reset binary variable indexer in each market
				for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
					
					// characteristic part
					u0[base] = cblas_ddot( data->K , (data->XC)[m] + (data->K)*j , 1 , model->cc , 1 );
					
					// binary variable part
					if( data->B > 0 ) {
						for( b = 0 ; b < ((data->Bo)[m])[j] ; b++ ) {
							u0[base] += (model->bc)[ ((data->XB)[m])[ bbase++ ] ];
						}
					}
					
					// dummy variable part
					dbase = 0;
					for( d = 0 ; d < data->D ; d++ ) {
						u0[base] += (model->dc)[ ((data->XD)[m])[ (data->D)*j + d] + dbase ];
						dbase += (data->Ld)[d];
					}
					
					x0[x_v_start+base] = log( ((data->s)[m])[j] ) - u0[base];
					
					// increment base
					base++;
				}
			}
			
			
			break;
			
		default: 
			for( n = 0 ; n < N ; n++ ) { x0[n] = 20.0 * unur_sample_cont(urngen) - 10.0; }
			break;
			
	}
	
}

void mnl_gmm_linmap( integer N , integer M , integer Annz , integer * Arows , integer * Acols , double * Adata )
{
	int d, l, k, b, m, j, Abase, jbase, bbase, ibase;
	
	// linear part nonzero pattern and elements
	
	Abase = 0;
	
	// dummy variable constraints (sum-to-one within variables)
	
	ibase = 0;
	for( d = 0 ; d < data->D ; d++ ) {
		for( l = 0 ; l < (data->Ld)[d] ; l++ ) {
			Arows[Abase] = dv_r_start + d;
			Acols[Abase] = g_v_start + ibase + l;
			Adata[Abase] = 1.0;
			Abase++;
		}
		ibase += (data->Ld)[d];
	}
	
	// utility equation constraints: characteristic terms
	
	if( data->K > 0 ) { // don't do any of this if K == 0
		
		if( scalechar ) {
			
			// find characteristics range and midpoint to modify characteristics in problem
			for( k = 0 ; k < data->K ; k++ ) {
				XCL[k] = 1.0e20; XCU[k] = -1.0e20;
				for( m = 0 ; m < data->M ; m++ ) {
					for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
						XCL[k] = MIN( ((data->XC)[m])[ (data->K)*j + k ] , XCL[k] );
						XCU[k] = MAX( ((data->XC)[m])[ (data->K)*j + k ] , XCU[k] );
					}
				}
				XCM[k] = 0.5 * ( XCU[k] - XCL[k] );
			}
			
			jbase = 0;
			for( m = 0 ; m < data->M ; m++ ) {
				
				// create linear characteristics part of utility equations with transformed coordinates
				for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
					for( k = 0 ; k < data->K ; k++ ) {
						Arows[Abase] = ue_r_start + jbase;
						Acols[Abase] = b_v_start  + k;
						Adata[Abase] = ( ((data->XC)[m])[ (data->K)*j + k ] - XCL[k] ) / XCM[k] - 1.0;
						Abase++;
					}
					jbase++;
				}
			}
			
		} else {
			
			jbase = 0;
			for( m = 0 ; m < data->M ; m++ ) {
				
				// create linear characteristics part of utility equations with transformed coordinates
				for( j = 0 ; j < (data->Jm)[m] ; j++ ) {
					for( k = 0 ; k < data->K ; k++ ) {
						Arows[Abase] = ue_r_start + jbase;
						Acols[Abase] = b_v_start  + k;
						Adata[Abase] = ((data->XC)[m])[ (data->K)*j + k ]; 
						Abase++;
					}
					jbase++;
				}
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
					printf(" ");
					Arows[Abase] = ue_r_start + jbase;
					Acols[Abase] = p_v_start  + ((data->XB)[m])[ bbase ];
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
				for( d = 0 ; d < data->D ; d++ ) {
					Arows[Abase] = ue_r_start + jbase;
					Acols[Abase] = g_v_start  + ibase + ((data->XD)[m])[ (data->D)*j + d ];
					Adata[Abase] = 1.0;
					Abase++;
					ibase += (data->Ld)[d];
				}
				jbase++;
			}
		}
	}
	
	// utility equation constraints: utility variable terms
	
	for( j = 0 ; j < data->J ; j++ ) {
		Arows[Abase] = ue_r_start + j;
		Acols[Abase] = u_v_start  + j;
		Adata[Abase] = - 1.0;
		Abase++;
	}
	
	// utility residual constraints: utility variable terms
	
	for( j = 0 ; j < data->J ; j++ ) {
		Arows[Abase] = re_r_start + j;
		Acols[Abase] = u_v_start  + j;
		Adata[Abase] = 1.0;
		Abase++;
	}
	
	// utility residual constraints: utility residual terms
	
	for( j = 0 ; j < data->J ; j++ ) {
		Arows[Abase] = re_r_start + j;
		Acols[Abase] = x_v_start  + j;
		Adata[Abase] = 1.0;
		Abase++;
	}
	
	// enforce FORTRAN-style indexing
	for( d = 0 ; d < Annz ; d++ ) { Arows[d] += 1; Acols[d] += 1; }
}

void mnl_gmm_structure( integer N , integer M , integer Gnnz , integer * Grows , integer * Gcols )
{
	int j;
	
	// nonlinear part Jacobian sparsity pattern
	for( j = 0 ; j < data->J ; j++ ) {
		Grows[j] = ob_r_start;
		Gcols[j] = x_v_start + j;
	}
	
	// enforce FORTRAN-style indexing
	for( j = 0 ; j < Gnnz ; j++ ) { Grows[j] += 1; Gcols[j] += 1; }
	
}

int mnl_gmm_callback_(integer		*Status,	// SNOPT status code
					  integer		*N,			// number of variables
					  doublereal	*x,			// current variable values
					  integer		*needF,		// 0 if f(x) is not needed, > 0 if it is
					  integer		*FN,		// length of the vector of objective and constraint values
					  doublereal	*Fvals,		// values (to be calculated) for objective and constraint values
					  integer		*needG,     // 0 if G(x) not needed, > 0 if it is
					  integer		*Gnnz,		// length of arrays iGvar and jGfun
					  doublereal	*Gvals,		// derivative values (MMF format)
					  char			*cu,		// character workspace
					  integer		*lencu,		// length of character workspace
					  integer		*iu,		// integer workspace
					  integer		*leniu,		// length of integer workspace
					  doublereal	*ru,		// double workspace
					  integer		*lenru )	// length of double workspace
{
	if( needF[0] > 0 ) {
		Fvals[0] = 0.5 * cblas_ddot( data->J , x + x_v_start , 1 , x + x_v_start , 1 );
	}
	
	if( needG[0] > 0 ) { 
		cblas_dcopy( data->J , x + x_v_start , 1 , Gvals , 1 ); 
	}
	
    return 0;
}

// set-up, solve, and post-process SNOPT formulation
int mnl_gmm_snopt(MNL_DATA * estdata , 
				  MNL_MODEL * estmodel , 
				  int trials , 
				  double * gmmres , 
				  double * rates , 
				  int scale , 
				  int checkders , 
				  char ic , 
				  double * OptTol , 
				  double * FeasTol , 
				  char * logfn )
{
	integer			Start = 0; // Cold = 0, Basis = 1, Warm = 2;
	
	integer			l;
	
    integer          N, M, FN;
	
	integer			 Annz;
	double			*Adata;
	integer			*Arows, *Acols;
	
	integer			 Gnnz;
    double			*Gdata;
	integer			*Grows, *Gcols;
	
	integer			*xState;
	double			*x, *xLoBnds, *xUpBnds, *lambdaB;
	
	integer			*FState;
	double			*Fvals, *FLoBnds, *FUpBnds, *lambdaC;
	
	double			ObjAdd; // constant to add to objective
	integer			ObjRow; // row of objective in combined function
	integer			INFO;	// 
	
	integer			minrw, miniw, mincw;
	
	// USER workspace
	
	// real (double) workspace
	integer			lenru = 500;
	double			ru[8*500];
	
	// integer workspace
	integer			leniu = 500;
	integer			iu[8*500]; 
	
	// char workspace
	integer			lencu = 500;
	char			cu[8*500];
	
	// SNOPT workspace
	// we initialize to 500, and re-allocate below
	
	// real (double) workspace
	integer			lenrw = 500;
	double			*rw;
	
	// integer workspace
	integer			leniw = 500;
	integer			*iw;
	
	// char workspace
	integer			lencw = 500;
	char			*cw;
	
	integer			nxName = 1; // do not use variable names
	char			xNames[1*8];
	
	integer			nFName = 1; // do not use constraint names
	char			FNames[1*8];
	
	integer			iPrint = 9; // "unit number for the print file"
	integer			iSumm  = 6; // "unit number for the Summary file"
	
	integer			prnt_len; // 
	integer			iSpecs = 4,  spec_len;
	
	char			Prob[200];
	char			printname[200];
	char			specname[200];
	
	integer			nS, nInf, npname = 1;
	double			sInf;
	
	integer			iSum, iPrt, strOpt_len;
	char			strOpt[200];
	
	int besti;
	double bestll;
	double * bestx;
	
	clock_t ticks;
	
	int n, k, b, d, t, tr, base;
	
	FILE * resfp = NULL;
	FILE * logfp = NULL;
	
	UNUR_DISTR * distr;
	UNUR_PAR *   param;
	
	// The RNGSTREAMS library sets a package seed. 
	unsigned long seed[] = {111u, 222u, 333u, 444u, 555u, 666u};
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	if( estdata == NULL || estmodel == NULL ) { return -1; }
	
	data = estdata; model = estmodel;
	
	scalechar = scale; initcond = ic;
	
	if( ic == 'm' ) { 
		
		// verify model-data consistency to use model initial conditions
		if( model->K != data->K ) { return -1; }
		if( model->B != data->B ) { return -2; }
		if( model->D != data->D ) { return -3; }
		for( d = 0 ; d < model->D ; d++ ) {
			if( (model->Ld)[d] > (data->Ld)[d] ) {  return -4; }
		}
		if( model->og != data->og ) { return -5; }
		
	}
	
    Us = ( double * ) calloc ( data->M , sizeof( double ) );
    E  = ( double * ) calloc ( data->J , sizeof( double ) );
    Es = ( double * ) calloc ( data->M , sizeof( double ) );
	
	if( data->K > 0 && scalechar ) {
		XCL = ( double * ) calloc ( data->K , sizeof( double ) );
		XCU = ( double * ) calloc ( data->K , sizeof( double ) );
		XCM = ( double * ) calloc ( data->K , sizeof( double ) );
	} else { XCL = NULL; XCU = NULL; XCM = NULL; }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * CREATE RANDOM NUMBER GENERATORS * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	RngStream_SetPackageSeed(seed);
	
	// Make an object for uniform random number generator.
	unrng = unur_urng_rngstream_new("unrng");
	if( unrng == NULL ) { 
		printf("ERROR - Uniform generator could not be constructed.\n");
		exit(EXIT_FAILURE); 
	}
	
	// use predefined distribution - uniform
	distr = unur_distr_uniform(NULL, 0);
	
	// use "auto" method (why, not sure)
	param = unur_auto_new(distr);
	
	// Set uniform generator in parameter object
	unur_set_urng(param, unrng);
	
	// Create the uniform random number generator object.
	urngen = unur_init(param); // param is "destroyed" here
	if( urngen == NULL ) { 
		printf("ERROR - Uniform generator could not be constructed.\n");
	}
	
	// free distribution object
	unur_distr_free(distr);
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	mnl_gmm_size( &N , &M , &Annz , &Gnnz ); FN = M + 1;
	
	logfp = fopen( logfn , "a" );
	if( logfp != NULL ) {
		fprintf(logfp,"Problem Sizes:\n");
		fprintf(logfp,"  %i variables\n",(int)N);
		fprintf(logfp,"  %i constraints\n",(int)M);
		fprintf(logfp,"  %i linear part nonzeros (%0.1f%% dense)\n",(int)Annz,100.0*(double)Annz/((double)(N*M)));
		fprintf(logfp,"  %i Jacobian nonzeros (%0.1f%% dense)\n",(int)Gnnz,100.0*(double)Gnnz/((double)(N*M)));
		fclose( logfp );
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// problem variables
    x			 = (double *) calloc (N, sizeof(double));
    bestx    	 = (double *) calloc (N, sizeof(double));
    xState		 = (integer *)calloc (N, sizeof(integer)); // initialized to zero for no information
    xLoBnds      = (double *) calloc (N, sizeof(double));
    xUpBnds      = (double *) calloc (N, sizeof(double));
    lambdaB      = (double *) calloc (N, sizeof(double)); // bounds multipliers
	
	// combined function (objective and constraints)
    Fvals		 = (double *) calloc (FN, sizeof(double)); // initializes to zero
    FState       = (integer *)calloc (FN, sizeof(integer));
    FLoBnds      = (double *) calloc (FN, sizeof(double)); // initializes to zero
    FUpBnds      = (double *) calloc (FN, sizeof(double)); // initializes to zero
    lambdaC      = (double *) calloc (FN, sizeof(double)); // constraint multipliers
	
	// linear part of the objective and constraints
    Adata		 = (double *) calloc (Annz, sizeof(double));
    Arows		 = (integer *)calloc (Annz, sizeof(integer));
    Acols		 = (integer *)calloc (Annz, sizeof(integer));
	
	// Jacobian of the nonlinear part of objective and constraints
    Gdata		 = (double *) calloc (Gnnz, sizeof(double));
    Grows		 = (integer *)calloc (Gnnz, sizeof(integer));
    Gcols		 = (integer *)calloc (Gnnz, sizeof(integer));
	
	// initial SNOPT workspace; resized below
	cw			 = (char*)   calloc(8*lencw,sizeof(char   ));
	iw			 = (integer*)calloc(  leniw,sizeof(integer));
	rw			 = (double*) calloc(  lenrw,sizeof(double ));
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * INITIALIZE FORTRAN-STYLE FILE REFERENCES  * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// label spec (options) file, using SNOPT's FORTRAN utilities
	sprintf(specname ,   "%s", "gmm_snopt.spc"); spec_len = strlen(specname);
	
	// open Print file, using SNOPT's FORTRAN utilities
	sprintf(printname,   "%s", "gmm_snopt.out"); prnt_len = strlen(printname);
	snopenappend_( &iPrint, printname,   &INFO, prnt_len );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * INITIALIZE SNOPT  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	printl( logfn , "Starting SNOPT...\n" );
	
	sninit_( &iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500 );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * SNOPT MEMORY ALLOCATION * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	snmema_(&INFO,
			&FN,
			&N,
			&nxName, &nFName,
			&Annz, &Gnnz, 
			&mincw, &miniw, &minrw, 
			cw, &lencw, 
			iw, &leniw, 
			rw, &lenrw, 
			8*500);
	
	// if memory was NOT sized successfully, 
	if( INFO != 104 ) {
		
		printl( logfn , "WARNING:: SNOPT could not estimate memory requirements correctly.\n" );
		
	} else {
		
		printl( logfn , "SNOPT estimated memory requirements.\n" );
		// fprintf(logfp,"SNOPT estimated memory requirements: %i, %i, %i.\n",(int)mincw,(int)miniw,(int)minrw);
		
		// re-initializing SNOPT workspace, if needed
		
		if( lencw < mincw ) { 
			lencw = mincw; 
			cw = (char*)realloc(cw, 8*lencw*sizeof(char));
		}
		
		if( leniw < miniw ) {
			leniw = miniw;
			iw = (integer*)realloc(iw, leniw*sizeof(integer));
		}
		
		if( lenrw < minrw ) {
			lenrw = minrw;
			rw = (double*) realloc(rw, lenrw*sizeof(double));
		}
		
		printl( logfn , "Re-initializating SNOPT.\n" );
		// fprintf(logfp,"Re-initializating SNOPT with sizes (%li,%li,%li)\n",lencw,leniw,lenrw);
		
		sninit_( &iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500 );
		
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * INITIALIZE SNOPT OPTIONS  * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// options
	snfilewrapper_(specname, 
				   &iSpecs, 
				   &INFO, 
				   cw, &lencw,
				   iw, &leniw, 
				   rw, &lenrw, 
				   spec_len,
				   8*lencw);
	
	if( INFO != 101 ) {
		printl( logfn , "WARNING: trouble reading specs file. Using default options.\n" );
    }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * INITIALIZE SNOPT PROBLEM DATA * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	Start = 0;													// cold start
	
	strcpy(Prob,"MNLGMM");										// Problem name
	
	ObjRow = 1; ObjAdd = 0.0;									// objective information
	
	mnl_gmm_bounds( N , M , xLoBnds , xUpBnds , FLoBnds , FUpBnds ); // define bounds
	printl( logfn , "defined bounds...\n" );
	
	mnl_gmm_linmap( N , M , Annz , Arows , Acols , Adata );		// define linear part
	printl( logfn , "defined linear map...\n" );
	
	mnl_gmm_structure( N , M , Gnnz , Grows , Gcols );			// define nonlinear Jacobian structure
	printl( logfn , "defined Jacobian structure...\n" );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 
	 printf("Variables (%li):\n",N);
	 for( n = 0 ; n < N ; n++ ) {
	 printf("  %0.4f <= x_%i <= %0.4f\n",xLoBnds[n],n+1,xUpBnds[n]);
	 }
	 
	 printf("Combined Map: \n");
	 for( n = 0 ; n < FN ; n++ ) {
	 printf("  %0.8f <= F(%i) <= %0.8f\n",FLoBnds[n],n+1,FUpBnds[n]);
	 }
	 
	 printf("Linear Part: \n");
	 for( n = 0 ; n < Annz ; n++ ) {
	 printf("  A(%i,%i) = %0.8f\n",(int)Arows[n],(int)Acols[n],Adata[n]);
	 }
	 
	 printf("Nonlinear Part Jacobian Structure: \n");
	 for( n = 0 ; n < Gnnz ; n++ ) {
	 printf("  Jdata(%i) ~ J(%i,%i)\n",n+1,(int)Grows[n],(int)Gcols[n]);
	 }
	
	// exit(1);
	 
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	printl( logfn , "\n\nMNL:: Starting...\n" );
	
	if( rates != NULL ) { 
		rates[0] = 0.0; rates[1] = 0.0; rates[2] = 0.0; rates[3] = 0.0; 
	}
	
	besti = 0; bestll = 1.0e20; // guaranteed to replace bestll if successful
	
	for( t = 0 ; t < trials ; t++ ) {
		
		// get initial condition
		mnl_gmm_initcond( N , x );
		
		// optional derivative check
		if( checkders ) {
			snopta_eval_G_check(N, FN,
								mnl_gmm_callback_,
								Gnnz, Grows, Gcols, 
								x, xLoBnds, xUpBnds,
								cu, &lencu, 
								iu, &leniu,
								ru, &lenru);
			exit(1);
		}
		
		// (timed) solve
		ticks = clock();
		
		snopta_(&Start,								// 0: Cold, 1: Basis, 2: Warm
				&FN,								// FN = M + 1
				&N,									// number of variables
				&nxName,							// 1 if no names used, otherwise == N
				&nFName,							// 1 if no names used, otherwise == FN
				&ObjAdd,							// scalar to add to objective (for printing)
				&ObjRow,							// row of the objective in F (FORTRAN-Style)
				Prob,								// problem name
				mnl_gmm_callback_,					// combined function F(x) needed by snopta_()
				Arows, Acols, &Annz, &Annz, Adata,	// sparse "A" matrix such that F(x) = userfun_(x) + Ax
				Grows, Gcols, &Gnnz, &Gnnz,			// Jacobian structure for G(x) = DF(x)
				xLoBnds, xUpBnds, xNames,			// lower bounds, upper bounds, and names for x variables
				FLoBnds, FUpBnds, FNames,			// lower bounds, upper bounds, and names for F values
				x, xState, lambdaB,					// x values, "states" (see pg. 18), and associated dual variables (multipliers)
				Fvals, FState, lambdaC,				// F values, "states" (see pg. 18?), and associated dual variables (multipliers)
				&INFO,								// result of call to snopta_(). See docs, pg. 19 for details
				&mincw, &miniw, &minrw,				// minimum values of SNOPT workspace sizes for snopta_()
				&nS, &nInf, &sInf,					// see docs, pg. 18 & 20
				cu, &lencu,							// character user workspace
				iu, &leniu,							// integer user workspace
				ru, &lenru,							// real user workspace
				cw, &lencw,							// character SNOPT workspace (at leat 500 + N + NF if names used)
				iw, &leniw,							// integer SNOPT workspace; minimum values given by snmema_()
				rw, &lenrw,							// real SNOPT workspace; minimum values given by snmema_()
				npname, 8*nxName, 8*nFName,
				8*500, 
				8*500);
		
		ticks = clock() - ticks;
		
		// print out rates
		switch( INFO ) { 
				
			case 1:	 
				if( rates != NULL ) { rates[0]++; } 
				
				besti = 1;
				if( Fvals[0] < bestll ) {
					bestll = Fvals[0]; 
					cblas_dcopy( N , x , 1 , bestx , 1 );
				}
				
				break;
				
			case 3:  if( rates != NULL ) { rates[1]++; } break;
			case 41: if( rates != NULL ) { rates[2]++; } break;
			default: if( rates != NULL ) { rates[3]++; } break;
		}
		
		// write out solve status
		logfp = fopen( logfn , "a" );
		if( logfp != NULL ) {
			if( INFO != 1 ) {
				fprintf(logfp,"MNL:: SNOPT failed to solve the problem, final status = %d\n", (int)INFO);
			} else {
				
				fprintf(logfp,"MNL:: SNOPT successful!\n");
				fprintf(logfp,"      GMM residual: %0.6f\n", Fvals[0]);
				fprintf(logfp,"      Took %0.6f seconds\n", ticks/(double)CLOCKS_PER_SEC);
				
				// print out solution to log file
				if( data->K > 0 ) {
					fprintf(logfp,"      beta coefficients: %0.6e ",x[b_v_start]);
					for( k = 1 ; k < data->K ; k++ ) {
						fprintf(logfp,", %0.6e ",x[b_v_start+k]);
					}
					fprintf(logfp,"\n");
				} else { fprintf(logfp,"      no beta coefficients\n"); }
				
				if( data->B > 0 ) {
					fprintf(logfp,"      rho coefficients: %0.6e ",x[p_v_start]);
					for( b = 1 ; b < data->B ; b++ ) {
						fprintf(logfp,", %0.6e ",x[p_v_start+b]);
					}
					fprintf(logfp,"\n");
				} else { fprintf(logfp,"      no rho coefficients\n"); }
				
				if( data->D > 0 ) {
					fprintf(logfp,"      gamma coefficients: ");
					b = 0;
					for( d = 0 ; d < data->D ; d++ ) {
						fprintf(logfp," ( %0.6e ",x[g_v_start+b++]);
						for( l = 1 ; l < (data->Ld)[d] ; l++ ) {
							fprintf(logfp,", %0.6e ",x[g_v_start+b++]);
						}
						fprintf(logfp,") ");
					}
					fprintf(logfp,"\n");
				} else { fprintf(logfp,"      no gamma coefficients\n"); }
				
			}
			
			fprintf(logfp,"\n");
			
			fclose( logfp );
		}
		
	}
	
	if( rates != NULL ) {
		rates[0] /= 1.0*trials; 
		rates[1] /= 1.0*trials; 
		rates[2] /= 1.0*trials; 
		rates[3] /= 1.0*trials;
	}
	
	// parse multistart results
	
	switch( besti ) {
			
		case 1:
			
			// refine best solution? Only if we are passed both tolerance variables
			if( OptTol != NULL && FeasTol != NULL ) {
				
				strcpy( strOpt , "Major optimality tolerance" ); 
				strOpt_len = strlen( strOpt );
				sngetr_(strOpt, 
						OptTol,
						&INFO,
						cw, &lencw,
						iw, &leniw, 
						rw, &lenrw,
						strOpt_len, 
						8*strOpt_len);
				printf("maj opt tol: %2e\n",OptTol[0]);
				
				strcpy( strOpt , "Major feasibility tolerance" ); 
				strOpt_len = strlen( strOpt );
				sngetr_(strOpt, 
						FeasTol,
						&INFO,
						cw, &lencw,
						iw, &leniw, 
						rw, &lenrw,
						strOpt_len, 
						8*strOpt_len);
				printf("maj fea tol: %2e\n",FeasTol[0]);
				
				cblas_dcopy( N , bestx , 1 , x , 1 );
				INFO = 1;
				while( INFO == 1 ) {
					
					OptTol[0] *= 0.1; FeasTol[0] *= 0.1;
					
					strcpy( strOpt , "Major optimality tolerance" ); 
					strOpt_len = strlen( strOpt );
					snsetr_(strOpt, 
							OptTol,
							&iPrint, &iSumm, &INFO , 
							cw, &lencw,
							iw, &leniw, 
							rw, &lenrw,
							strOpt_len, 
							8*strOpt_len);
					
					strcpy( strOpt , "Major feasibility tolerance" ); 
					strOpt_len = strlen( strOpt );
					snsetr_(strOpt, 
							FeasTol,
							&iPrint, &iSumm, &INFO , 
							cw, &lencw,
							iw, &leniw, 
							rw, &lenrw,
							strOpt_len, 
							8*strOpt_len);
					
					// attempt a solve
					snopta_(&Start,								// 0: Cold, 1: Basis, 2: Warm
							&FN,								// FN = M + 1
							&N,									// number of variables
							&nxName,							// 1 if no names used, otherwise == N
							&nFName,							// 1 if no names used, otherwise == FN
							&ObjAdd,							// scalar to add to objective (for printing)
							&ObjRow,							// row of the objective in F (FORTRAN-Style)
							Prob,								// problem name
							mnl_gmm_callback_,					// combined function F(x) needed by snopta_()
							Arows, Acols, &Annz, &Annz, Adata,	// sparse "A" matrix such that F(x) = userfun_(x) + Ax
							Grows, Gcols, &Gnnz, &Gnnz,			// Jacobian structure for G(x) = DF(x)
							xLoBnds, xUpBnds, xNames,			// lower bounds, upper bounds, and names for x variables
							FLoBnds, FUpBnds, FNames,			// lower bounds, upper bounds, and names for F values
							x, xState, lambdaB,					// x values, "states" (see pg. 18), and associated dual variables (multipliers)
							Fvals, FState, lambdaC,				// F values, "states" (see pg. 18?), and associated dual variables (multipliers)
							&INFO,								// result of call to snopta_(). See docs, pg. 19 for details
							&mincw, &miniw, &minrw,				// minimum values of SNOPT workspace sizes for snopta_()
							&nS, &nInf, &sInf,					// see docs, pg. 18 & 20
							cu, &lencu,							// character user workspace
							iu, &leniu,							// integer user workspace
							ru, &lenru,							// real user workspace
							cw, &lencw,							// character SNOPT workspace (at leat 500 + N + NF if names used)
							iw, &leniw,							// integer SNOPT workspace; minimum values given by snmema_()
							rw, &lenrw,							// real SNOPT workspace; minimum values given by snmema_()
							npname, 8*nxName, 8*nFName,
							8*500, 
							8*500);
					
					if( INFO == 1 ) { 
						cblas_dcopy( N , x , 1 , bestx , 1 );
						bestll = Fvals[0];
					} else {
						OptTol[0] *= 10.0; FeasTol[0] *= 10.0;
					}
					
				}
				
			}
			
			// success: store solve data in model structure. We presume the 
			// model structure is formatted correctly, relative to the given
			// data. 
			
			// return likelihood value
			if( gmmres != NULL ) { gmmres[0] = bestll; }
			
			// free data in model structure, if it is "on" (might have been used for initial conditions)
			if( model->on == 'y' ) { mnl_model_free( model ); }
			mnl_model_alloc( model , data->K , data->B , data->D , data->Ld );
			
			// store characteristic coefficients
			if( data->K > 0 ) {
				if( scalechar ) { // transform coefficients out of unit ([-1,1]) space
					for( k = 0 ; k < data->K ; k++ ) { bestx[b_v_start+k] /= XCM[k]; }
				}
				cblas_dcopy( model->K , bestx + b_v_start , 1 , model->cc , 1 );
			}
			
			// store binary variable "penalties"
			if( data->B > 0 ) {
				cblas_dcopy( model->B , bestx + p_v_start , 1 , model->bc , 1 );
			}
			
			// store dummy variable coefficients
			if( data->D > 0 ) {
				cblas_dcopy( model->L , bestx + g_v_start , 1 , model->dc , 1 );
			}
			
			break;
			
		default: break;
			
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	free(x);
    free(xState);
    free(xLoBnds);
    free(xUpBnds);
	free(lambdaB);
	
	free(Fvals);
    free(FState);
    free(FLoBnds);
    free(FUpBnds);
	free(lambdaC);
	
	free(Adata);
    free(Arows);
    free(Acols);
	
	free(Gdata);
    free(Grows);
    free(Gcols);
	
	// SNOPT workspace
	free(cw);
	free(iw);
	free(rw);
	
	// close print and specs files
	snclose_( &iPrint );
	snclose_( &iSpecs );
	
	free(bestx);
	
	if( scalechar ) {
		if( data->K > 0 ) { free( XCL ); free( XCU ); free( XCM ); }
	}
	
	free(Us);
	free(E);
	free(Es);
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	return besti;
	
}