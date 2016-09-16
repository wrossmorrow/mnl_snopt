/*
 *  mnl_snopt.h
 *  
 *
 *  Created by W. Ross Morrow on 10/23/13.
 *  Copyright 2013 Iowa State University. All rights reserved.
 *
 */

#ifndef MNLSNOPT_H
#define MNLSNOPT_H

#include "mnl_data.h"

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
				  char * logfn );

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
				  char * logfn );

#endif