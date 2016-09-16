/*
 *  rcl_snopt.h
 *  
 *
 *  Created by W. Ross Morrow on 11/2/13.
 *  Copyright 2013 Iowa State University. All rights reserved.
 *
 */

#ifndef RCLSNOPT_H
#define RCLSNOPT_H

#include "mnl_data.h"

int rcl_mle_snopt(MNL_DATA * estdata , 
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

int rcl_gmm_snopt(MNL_DATA * estdata , 
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