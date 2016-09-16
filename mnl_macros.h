/*
 *  mnl_utility.h
 *  
 *
 *  Created by W. Ross Morrow on 10/23/13.
 *  Copyright 2013 Iowa State University. All rights reserved.
 *
 */

#ifndef XPS_UTILITY_H
#define XPS_UTILITY_H

#include <math.h>

// these macros are copies from Todd Munson's PATH software Macros.h header

#if !defined(MAX)
#  define MAX(a,b)  ((a) > (b) ? (a) : (b))
#endif

#if !defined(MIN)
#  define MIN(a,b)  ((a) < (b) ? (a) : (b))
#endif

#if !defined(MID)
#  define MID(a,b,c)  (((a) < (b)) ?                   \
(((b) < (c)) ? (b) : (c)) :     \
(((a) < (c)) ? (a) : (c)))
#endif

#if !defined(ABS)
#  define ABS(a)    (fabs(a))
#endif

#if !defined(SGN)
#  define SGN(a)    ((a) > 0 ? (1) : (-1))
#endif

#endif