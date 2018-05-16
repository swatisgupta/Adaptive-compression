#ifndef  __ACOMPS_UTIL_H__
#define  __ACOMPS_UTIL_H__ 


#include <inttypes.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <assert.h>
#include "plot.h"
#include "acomp_types.h" 

extern char *c_var;
#define IS_REAL 1
typedef enum {
    FALSE = 0,
    TRUE 
} bool_t;	

//#define DEBUG
//#define PLOT    
//#define PRINT_STATS
//#define VERIFY

#ifdef PRINT_STATS
 #define STAT_PRINTF(f_, ...) { printf((f_), __VA_ARGS__); fflush(stdout); } 
#else
 #define STAT_PRINTF(f_, ...) 
#endif
  
#ifdef DEBUG
 #define DB_PRINTF(f_, ...)  { printf((f_), __VA_ARGS__); fflush(stdout); } 
#else
 #define DB_PRINTF(f_, ...) 
#endif

#endif
