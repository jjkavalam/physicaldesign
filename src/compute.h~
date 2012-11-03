#ifndef _PLACER_H_
#define _PLACER_H_

#include <gsl/gsl_multimin.h>
#include <stdio.h>
#include <stdlib.h>
#include "parser/bookshelf_IO.h"
#include "blas_sparse.h"
#include "spblasi_matrix_double.h"
#include "spblasi_table.h"
#include <string.h>


// for x: Indices start at 1 
extern double *D, *B, *d_x, *u, *x, *c;
extern blas_sparse_matrix Z, C;
extern int nRegions;

extern int *Partition;

typedef struct {
// Represent the floorplan
// with centre of each region
int currLevel;
} REGIONS;
  
typedef struct {
// Represents the position coordinates of each block
} PLACEMENT;

void globalPlace();  
  
double WireLength();

void initPartitions();

void initMatrices(int level);
void setupCd_x();
void setupD();
void setupB();
void setupZ();
void setupu();
void setupx();
void setupc();

void destroyMatrices();
void destroyPartitions();

double func(const gsl_vector *v, void *params);
void dfunc(const gsl_vector *v, void *params, gsl_vector *df);
void fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df);

#endif
