#ifndef PLACER_H
#define PLACER_H

#include <gsl/gsl_multimin.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "parser/bookshelf_IO.h"
#include "../utils/spblas_0_8/csr_double.h"

extern double *C, *D, *B, *Z, *d_x, *u, *x, *c;
extern int nRegions;

void initMatrices(int level);
void setupC();
void setupD();
void setupB();
void setupZ();
void setupd_x();
void setupu();
void setupx();
void setupc();

double func(const gsl_vector *v, void *params);
void dfunc(const gsl_vector *v, void *params, gsl_vector *df);
void fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df);

#endif
