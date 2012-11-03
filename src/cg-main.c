#include "placer.h"

#define MAX_ITER 100

int nRegions, level;	

main() {
	nRegions = 1;
	level = 0;

	initMatrices(0);	
	
	/* CG */
	size_t iter = 0;
	int status;

	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	double p[1] = {0.0};
	gsl_vector *v;
	gsl_multimin_function_fdf function;

	// n = number of variables
	function.n = 2;
	function.f = &func;
	function.df = &dfunc;
	function.fdf = &fdf;
 	function.params = (void *)p;
	
	int size = movableNodes - nRegions;
	v = gsl_vector_alloc(size);

	// Assume the .nodes file puts movableNodes before terminals
	for (int i = 0; i < size; i++)
		gsl_vector_set(v, i, x[i]);
	
	// Other option is _fr
	T = gsl_multimin_fdfminimizer_conjugate_pr;
	s = gsl_multimin_fdfminimizer_alloc(T, size);

	// The last two arguments are step-size and tolerance
	gsl_multimin_fdfminimizer_set(s, &function, v, 0.0001, 0.01);
	do {
		status = gsl_multimin_fdfminimizer_iterate(s);
		if (status)
			break;

		status = gsl_multimin_test_gradient(s->gradient, 0.01);
		printf("%d: %f %f %f\n", iter, gsl_vector_get(s->x, 0),	gsl_vector_get(s->x, 1), s->f);
		if (status == GSL_SUCCESS) {
			printf("%d: %f %f %f\n", iter, gsl_vector_get(s->x, 0),	gsl_vector_get(s->x, 1), s->f);
		}
	} while (status == GSL_CONTINUE && iter++ < MAX_ITER);

	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(v);
}
