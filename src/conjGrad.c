#include <stdio.h>
#include <gsl/gsl_multimin.h>

#define MAX_ITER 100

double func(const gsl_vector *v, void *params) {
	double x, y;
	double *p = (double *)params;

	x = gsl_vector_get(v, 0);
	y = gsl_vector_get(v, 1);

	return (x*x + p[0]*p[0]/(x*x) + 2*p[0] + y*y + p[1]*p[1]/(y*y) + 2*p[1]);
}

void dfunc(const gsl_vector *v, void *params, gsl_vector *df) {
	double x, y;
	double *p = (double *)params;

	x = gsl_vector_get(v, 0);
	y = gsl_vector_get(v, 1);

	gsl_vector_set(df, 0, 2*x - 2*p[0]*p[0]/(x*x*x));
	gsl_vector_set(df, 1, 2*y - 2*p[1]*p[1]/(y*y*y));
}

void fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df) {
	*f = func(v, params);
	dfunc(v, params, df);
}

main() {
	size_t iter = 0;
	int status;

	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	double p[3] = {1.0, 1.0, 0.0};
	gsl_vector *v;
	gsl_multimin_function_fdf function;

	// n = number of variables
	function.n = 2;
	function.f = &func;
	function.df = &dfunc;
	function.fdf = &fdf;
 	function.params = (void *)p;

	v = gsl_vector_alloc(2);
	gsl_vector_set(v, 0, 1000);
	gsl_vector_set(v, 1, 50);
	
	// Other option is _fr
	T = gsl_multimin_fdfminimizer_conjugate_pr;
	s = gsl_multimin_fdfminimizer_alloc(T, 2);

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
