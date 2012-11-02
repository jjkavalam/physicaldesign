void setupC() {
	blas_sparse_matrix C;
	C = BLAS_duscr_begin(movableNodes, movableNodes);
//	for (i = 0; i < ; i++)
//		BLAS_duscr_insert_entry(A, val[i], indx[i], jndx[i]);

	for (int n = 0; n < numNets; n++) {
		int nPerNode = netlistIndex[n+1] - netlistIndex[n];
		double e = 2.0/nPerNode;
		for (int t = 0; t < nPerNode; t++) {
			mu = netlist[netlistIndex[n] + t];
				
		}
	}

}

void initMatrices(int level) {
	setupC();
	setupD();
	setupB();
	setupZ();
	setupd_x()
	setupu();
	setupx();
	setupc();
}

double func(const gsl_vector *v, void *params) {
	// 0.5*x_T*Z_T*C*Zx + c_T*x
	
	double *x;
	int size = movableNodes - nRegions;
	
	x = (double*) malloc(sizeof(double)*size);

	for (int i = 0; i < size; i++)
		x[i] = gsl_vector_get(v, i);

	// use blas to multiply
	BLAS_dusmv(blas_no_trans, 1, Z, x, 1, tmp1, 1);
	BLAS_dusmv(blas_no_trans, 1, C, tmp1, 1, tmp2, 1);
	BLAS_dusmv(blas_trans, 0.5, Z, tmp2, 1, tmp3, 1);
	for (int i = 0; i < size; i++) {
		tmp += tmp3[i]*x[i];
		cost_x = tmp + c[i]*x[i];
	}

	return cost_x;
}

void dfunc(const gsl_vector *v, void *params, gsl_vector *df) {
	// Z_T*C*Z*x + c

	double *x;
	int size = movableNodes - nRegions;
	
	x = (double*) malloc(sizeof(double)*size);
	
	for (int i = 0; i < size; i++) 
		x[i] = gsl_vector_get(v, i);

	// TODO: compare x with that in func to reuse Z_T*C*Z*x 
	BLAS_dusmv(blas_no_trans, 1, Z, x, 1, tmp1, 1);
	BLAS_dusmv(blas_no_trans, 1, C, tmp1, 1, tmp2, 1);
	BLAS_dusmv(blas_trans, 1, Z, tmp2, 1, tmp3, 1);
		
	for (int i = 0; i < size; i++) 
		gsl_vector_set(df, i, tmp3[i] + c[i]);
}

void fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df) {
	*f = func(v, params);
	dfunc(v, params, df);
}
