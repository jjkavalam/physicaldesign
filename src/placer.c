#include "placer.h"
#include "compute.h"

// cT will be computed before each CG run
// and it will be available as a global variable to func, dfunc  
VECT cT;

// update cT
void updatecT();


PLACEMENT globalPlace(Partition, Regions){
    // local variables
    int i, j, flag;
    // Find the index of the largest area block in each partition
    int* depBlocks;
    int* isDep;

    // Setup the matrices based on the new partition
    setupDBxi(Partition,depBlocks);

    // update cT which is fixed through the run
    updatecT();

    // setup the initial positions of the block to centre of each region
    VECT xi;        
    for (i = 0;i<movableNodes;i++){
        xi.val[i] = Regions.centres[Partition[i]];
    }

    // setup the Centre of Gravity constraits vector u
    VECT u = Regions.u;
        
    // run the CG loop
    runCG();

    // compute the positions of the dependent blocks as well
    VECT v1 = computeBv(xi);    
    VECT v2 = computeinvDv(v1);
    // negate vector v2
    VECT v3 = vectscale(v2, -1);
    VECT v4 = computeinvDv(u);
    VECT xd = vectsum(v3,v4);
        
    // prepare the struct for return
    // use the same ordering as given by the partitioner
    double* xCoords = (double*) malloc(sizeof(double)*movableNodes);
    for (i=0,j=0;i<movableNodes;i++){
        if (isDep[i]){
            xCoords[i] = xd[j];
            j++;
        }
        else
        {
            xCoords[i] = xi[i-j];
        }
    }
    PLACEMENT P;
    P.xCoords = xCoords;

    return P;

}

void destroyMatrices(){

    // Destroy BLAS data types
    BLAS_usds(C);
    //BLAS_usds(Z); 
    
    destroyxy();
    destroydxdy();   
}

void freeVECT(VECT v){
    free(v.val);    
}

void initPartitions(){
    int i;
    Partition = (int*)malloc(sizeof(int)*movableNodes);
    for (i=0;i<movableNodes;i++) 
        Partition[i] = 0;
        
}
void destroyPartitions(){
    free(Partition);    
}

PLACEMENT initialPlacement(){  
    // setup xy
    x = (double*)malloc(sizeof(double)*(movableNodes+numTerminals));
    memset((void*)x,0,sizeof(double)*(movableNodes+numTerminals));

    y = (double*)malloc(sizeof(double)*(movableNodes+numTerminals));
    memset((void*)y,0,sizeof(double)*(movableNodes+numTerminals));

    // Level 0 initialization to the values in the .pl file        
    int i;
    for(i=0;i<(movableNodes+numTerminals);i++){
        x[i] = xCellCoord[i+1];
        y[i] = yCellCoord[i+1];
    }    

    // setup the fixed matrices  
    setupCdxdy(); 

}
destroyPlacement(PLACEMENT P){
    free(P.x);
    free(P.y);
}

double func(const gsl_vector *v, void *params) {

	// 0.5*x_T*Z_T*C*Zx + c_T*x
	// params is the number of variables
	
	double Cost, sum1, sum2;
	
	int size = params[0];
	
    VECT Vx, Vy;
    
    Vx.val = (double*) malloc(sizeof(double)*size);
    Vy.val = (double*) malloc(sizeof(double)*size);
    
	for (int i = 0; i < size; i++)
		x[i] = gsl_vector_get(v, i);

    VECT v6 = computeZTCZx(Vx);
    sum1 = dotxi(v6);
    sum2 = dotxi(cT);
    Cost = 0.5 * sum1 + sum2;    

    freeVECT(Vx);
    freeVECT(Vy);
    
    return Cost;
}

void dfunc(const gsl_vector *v, void *params, gsl_vector *df) {
	// Z_T*C*Z*x + c

	int size = params[0];
	
    VECT Vx, Vy;
    
    Vx.val = (double*) malloc(sizeof(double)*size);
    Vy.val = (double*) malloc(sizeof(double)*size);
    
	for (int i = 0; i < size; i++)
		x[i] = gsl_vector_get(v, i);
	
    // A possible approximation is to use the already computed 
    VECT g = vectsum(computeZTCZx(v),cT); 
		
	for (int i = 0; i < size; i++) 
		gsl_vector_set(df, i, g.val[i]);

    freeVECT(Vx);
    freeVECT(Vy);
    
    return;
		
}

void fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df) {
	*f = func(v, params);
	dfunc(v, params, df);
}

void runCG(){
	// Conjugate Gradient 
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
