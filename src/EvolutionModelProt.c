/*
    'msd' detects molecular signatures of phenotypic convergence / 
    'enr' computes GO terms enrichments of a list of genes / 

    Copyright (C) 2017 Gilles DIDIER

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex.h>
#include "EvolutionModelProt.h"

static void setInitGeneratorRef(double *mat, double *pi);
static void setGenerator(TypeProtData *data);
static void setMatricesRef();
static double getTransitionProt(TypeEvolutionModel *model, int x, int y);
static double getLogTransitionProt(TypeEvolutionModel *model, int x, int y);
static double getSteadyProt(TypeEvolutionModel *model, int x);
static double getLogSteadyProt(TypeEvolutionModel *model, int x);
static double getInitialProt(TypeEvolutionModel *model, int x);
static double getLogInitialProt(TypeEvolutionModel *model, int x);
static int setTimeProt(TypeEvolutionModel *model, double t, void *user);
static void initRef();

double matWAG[400] = {-1.894441e+01, 5.515710e-01, 5.098480e-01, 7.389980e-01, 1.027040e+00, 9.085980e-01, 1.582850e+00, 1.416720e+00, 3.169540e-01, 1.933350e-01, 3.979150e-01, 9.062650e-01, 8.934960e-01, 2.104940e-01, 1.438550e+00, 3.370790e+00, 2.121110e+00, 1.131330e-01, 2.407350e-01, 2.006010e+00,
5.515710e-01, -1.913622e+01, 6.353460e-01, 1.473040e-01, 5.281910e-01, 3.035500e+00, 4.391570e-01, 5.846650e-01, 2.137150e+00, 1.869790e-01, 4.976710e-01, 5.351420e+00, 6.831620e-01, 1.027110e-01, 6.794890e-01, 1.224190e+00, 5.544130e-01, 1.163920e+00, 3.815330e-01, 2.518490e-01,
5.098480e-01, 6.353460e-01, -2.595825e+01, 5.429420e+00, 2.652560e-01, 1.543640e+00, 9.471980e-01, 1.125560e+00, 3.956290e+00, 5.542360e-01, 1.315280e-01, 3.012010e+00, 1.982210e-01, 9.616210e-02, 1.950810e-01, 3.974230e+00, 2.030060e+00, 7.191670e-02, 1.086000e+00, 1.962460e-01,
7.389980e-01, 1.473040e-01, 5.429420e+00, -1.816622e+01, 3.029490e-02, 6.167830e-01, 6.174160e+00, 8.655840e-01, 9.306760e-01, 3.943700e-02, 8.480470e-02, 4.798550e-01, 1.037540e-01, 4.673040e-02, 4.239840e-01, 1.071760e+00, 3.748660e-01, 1.297670e-01, 3.257110e-01, 1.523350e-01,
1.027040e+00, 5.281910e-01, 2.652560e-01, 3.029490e-02, -8.236647e+00, 9.881790e-02, 2.135200e-02, 3.066740e-01, 2.489720e-01, 1.701350e-01, 3.842870e-01, 7.403390e-02, 3.904820e-01, 3.980200e-01, 1.094040e-01, 1.407660e+00, 5.129840e-01, 7.170700e-01, 5.438330e-01, 1.002140e+00,
9.085980e-01, 3.035500e+00, 1.543640e+00, 6.167830e-01, 9.881790e-02, -2.638536e+01, 5.469470e+00, 3.300520e-01, 4.294110e+00, 1.139170e-01, 8.694890e-01, 3.894900e+00, 1.545260e+00, 9.992080e-02, 9.333720e-01, 1.028870e+00, 8.579280e-01, 2.157370e-01, 2.277100e-01, 3.012810e-01,
1.582850e+00, 4.391570e-01, 9.471980e-01, 6.174160e+00, 2.135200e-02, 5.469470e+00, -2.218592e+01, 5.677170e-01, 5.700250e-01, 1.273950e-01, 1.542630e-01, 2.584430e+00, 3.151240e-01, 8.113390e-02, 6.823550e-01, 7.049390e-01, 8.227650e-01, 1.565570e-01, 1.963030e-01, 5.887310e-01,
1.416720e+00, 5.846650e-01, 1.125560e+00, 8.655840e-01, 3.066740e-01, 3.300520e-01, 5.677170e-01, -8.574782e+00, 2.494100e-01, 3.045010e-02, 6.130370e-02, 3.735580e-01, 1.741000e-01, 4.993100e-02, 2.435700e-01, 1.341820e+00, 2.258330e-01, 3.369830e-01, 1.036040e-01, 1.872470e-01,
3.169540e-01, 2.137150e+00, 3.956290e+00, 9.306760e-01, 2.489720e-01, 4.294110e+00, 5.700250e-01, 2.494100e-01, -2.147922e+01, 1.381900e-01, 4.994620e-01, 8.904320e-01, 4.041410e-01, 6.793710e-01, 6.961980e-01, 7.401690e-01, 4.733070e-01, 2.625690e-01, 3.873440e+00, 1.183580e-01,
1.933350e-01, 1.869790e-01, 5.542360e-01, 3.943700e-02, 1.701350e-01, 1.139170e-01, 1.273950e-01, 3.045010e-02, 1.381900e-01, -2.069729e+01, 3.170970e+00, 3.238320e-01, 4.257460e+00, 1.059470e+00, 9.992880e-02, 3.194400e-01, 1.458160e+00, 2.124830e-01, 4.201700e-01, 7.821300e+00,
3.979150e-01, 4.976710e-01, 1.315280e-01, 8.480470e-02, 3.842870e-01, 8.694890e-01, 1.542630e-01, 6.130370e-02, 4.994620e-01, 3.170970e+00, -1.742991e+01, 2.575550e-01, 4.854020e+00, 2.115170e+00, 4.158440e-01, 3.447390e-01, 3.266220e-01, 6.653090e-01, 3.986180e-01, 1.800340e+00,
9.062650e-01, 5.351420e+00, 3.012010e+00, 4.798550e-01, 7.403390e-02, 3.894900e+00, 2.584430e+00, 3.735580e-01, 8.904320e-01, 3.238320e-01, 2.575550e-01, -2.265861e+01, 9.342760e-01, 8.883600e-02, 5.568960e-01, 9.671300e-01, 1.386980e+00, 1.375050e-01, 1.332640e-01, 3.054340e-01,
8.934960e-01, 6.831620e-01, 1.982210e-01, 1.037540e-01, 3.904820e-01, 1.545260e+00, 3.151240e-01, 1.741000e-01, 4.041410e-01, 4.257460e+00, 4.854020e+00, 9.342760e-01, -2.112807e+01, 1.190630e+00, 1.713290e-01, 4.939050e-01, 1.516120e+00, 5.157060e-01, 4.284370e-01, 2.058450e+00,
2.104940e-01, 1.027110e-01, 9.616210e-02, 4.673040e-02, 3.980200e-01, 9.992080e-02, 8.113390e-02, 4.993100e-02, 6.793710e-01, 1.059470e+00, 2.115170e+00, 8.883600e-02, 1.190630e+00, -1.573167e+01, 1.614440e-01, 5.459310e-01, 1.719030e-01, 1.529640e+00, 6.454280e+00, 6.498920e-01,
1.438550e+00, 6.794890e-01, 1.950810e-01, 4.239840e-01, 1.094040e-01, 9.333720e-01, 6.823550e-01, 2.435700e-01, 6.961980e-01, 9.992880e-02, 4.158440e-01, 5.568960e-01, 1.713290e-01, 1.614440e-01, -9.886447e+00, 1.613280e+00, 7.953840e-01, 1.394050e-01, 2.160460e-01, 3.148870e-01,
3.370790e+00, 1.224190e+00, 3.974230e+00, 1.071760e+00, 1.407660e+00, 1.028870e+00, 7.049390e-01, 1.341820e+00, 7.401690e-01, 3.194400e-01, 3.447390e-01, 9.671300e-01, 4.939050e-01, 5.459310e-01, 1.613280e+00, -2.507035e+01, 4.378020e+00, 5.237420e-01, 7.869930e-01, 2.327390e-01,
2.121110e+00, 5.544130e-01, 2.030060e+00, 3.748660e-01, 5.129840e-01, 8.579280e-01, 8.227650e-01, 2.258330e-01, 4.733070e-01, 1.458160e+00, 3.266220e-01, 1.386980e+00, 1.516120e+00, 1.719030e-01, 7.953840e-01, 4.378020e+00, -1.979670e+01, 1.108640e-01, 2.911480e-01, 1.388230e+00,
1.131330e-01, 1.163920e+00, 7.191670e-02, 1.297670e-01, 7.170700e-01, 2.157370e-01, 1.565570e-01, 3.369830e-01, 2.625690e-01, 2.124830e-01, 6.653090e-01, 1.375050e-01, 5.157060e-01, 1.529640e+00, 1.394050e-01, 5.237420e-01, 1.108640e-01, -9.853066e+00, 2.485390e+00, 3.653690e-01,
2.407350e-01, 3.815330e-01, 1.086000e+00, 3.257110e-01, 5.438330e-01, 2.277100e-01, 1.963030e-01, 1.036040e-01, 3.873440e+00, 4.201700e-01, 3.986180e-01, 1.332640e-01, 4.284370e-01, 6.454280e+00, 2.160460e-01, 7.869930e-01, 2.911480e-01, 2.485390e+00, -1.890794e+01, 3.147300e-01,
2.006010e+00, 2.518490e-01, 1.962460e-01, 1.523350e-01, 1.002140e+00, 3.012810e-01, 5.887310e-01, 1.872470e-01, 1.183580e-01, 7.821300e+00, 1.800340e+00, 3.054340e-01, 2.058450e+00, 6.498920e-01, 3.148870e-01, 2.327390e-01, 1.388230e+00, 3.653690e-01, 3.147300e-01, -2.005557e+01};

double piWAG[20] = {8.662790e-02, 4.397200e-02, 3.908940e-02, 5.704510e-02, 1.930780e-02, 3.672810e-02, 5.805890e-02, 8.325180e-02, 2.443130e-02, 4.846600e-02, 8.620900e-02, 6.202860e-02, 1.950270e-02, 3.843190e-02, 4.576310e-02, 6.951790e-02, 6.101270e-02, 1.438590e-02, 3.527420e-02, 7.089560e-02};

typedef struct TMP_PROT_DATA {
    double initial[CARDINAL_PROT], Q[CARDINAL_PROT2], D[CARDINAL_PROT], P[CARDINAL_PROT2], Pi[CARDINAL_PROT2];
    int done;
} TypeTmpProtData;

TypeTmpProtData ref = {.done = 0};

void setInitGeneratorRef(double *mat, double *pi) {
	int i, j;
	for(i=0; i<CARDINAL_PROT; i++)
		ref.initial[i] = pi[i];
	for(i=0; i<CARDINAL_PROT; i++) {
		ref.Q[i*CARDINAL_PROT+i] = 0.;
		for(j=0; j<i; j++) {
			ref.Q[i*CARDINAL_PROT+j] = mat[i*CARDINAL_PROT+j]*pi[j];
			ref.Q[i*CARDINAL_PROT+i] -= ref.Q[i*CARDINAL_PROT+j];
		}
		for(j=i+1; j<CARDINAL_PROT; j++) {
			ref.Q[i*CARDINAL_PROT+j] = mat[i*CARDINAL_PROT+j]*pi[j];
			ref.Q[i*CARDINAL_PROT+i] -= ref.Q[i*CARDINAL_PROT+j];
		}
	}
}

void setGenerator(TypeProtData *data) {
	data->initial = ref.initial;
	data->Q = ref.Q;
	data->D = ref.D;
	data->P = ref.P;
	data->Pi = ref.Pi;
}

void setMatricesRef() {
	int i, j, ind, s;
	gsl_matrix *mtmp = gsl_matrix_alloc(CARDINAL_PROT, CARDINAL_PROT);
	for(i=0, ind = 0; i<CARDINAL_PROT; i++)
		for(j=0; j<CARDINAL_PROT; j++)
			gsl_matrix_set(mtmp, i, j, ref.Q[ind++]);
	gsl_vector_complex *eval = gsl_vector_complex_alloc(CARDINAL_PROT);
	gsl_matrix_complex *evec = gsl_matrix_complex_alloc(CARDINAL_PROT, CARDINAL_PROT);
	gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(CARDINAL_PROT);
	gsl_eigen_nonsymmv(mtmp, eval, evec, w);
	gsl_eigen_nonsymmv_free(w);
	for(i=0; i<CARDINAL_PROT; i++)
		ref.D[i] = gsl_vector_complex_get(eval, i).dat[0];
	gsl_vector_complex_free(eval);
	for(i=0, ind=0; i<CARDINAL_PROT; i++) {
		for(j=0; j<CARDINAL_PROT; j++) {
			ref.P[ind] = gsl_matrix_complex_get(evec, i, j).dat[0];
			gsl_matrix_set(mtmp, i, j, ref.P[ind]);
			ind++;
		}
	}
	gsl_matrix_complex_free(evec);
	gsl_matrix *inverse = gsl_matrix_alloc(CARDINAL_PROT, CARDINAL_PROT);
	gsl_permutation *perm = gsl_permutation_alloc(CARDINAL_PROT);
	gsl_linalg_LU_decomp(mtmp, perm, &s);
	gsl_linalg_LU_invert(mtmp, perm, inverse);
	for(i=0, ind=0; i<CARDINAL_PROT; i++)
		for(j=0; j<CARDINAL_PROT; j++)
			ref.Pi[ind++] = gsl_matrix_get(inverse, i, j);
	gsl_matrix_free (mtmp);
	gsl_matrix_free (inverse);
	gsl_permutation_free(perm);
}

void expMatrix(void *data, double t, double *res) {
	double vtmp[CARDINAL_PROT], mtmp[CARDINAL_PROT2];
	int i, j, ind=0;
	
	for(i=0; i<CARDINAL_PROT; i++)
		vtmp[i] =  exp(((TypeProtData*)data)->D[i]*t);
	for(i=0; i<CARDINAL_PROT; i++)
		for(j=0; j<CARDINAL_PROT; j++) {
			mtmp[ind] = vtmp[j]*((TypeProtData*)data)->P[ind];
			ind++;
		}
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, CARDINAL_PROT, CARDINAL_PROT, CARDINAL_PROT, 1.0, mtmp, CARDINAL_PROT, ((TypeProtData*)data)->Pi, CARDINAL_PROT, 0.0, res, CARDINAL_PROT);
/*lda taille physique des lignes de A
 * M nbre de lignes de A
 * N nbre de colonnes de B
 * K nbre de colonnes de A et de lignes de B
cblas_dgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double * A, const int lda, const double * B, const int ldb, const double beta, double * C, const int ldc)
calcule C<- alphaAB+betaC*/
}

double getTransitionProt(TypeEvolutionModel *model, int x, int y) {
	if(x<CARDINAL_PROT && y<CARDINAL_PROT)
		return ((TypeProtData*)model->data)->trans[x*CARDINAL_PROT+y];
	return 1.;
}

double getLogTransitionProt(TypeEvolutionModel *model, int x, int y) {
	if(x<CARDINAL_PROT && y<CARDINAL_PROT)
		return log(((TypeProtData*)model->data)->trans[x*CARDINAL_PROT+y]);
	return 0.;
}

double getSteadyProt(TypeEvolutionModel *model, int x) {
	return exp(((TypeProtData*)model->data)->Q[x*CARDINAL_PROT+x]*model->t);
}

double getLogSteadyProt(TypeEvolutionModel *model, int x) {
	return ((TypeProtData*)model->data)->Q[x*CARDINAL_PROT+x]*model->t;
}

double getInitialProt(TypeEvolutionModel *model, int x) {
	return ((TypeProtData*)model->data)->initial[x];
}

double getLogInitialProt(TypeEvolutionModel *model, int x) {
	return log(((TypeProtData*)model->data)->initial[x]);
}

int setTimeProt(TypeEvolutionModel *model, double t, void *user) {
	model->t = t;
	expMatrix(model->data, t, ((TypeProtData*)model->data)->trans);
	return 1.;
}

void freeModelProt(TypeEvolutionModel *model) {
	free((void*)((TypeProtData*)model->data)->trans);
	free(model->data);
	free((void*)model);
}

void initRef() {
	setInitGeneratorRef(matWAG, piWAG);
	setMatricesRef();
	ref.done = 1;
}

TypeEvolutionModel *getEvolutionModelProt() {
	TypeEvolutionModel *model;
	if(ref.done == 0)
		initRef();
	model = (TypeEvolutionModel *) malloc(sizeof(TypeEvolutionModel));
	model->cardinal = CARDINAL_PROT;
	model->data = malloc(sizeof(TypeProtData));
	((TypeProtData*)model->data)->trans = (double*) malloc(CARDINAL_PROT2*sizeof(double));
	((TypeProtData*)model->data)->user = NULL;
	setGenerator(model->data);
	setTimeProt(model, 1., NULL);
	model->freeModel = freeModelProt;
	model->setTime = setTimeProt;
	model->getTransition = getTransitionProt;
	model->getSteady = getSteadyProt;
	model->getInitial = getInitialProt;
	model->getLogTransition = getLogTransitionProt;
	model->getLogSteady = getLogSteadyProt;
	model->getLogInitial = getLogInitialProt;
	return model;
}
