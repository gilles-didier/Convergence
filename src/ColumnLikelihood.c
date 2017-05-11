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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "ColumnLikelihood.h"

static double *getConditionedLike(TypeTree *tree, TypeEvolutionModel *model, int *column, int node);
static void drawColumnRec(TypeTree *tree, TypeEvolutionModel *model, void *rand, int *column, int node);
static gsl_rng *getRandGen(void *RandG);

static gsl_rng *randGen = NULL;

void initializeDrawColumn() {
	randGen = gsl_rng_alloc(gsl_rng_taus);
}

void freeDrawColumn() {
	if(randGen != NULL)
		gsl_rng_free(randGen);
}

double *getGamma(double alpha, int nCat) {  
    int i;
    double *res = (double*) malloc((nCat+1)*sizeof(double));
	for(i=0; i<nCat; i++)
		res[i] = gsl_cdf_gamma_Pinv(((double)i+0.5)/((double)nCat), alpha, alpha);
	res[i] = END_DOUBLE;
	return res;
}

double getGammaLikelihoodColumn(TypeTree *tree, TypeEvolutionModel **model, double *gamma, int *column) {
	int g;
	double like=0.;
	for(g=0; gamma[g] != END_DOUBLE; g++)
		like += getLikelihoodColumn(tree, model[g], column);
	return like/((double)g);
}

double getLikelihoodColumn(TypeTree *tree, TypeEvolutionModel *model, int *column) {
	double *clikes, like=0.;
	int a;
	clikes = getConditionedLike(tree, model, column, tree->root);
	for(a=0; a<model->cardinal; a++)
		like += model->getInitial(model, a)*clikes[a];
	free((void*)clikes);
	return like;	
}

double *getConditionedLike(TypeTree *tree, TypeEvolutionModel *model, int *column, int node) {
	double *prob;
	int a, b;
	prob = (double*) malloc(model->cardinal*sizeof(double));
	if(tree->node[node].child == NOSUCH) {
		if(column[node]>=0 && column[node]<model->cardinal) {
			for(a=0; a<column[node]; a++)
				prob[a] = 0.;
			prob[column[node]] = 1.;
			for(a=column[node]+1; a<model->cardinal; a++)
				prob[a] = 0.;
		} else
			for(a=0; a<model->cardinal; a++)
				prob[a] = 1.;
	} else {
		double *pchild;
		int c = tree->node[node].child;
		pchild = getConditionedLike(tree, model, column, c);
		model->setTime(model, tree->time[c], &c);
		for(a=0; a<model->cardinal; a++) {
			prob[a] = 0.;
			for(b=0; b<model->cardinal; b++)
				prob[a] += model->getTransition(model, a, b)*pchild[b];
		}
		free((void*)pchild);
		for(c=tree->node[tree->node[node].child].sibling; c!=NOSUCH; c=tree->node[c].sibling) {
			pchild = getConditionedLike(tree, model, column, c);
			model->setTime(model, tree->time[c], &c);
			for(a=0; a<model->cardinal; a++) {
				double tmp = 0.;
				for(b=0; b<model->cardinal; b++)
					tmp += model->getTransition(model, a, b)*pchild[b];
				prob[a] *= tmp;
			}
			free((void*)pchild);
		}
	}
	return prob;	
}

void drawColumnRec(TypeTree *tree, TypeEvolutionModel *model, void *rand, int *column, int node) {
    if (tree->node[node].child != NOSUCH) {
		int c;
		for(c=tree->node[node].child; c!=NOSUCH; c=tree->node[c].sibling) {
			model->setTime(model, tree->time[c], &c);
			column[c] = drawTransitionFrom(model, rand, column[node]);
			drawColumnRec(tree, model, rand, column, c);
		}
	}
}

void drawColumn(TypeTree *tree, TypeEvolutionModel *model, void *randG, int *column) {	
	gsl_rng *rg = getRandGen(randG);
	column[tree->root] = drawInitial(model, randG);
	drawColumnRec(tree, model, rg, column, tree->root);
}

void drawColumnUnknown(TypeTree *tree, TypeEvolutionModel *model, double *unknown, void *randG, int *column) {
	int n;
	gsl_rng *rg = getRandGen(randG);
	drawColumn(tree, model, rg, column);
	for(n=0; tree->node[n].child==NOSUCH; n++)
		if(gsl_ran_bernoulli(rg, unknown[n]))
			column[n] = END_INT;
}

gsl_rng *getRandGen(void *randG) {
	gsl_rng *rg;
	if(randG != NULL) {
		rg = (gsl_rng *) randG;
	} else {
		if(randGen == NULL)
			initializeDrawColumn();
		rg = randGen;
	}
	return rg;
}
