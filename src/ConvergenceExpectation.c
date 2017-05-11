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
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Utils.h"
#include "ConvergenceExpectation.h"
#include "ColumnLikelihood.h"
#include "NLOpt.h"
#include "EvolutionModelProtStored.h"

typedef struct CONVERGENCE_EXPECTATION {
	double **B;
	double *C;
	int size;
} TypeConvergenceExpectation;

static TypeConvergenceExpectation *pendingConvergenceExpectation(TypeEvolutionModel *model, int AA, double t, int node, TypeConvergenceExpectation *ce);
static TypeConvergenceExpectation *fillConvergenceExpectation(TypeTree *tree, int *character, TypeEvolutionModel *model, int *column, int AA, int node);
static void freeConvergenceExpectation(TypeConvergenceExpectation *ce);
static TypeConvergenceExpectation* initConvergenceExpectation(int cardinal, int size);
static double toOpt(unsigned k, const double *x, double *grad, void *data);

typedef struct OPT_DATA { 
	TypeEvolutionModel *model;
	TypeTree *tree; 
	int *character, *column, symbol;
	double ce;
} TypeOptData;

double toOpt(unsigned k, const double *x, double *grad, void *data) {
	double *timeSave, res;
	int n, a[2]={((TypeOptData *)data)->symbol, END_INT};
	timeSave = ((TypeOptData *)data)->tree->time;
	((TypeOptData *)data)->tree->time = (double*) malloc(((TypeOptData *)data)->tree->size*sizeof(double));
	for(n=0; n<((TypeOptData *)data)->tree->size; n++)
		((TypeOptData *)data)->tree->time[n] = x[0]*timeSave[n];
	res = fabs(getConvergenceIndexList(((TypeOptData *)data)->tree, ((TypeOptData *)data)->character, ((TypeOptData *)data)->model, ((TypeOptData *)data)->column, a)-((TypeOptData *)data)->ce);
	free((void*)((TypeOptData *)data)->tree->time);
	((TypeOptData *)data)->tree->time = timeSave;
   return res;
}

double calibRate(TypeTree *tree, int *character, TypeEvolutionModel *model, void *randG, int a, double ce) {
	TypeOptData data;
	TypeNLOptOption *option=getNLOption();
	nlopt_opt opt;
	double min, x=1., rate=1.;
	int status, n, t;
	data.model = model;
	data.tree = tree;
	data.character = character;
	data.symbol = a;
	data.ce = ce;
	data.column = (int*) malloc(tree->size*sizeof(int));
	for(n=0; n<tree->size; n++)
		data.column[n] = data.symbol;
	opt = nlopt_create(NLOPT_ALGO, 1); 
	nlopt_set_lower_bounds1(opt, log(0.75)/(getMinSteady(model)*getMaxLength(tree)));
	nlopt_set_upper_bounds1(opt, log(0.25)/(getMinSteady(model)*getMinLength(tree)));
	nlopt_set_min_objective(opt, toOpt, &data);
	nlopt_set_xtol_abs1(opt, option->tolOptim);
	nlopt_set_maxeval(opt, option->maxIter);
	min = HUGE_VAL;
	for(t=0; t<option->trials; t++) {
		double tmp;
		fillRandStartNLOpt(opt, randG, &x);
		if(((status = nlopt_optimize(opt, &x, &tmp))>=0)) {
			if(tmp < min) {
				rate = x;
				min = tmp;
			}
		 }
	}
	if(min == HUGE_VAL) {
		fprintf(stderr, "Numerical optimisation failed!\n");
		exit(1);
	}
	free((void*)data.column);
	nlopt_destroy(opt);
//	printf("res %.2le/%.2le\n", min, ce);
	return rate;	
}

int *getSymbols(TypeTree *tree, int *character, int cardinal, int *column) {
	int a, i, j, nConv, *al, *conv, *occ;
	nConv = 0;
	for(i=0; tree->node[i].child == NOSUCH; i++)
		if(character[i])
			nConv++;
	conv = (int*) malloc((nConv+1)*sizeof(int));
	for(i=0, j=0; tree->node[i].child == NOSUCH; i++)
		if(character[i])
			conv[j++] = i;
	occ = (int*) malloc(cardinal*sizeof(int));
	for(i=0, j=0; i<cardinal; i++)
		occ[i] = 0;
	al = (int*) malloc((nConv+1)*sizeof(int));
	for(i=0, a=0; i<nConv; i++) {
		occ[column[conv[i]]]++;
		if(occ[column[conv[i]]] == 1)
			al[a++] = column[conv[i]];
	}
	al[a] = END_INT;
	free((void*)occ);
	free((void*)conv);
	return al;
}

double getConvergenceIndex(TypeTree *tree, int *character, TypeEvolutionModel *model, int *column) {
	int *al;
	double ci;
	al = getSymbols(tree, character, model->cardinal, column);
	ci = getConvergenceIndexList(tree, character, model, column, al);
	free((void*)al);
	return ci;
}

double getConvergenceIndexList(TypeTree *tree, int *character, TypeEvolutionModel *model, int *column, int *al) {
    int a;
    double ci=0.;
	for(a=0; al[a] != END_INT; a++) {
		double tmp = getConvergenceExpectation(tree, character, model, column, al[a]);
		if(tmp>ci)
			ci = tmp;
	}
	return ci/getLikelihoodColumn(tree, model, column);
}

double getConvergenceExpectation(TypeTree *tree, int *character, TypeEvolutionModel *model, int *column, int AA) {
	int k, z;
	double res = 0.;
	TypeConvergenceExpectation *p = fillConvergenceExpectation(tree, character, model, column, AA, tree->root);  
	for(k=0; k<=p->size ;k++) {
		double sum = 0.;
		for(z=0; z<model->cardinal; z++)
			sum += model->getInitial(model, z)*p->B[k][z];
		res += k*sum+(k+1)*model->getInitial(model, AA)*p->C[k];
	}
	freeConvergenceExpectation(p);
	return res;
}

TypeConvergenceExpectation *pendingConvergenceExpectation(TypeEvolutionModel *model, int AA, double t, int node, TypeConvergenceExpectation *ce) {	
	int y, z, k;
	TypeConvergenceExpectation *pce = initConvergenceExpectation(model->cardinal, ce->size);
	model->setTime(model, t, &node);
	double steadAA = model->getSteady(model, AA);
	for(y=0; y<AA; y++) {
		double transyAA = model->getTransition(model, y, AA);
		pce->B[0][y] = 0.;
		for(k=1; k<=ce->size; k++)
			pce->B[k][y] = (ce->C[k-1])*transyAA;
		for(z=0; z<model->cardinal; z++) {
			double transyz = model->getTransition(model, y, z);
			pce->B[0][y] += ce->B[0][z]*transyz;
			for(k=1; k<=ce->size; k++)
				pce->B[k][y] += ce->B[k][z]*transyz;
		}
	}
	pce->B[0][AA] = 0.;
	for(k=1; k<=ce->size; k++)
		pce->B[k][AA] = ce->C[k-1]*(model->getTransition(model, AA, AA)-steadAA);
	for(z=0; z<model->cardinal; z++) {
		double transAAz = model->getTransition(model, AA, z);
		pce->B[0][AA] += ce->B[0][z]*transAAz;
		for(k=1; k<=ce->size; k++)
			pce->B[k][AA] += ce->B[k][z]*transAAz;
	}
	for(y=AA+1; y<model->cardinal; y++) {
		double transyAA = model->getTransition(model, y, AA);
		pce->B[0][y] = 0.;
		for(k=1; k<=ce->size; k++)
			pce->B[k][y] = (ce->C[k-1])*transyAA;
		for(z=0; z<model->cardinal; z++) {
			double transyz = model->getTransition(model, y, z);
			pce->B[0][y] += ce->B[0][z]*transyz;
			for(k=1; k<=ce->size; k++)
				pce->B[k][y] += ce->B[k][z]*transyz;
		}
	}
	for(k=0; k<=ce->size; k++)
		pce->C[k] = ce->C[k]*steadAA;
	return pce;
}

TypeConvergenceExpectation *fillConvergenceExpectation(TypeTree *tree, int *character, TypeEvolutionModel *model, int *column, int AA, int node) {
	int k, y;
	TypeConvergenceExpectation *ret;
	if(tree->node[node].child == NOSUCH) {
		ret = initConvergenceExpectation(model->cardinal, character[node]?1:0);  
		if(column[node] ==  AA && character[node])
			ret->C[0] = 1.;
		else {
			if(column[node]<0 || column[node]>=model->cardinal) {
				int a = 0;
				for(a=0; a<model->cardinal; a++)
					ret->B[0][a] = 1.;
			} else
				ret->B[0][column[node]] = 1.;
		}
	} else {
		int d=tree->node[node].child,  g=tree->node[tree->node[node].child].sibling, k1, k2;
		TypeConvergenceExpectation *fd, *fg, *bd, *bg;
		fd = fillConvergenceExpectation(tree, character, model, column, AA, d);
		fg = fillConvergenceExpectation(tree, character, model, column, AA, g);
		bd = pendingConvergenceExpectation(model, AA, tree->time[d], d, fd);
		bg = pendingConvergenceExpectation(model, AA, tree->time[g], g, fg);
		ret = initConvergenceExpectation(model->cardinal, fd->size+fg->size);
		for(k=0; k<=ret->size; k++)
			for(y=0; y<model->cardinal; y++) 
				ret->B[k][y] = 0.;
		for(k1=0; k1<=fd->size; k1++)
			for(k2=0; k2<=fg->size; k2++) {
				for(y=0; y<AA ; y++) 
					ret->B[k1+k2][y] += bd->B[k1][y]*bg->B[k2][y];
				ret->B[k1+k2][AA] += bd->B[k1][AA]*bg->B[k2][AA];
				for(y=AA+1; y<model->cardinal; y++) 
					ret->B[k1+k2][y] += bd->B[k1][y]*bg->B[k2][y];
			}
		for(k=0; k<=ret->size; k++)
			ret->C[k] = 0.;
		for(k1=0; k1<=fd->size; k1++)
			for(k2=0; k2<=fg->size; k2++)
				ret->C[k1+k2] += bd->C[k1]*bg->C[k2]+bd->B[k1][AA]*bg->C[k2]+bd->C[k1]*bg->B[k2][AA];
		freeConvergenceExpectation(fd);
		freeConvergenceExpectation(fg);
		freeConvergenceExpectation(bd);
		freeConvergenceExpectation(bg);
    }
	return ret;
}

void freeConvergenceExpectation(TypeConvergenceExpectation *ce) {
	if (ce == NULL)
		return;
	if(ce->B != NULL) {
		int i;
		for(i=0; i<=ce->size; i++)
			if (ce->B[i] != NULL)
				free((void*)ce->B[i]);
		free(ce->B);
	}
	if(ce->C != NULL)
		free(ce->C);
	free(ce);
}

TypeConvergenceExpectation* initConvergenceExpectation(int cardinal, int size) {
	int i, j;
	TypeConvergenceExpectation *ret = (TypeConvergenceExpectation *) malloc(sizeof(TypeConvergenceExpectation));
	ret->size = size;
	ret->B = (double **) malloc((ret->size+1)*sizeof(double *));
	ret->C = (double *) malloc((ret->size+1)*sizeof(double));
	for(i=0; i<=ret->size; i++) {
		ret->B[i] = (double *) malloc(cardinal*sizeof(double));
		for(j=0; j<cardinal; j++)
			ret->B[i][j] = 0.;
		ret->C[i] = 0.;
	}
	return ret;
}
