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
#include <nlopt.h>
#include "AlignmentLikelihood.h"
#include "ColumnLikelihood.h"
#include "EvolutionModelProt.h"
#include "EvolutionModelProtStored.h"
#include "NLOpt.h"

typedef struct OPT_DATA { 
	TypeEvolutionModel *model;
	TypeTree *tree; 
	TypeAlignment *align;
	int nCat;
} TypeOptData;

static double toOpt(unsigned n, const double *x, double *grad, void *data);

double getGammaLogLikelihoodAlignment(TypeTree *tree, TypeEvolutionModel *model, double *gamma, TypeAlignment *align) {
	int g, nG, n, *column, pos;
	double logLike = 0.;
	TypeEvolutionModel **modelT;
	for(g=0; gamma[g] != END_DOUBLE; g++)
		;
	nG = g;
	modelT = (TypeEvolutionModel**) malloc(nG*sizeof(TypeEvolutionModel*));
	for(g=0; g<nG; g++)
		modelT[g] = getEvolutionModelProtStored(tree, gamma[g]);
	column = (int*) malloc(align->number*sizeof(int));
	for(pos=0; pos<align->size; pos++) {
		for(n=0; n<align->number; n++)
			column[n] = align->sequence[n][pos];
		logLike += log(getGammaLikelihoodColumn(tree, modelT, gamma, column));
	}
	free((void*)column);
	for(g=0; g<nG; g++)
		modelT[g]->freeModel(modelT[g]);
	free((void*)modelT);
	return logLike;
}

double getLikelihoodAlignment(TypeTree *tree, TypeEvolutionModel *model, TypeAlignment *align) {
	int pos, *column;
	double logLike = 0.;
	column = (int*) malloc(align->number*sizeof(int));
	for(pos=0; pos<align->size; pos++) {
		int n;
		for(n=0; n<align->number; n++)
			column[n] = align->sequence[n][pos];
		logLike += log(getLikelihoodColumn(tree, model, column));
	}
	free((void*)column);
	return exp(logLike);
}

double toOpt(unsigned n, const double *x, double *grad, void *data) {
	double *gamma, result;
	gamma = getGamma(x[0], ((TypeOptData *)data)->nCat);
	result = getGammaLogLikelihoodAlignment(((TypeOptData *)data)->tree, ((TypeOptData *)data)->model, gamma, ((TypeOptData *)data)->align);
	free((void*)gamma);
    return result;
}

double estimateGammaParameter(TypeTree *tree, TypeEvolutionModel *model, void *randG, TypeAlignment *align, int nCat) {
	double max, x, rate;
	TypeNLOptOption *option=getNLOption();
	TypeOptData data;
	int t, status;
    data.model = model;
    data.tree = tree;
    data.align = align;
    data.nCat = nCat;
	nlopt_opt opt;
	opt = nlopt_create(NLOPT_ALGO, 1); 
	nlopt_set_lower_bounds1(opt, 0.105);
	nlopt_set_upper_bounds1(opt, 100);
	nlopt_set_max_objective(opt, toOpt, &data);
	nlopt_set_xtol_rel(opt, option->tolOptim);
	nlopt_set_maxeval(opt, option->maxIter);
	rate = 1.;
	x = 1.;
	max = -HUGE_VAL;
	for(t=0; t<option->trials; t++) {
		double tmp;
		fillRandStartNLOpt(opt, randG, &x);
		if(((status = nlopt_optimize(opt, &x, &tmp))>=0)) {
//printf("Par %.2le\t%.2le (%.2le)\n", x, tmp, max);
			if(tmp > max) {
				rate = x;
				max = tmp;
			}
		 }
	}
//printf("Fin %.2le\t%.2le\n", rate, max);
	if(max == -HUGE_VAL) {
		fprintf(stderr, "Numerical optimisation failed!\n");
		exit(1);
	}
	nlopt_destroy(opt);
	return rate;
}
