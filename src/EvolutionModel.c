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
#include <gsl/gsl_randist.h>

#include "EvolutionModel.h"

static gsl_rng *getRandGen(void *RandG);

static gsl_rng *randGen = NULL;

void initializeDrawEvolution() {
	randGen = gsl_rng_alloc(gsl_rng_taus);
}

void freeDrawEvolution() {
	if(randGen != NULL)
		gsl_rng_free(randGen);
}

void freeModelBasic(TypeEvolutionModel *model) {
	if(model->data!=NULL)
		free(model->data);
	free((void*) model);
}

int setTimeBasic(TypeEvolutionModel *model, double t) {
	if(t>=0.) {
		model->t = t;
		return 1;
	} else
		return 0;
}

double getMinSteady(TypeEvolutionModel *model) {
	double min;
	int i;
	min = model->getLogSteady(model, 0)/model->t;
	for(i=1; i<model->cardinal; i++)
		if(model->getLogSteady(model, i)/model->t<min)
			min = model->getLogSteady(model, i)/model->t;
	return min;
}

double getMaxSteady(TypeEvolutionModel *model) {
	double max;
	int i;
	max = model->getLogSteady(model, 0)/model->t;
	for(i=1; i<model->cardinal; i++)
		if(model->getLogSteady(model, i)/model->t>max)
			max = model->getLogSteady(model, i)/model->t;
	return max;
}


int drawInitial(TypeEvolutionModel *model, void *randG) {
	double r, sum=0.;
	int d;
	gsl_rng *rg = getRandGen(randG);
	r = gsl_rng_uniform(rg);
	for(d=0; d<model->cardinal && (sum += model->getInitial(model, d))<r; d++)
		;
	return d;
}

int drawTransitionFrom(TypeEvolutionModel *model, void *randG, int l) {
	double r, sum=0.;
	int d;
	gsl_rng *rg = getRandGen(randG);
	r = gsl_rng_uniform(rg);
	for(d=0; d<model->cardinal && (sum += model->getTransition(model, l, d))<r; d++)
		;
	return d;
}


gsl_rng *getRandGen(void *randG) {
	gsl_rng *rg;
	if(randG != NULL) {
		rg = (gsl_rng *) randG;
	} else {
		if(randGen == NULL)
			initializeDrawEvolution();
		rg = randGen;
	}
	return rg;
}
