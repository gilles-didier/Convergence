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
#include <ctype.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "NLOpt.h"

#define TRIALS 5
#define MAX_ITER 50
#define TOLERANCE_OPTIM 0.001

#define TAG_TRI "TRI"
#define TAG_TOL "TOL"
#define TAG_ITE "ITE"
#define SIZE_TAG 20
#define SIZE_VAL 100

static gsl_rng *randGen = NULL;
static TypeNLOptOption option={.tolOptim=TOLERANCE_OPTIM,.trials=TRIALS,.maxIter=MAX_ITER};

static gsl_rng *getRandGen(void *randG);


void fprintNLoptOptionTag(FILE *f, TypeNLOptOption *option) {
    fprintf(f, ":%s %d\n", TAG_TRI, option->trials);
    fprintf(f, ":%s %lE\n", TAG_TOL, option->tolOptim);
    fprintf(f, ":%s %d\n", TAG_ITE, option->maxIter);
}

void initializeRandNLOpt() {
	randGen = gsl_rng_alloc(gsl_rng_taus);
}

void freeRandNLOpt() {
	if(randGen != NULL)
		gsl_rng_free(randGen);
}

void fscanNLoptOptionTag(FILE *f, TypeNLOptOption *option) {
    char c, tag[SIZE_TAG+1], val[SIZE_VAL+1];
    for(c=fgetc(f); c!=EOF && isspace(c); c=fgetc(f));
    while(c == ':') {
        int i;
        c=fgetc(f);
        for(i=0; c!=EOF && !isspace(c) && i<SIZE_TAG; c=fgetc(f))
            tag[i++] = c;
        tag[i] = '\0';
        if(i>=SIZE_TAG) {
            fprintf(stderr, "Error when reading an optimizer options file - Tag too long:\n%s...\n", tag);
            exit(1);
        }
        for(; c!=EOF && isspace(c); c=fgetc(f));
        for(i=0; c!=EOF && !isspace(c) && i<SIZE_VAL; c=fgetc(f))
            val[i++] = c;
        val[i] = '\0';
        if(i>=SIZE_VAL) {
            fprintf(stderr, "Error when reading an optimizer options file - value too long:\n%s...\n", val);
            exit(1);
        }
        if(strcmp(tag, TAG_TRI) == 0)
            option->trials = atoi(val);
        if(strcmp(tag, TAG_TOL) == 0)
            option->tolOptim = atof(val);
        if(strcmp(tag, TAG_ITE) == 0)
            option->maxIter = atoi(val);
        for(; c!=EOF && isspace(c); c=fgetc(f));
    }
}

void fprintNLoptOption(FILE *f, TypeNLOptOption *option) {
    fprintf(f, "Optimizer runs %d trials and stops with tolerance %.lE or after more than %d iterations.\n", option->trials, option->tolOptim, option->maxIter);
}

void sprintNLoptOption(char *buffer, TypeNLOptOption *option) {
    buffer += sprintf(buffer, "Optimizer runs %d trials and stops with tolerance %.lE or after more than %d iterations.\n", option->trials, option->tolOptim, option->maxIter);
}

void setNLOption(TypeNLOptOption opt) {
	option = opt;
}

TypeNLOptOption *getNLOption() {
	return &option;
}

void fillRandStartNLOpt(nlopt_opt opt, void *randG, double *x) {
	int n, i;
	double *lb, *ub;
	gsl_rng *rg = getRandGen(randG);
	n = nlopt_get_dimension(opt);
	lb = (double*) malloc(n*sizeof(double));
	ub = (double*) malloc(n*sizeof(double));
	nlopt_get_lower_bounds(opt, lb);
	nlopt_get_upper_bounds(opt, ub);
	for(i=0; i<n; i++) {
		if(isinf(ub[i]))
			ub[i] = 1000000.;
		if(isinf(lb[i]))
			lb[i] = -1000000.;
		x[i] = lb[i]+gsl_rng_uniform(rg)*(ub[i]-lb[i]);
	}
	free((void*)lb);
	free((void*)ub);
}

gsl_rng *getRandGen(void *randG) {
	gsl_rng *rg;
	if(randG != NULL) {
		rg = (gsl_rng *) randG;
	} else {
		if(randGen == NULL)
			initializeRandNLOpt();
		rg = randGen;
	}
	return rg;
}
