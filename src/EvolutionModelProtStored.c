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
#include "EvolutionModelProtStored.h"

static int setTimeProtStored(TypeEvolutionModel *model, double t, void *user);
static void freeModelProtStored(TypeEvolutionModel *model);

int setTimeProtStored(TypeEvolutionModel *model, double t, void *user) {
	((TypeProtStoredData*)((TypeProtData*)model->data)->user)->cur = *((int*)user);
	model->t = ((TypeProtStoredData*)((TypeProtData*)model->data)->user)->time[((TypeProtStoredData*)((TypeProtData*)model->data)->user)->cur];
	((TypeProtData*)model->data)->trans = ((TypeProtStoredData*)((TypeProtData*)model->data)->user)->table+((TypeProtStoredData*)((TypeProtData*)model->data)->user)->cur*CARDINAL_PROT2;
	return 1.;
}

void freeModelProtStored(TypeEvolutionModel *model) {
	free((void*) ((TypeProtStoredData*)((TypeProtData*)model->data)->user)->time);
	free((void*) ((TypeProtStoredData*)((TypeProtData*)model->data)->user)->table);
	free((void*)((TypeProtData*)model->data)->user);
	free((void*)((TypeProtData*)model->data));
	free((void*) model);
}

TypeEvolutionModel *getEvolutionModelProtStored(TypeTree *tree, double rate) {
	TypeEvolutionModel *model;
	int n;
	TypeProtStoredData *user;
	model = getEvolutionModelProt();
	free((void*)((TypeProtData*)model->data)->trans);
	user = (TypeProtStoredData*) malloc(sizeof(TypeProtStoredData));
	user->size = tree->size;
	user->table = (double*) malloc(tree->size*CARDINAL_PROT2*sizeof(double));
	user->time = (double*) malloc(tree->size*sizeof(double));
	for(n=0; n<tree->root; n++) {
		user->time[n] = rate*tree->time[n];
		expMatrix(((TypeProtData*)model->data), user->time[n], user->table+n*CARDINAL_PROT2);
	}
	for(n=tree->root+1; n<tree->size; n++) {
		user->time[n] = rate*tree->time[n];
		expMatrix(((TypeProtData*)model->data), user->time[n], user->table+n*CARDINAL_PROT2);
	}
	((TypeProtData*)model->data)->user = (void*) user;
	setTimeProtStored(model, user->time[tree->node[tree->root].child], (void*) &(tree->node[tree->root].child));
	model->freeModel = freeModelProtStored;
	model->setTime = setTimeProtStored;
	return model;
}
