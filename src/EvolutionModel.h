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




#ifndef EvolutionModelF
#define EvolutionModelF

typedef struct  EVOLUTION_MODEL {
	int cardinal;
	double t;
	void *data;
	void (*freeModel)(struct  EVOLUTION_MODEL*);
	int (*setTime)(struct  EVOLUTION_MODEL*, double, void *user);
	double (*getTransition)(struct  EVOLUTION_MODEL*, int, int);
	double (*getSteady)(struct  EVOLUTION_MODEL*, int);
	double (*getInitial)(struct  EVOLUTION_MODEL*, int);
	double (*getLogTransition)(struct  EVOLUTION_MODEL*, int, int);
	double (*getLogSteady)(struct  EVOLUTION_MODEL*, int);
	double (*getLogInitial)(struct  EVOLUTION_MODEL*, int);
} TypeEvolutionModel;

void freeModelBasic(TypeEvolutionModel *model);
int setTimeBasic(TypeEvolutionModel *model, double t);
double getMinSteady(TypeEvolutionModel *model);
double getMaxSteady(TypeEvolutionModel *model);
void initializeDrawEvolution();
void freeDrawEvolution();
int drawInitial(TypeEvolutionModel *model, void *ran);
int drawTransitionFrom(TypeEvolutionModel *model, void *ran, int l);

#endif
