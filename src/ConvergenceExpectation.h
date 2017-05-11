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




#ifndef ConvergenceExpectationF
#define ConvergenceExpectationF

#include "EvolutionModel.h"
#include "Tree.h"

#define CALIB_VALUE 1.1
#define CALIB_SYMBOL_PROT 10

double calibRate(TypeTree *tree, int *character, TypeEvolutionModel *model, void *randG, int a, double ce);
int *getSymbols(TypeTree *tree, int *character, int cardinal, int *column);
double getConvergenceIndex(TypeTree *tree, int *character, TypeEvolutionModel *model, int *column);
double getConvergenceIndexList(TypeTree *tree, int *character, TypeEvolutionModel *model, int *column, int *al);
double getConvergenceExpectation(TypeTree *tree, int *character, TypeEvolutionModel *model, int *column, int AA);

#endif
