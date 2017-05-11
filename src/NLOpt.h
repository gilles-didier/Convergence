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




#ifndef NLOptF
#define NLOptF

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <nlopt.h>

typedef struct NLOPT_OPTION {
    double tolOptim;
    int trials, maxIter;
} TypeNLOptOption;

void initializeRandNLOpt();
void freeRandNLOpt();
void fillRandStartNLOpt(nlopt_opt opt, void *randG, double *x);
void setNLOption(TypeNLOptOption opt);
TypeNLOptOption *getNLOption();
void fprintNLoptOption(FILE *f, TypeNLOptOption *option);
void sprintNLoptOption(char *buffer, TypeNLOptOption *option);
void fprintNLoptOptionTag(FILE *f, TypeNLOptOption *option);
void fscanNLoptOptionTag(FILE *f, TypeNLOptOption *option);

//#define NLOPT_ALGO NLOPT_GN_ISRES
//#define NLOPT_ALGO NLOPT_GN_ESCH
#define NLOPT_ALGO NLOPT_LN_BOBYQA
//#define NLOPT_ALGO NLOPT_LN_COBYLA
//#define NLOPT_ALGO NLOPT_AUGLAG


#ifdef __cplusplus
}
#endif

#endif
