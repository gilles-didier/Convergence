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




#ifndef StatAnnotationF
#define StatAnnotationF

#include <stdlib.h>
#include <stdio.h>

#include "Annotation.h"

#define INC_SIZE 100
#define MAX_CURRENT 1000000

#ifdef __cplusplus
extern "C" {
#endif

typedef struct TYPE_STAT_ANNOTATION {
    char **nameG, **nameA;
    int sizeG, sizeA, **mat, *startG, *itemG, *startA, *itemA;
} TypeStatAnnotation;

typedef struct TYPE_SIGNIFICANT {
    int ontology, effC, effO, effCO, *gene;
    double fisher, corrected;
} TypeSignificant;

typedef struct TYPE_SIGNIFICANT_TABLE {
    int size, total;
    TypeSignificant *table;
} TypeSignificantTable;

typedef struct TYPE_STAT_FISHER {
    int size;
    double *fisher, *corrected;
} TypeStatFisher;

TypeStatAnnotation *newStatAnnotation(TypeAnnotation *annot, char **name);
double getFisher(int n1, int n, int k, int t);
TypeSignificantTable *getSignificantTableList(TypeStatAnnotation *sap, int *list);
void fprintSignificantTable(FILE *f, TypeSignificantTable *sig, int totalGene, char **nameA, char **nameG, TypeOntologyInfo *info);
void fprintSignificantLine(FILE *f, TypeSignificant line, int totalGene, char **nameA, char **nameG, TypeOntologyInfo *info);
void fillCorrectedEmpirical(TypeStatAnnotation *sap, TypeSignificantTable *res, int nList, int nSim);
void freeSignificantTable(TypeSignificantTable *sig);
void freeStatAnnotation(TypeStatAnnotation *sap);

#ifdef __cplusplus
}
#endif

#endif
