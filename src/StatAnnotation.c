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
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>

#include "Utils.h"
#include "StatAnnotation.h"

static int compareSignificant(const void* a, const void* b);
static void fillRandomList(const gsl_rng *rg, int *list, int n, int max);

TypeStatAnnotation *newStatAnnotation(TypeAnnotation *annot, char **name) {
    int *indexA, i, j, size, n, ind, sum, number;
    TypeLexiTree *dict;
    TypeStatAnnotation *sap;
	
	if(name == NULL) {
		sap = (TypeStatAnnotation*) malloc(sizeof(TypeStatAnnotation));
		size = 0;
		sap->sizeG = annot->sizeG;
		sap->sizeA = annot->sizeA;
		sap->mat = (int**) malloc(sap->sizeG*sizeof(int*));
		sap->nameG = (char**) malloc(sap->sizeG*sizeof(char*));
		sap->nameA = (char**) malloc(sap->sizeA*sizeof(char*));
		for(i=0; i<sap->sizeG; i++)
			sap->mat[i] = (int*) malloc(sap->sizeA*sizeof(int));
		for(i=0; i<annot->sizeG; i++) {
			sap->nameG[i] = strdpl(annot->nameG[i]);
			sap->mat[i] = (int*) malloc(sap->sizeA*sizeof(int));
			for(j=0; j<sap->sizeA; j++)
				sap->mat[i][j] = annot->mat[i][j];
		}
		for(i=0; i<annot->sizeA; i++)
			sap->nameA[i] = strdpl(annot->nameA[i]);
		sum = 0;
		for(i=0; i<sap->sizeG; i++)
			for(j=0; j<sap->sizeA; j++)
				if(sap->mat[i][j])
					sum++;
		sap->startG = (int*) malloc((sap->sizeG+1)*sizeof(int));
		sap->itemG = (int*) malloc(sum*sizeof(int));
		ind = 0;
		for(i=0; i<sap->sizeG; i++) {
			sap->startG[i] = ind;
			for(j=0; j<sap->sizeA; j++)
				if(sap->mat[i][j])
					sap->itemG[ind++] = j;
		}
		sap->startG[sap->sizeG] = ind;
		sap->startA = (int*) malloc((sap->sizeA+1)*sizeof(int));
		sap->itemA = (int*) malloc(sum*sizeof(int));
		ind = 0;
		for(j=0; j<sap->sizeA; j++) {
			sap->startA[j] = ind;
			for(i=0; i<sap->sizeG; i++)
				if(sap->mat[i][j])
					sap->itemA[ind++] = i;
		}
		sap->startA[sap->sizeA] = ind;
		return sap;
    } else {
		dict = getDictFromTable(annot->nameG, annot->sizeG);
		for(number=0; name[number]!=NULL; number++)
			;
		sap = (TypeStatAnnotation*) malloc(sizeof(TypeStatAnnotation));
		indexA = (int*) malloc(annot->sizeG*sizeof(int));
		size = 0;
		for(i=0; i<annot->sizeG; i++)
			indexA[i] = END_INT;
		for(i=0; i<number; i++)
			if((n = findWordLexi(name[i], dict)) != END_INT)
				indexA[n] = size++;
		sap->sizeG = size;
		sap->sizeA = annot->sizeA;
		sap->mat = (int**) malloc(sap->sizeG*sizeof(int*));
		sap->nameG = (char**) malloc(sap->sizeG*sizeof(char*));
		sap->nameA = (char**) malloc(sap->sizeA*sizeof(char*));
		for(i=0; i<sap->sizeG; i++)
			sap->mat[i] = (int*) malloc(sap->sizeA*sizeof(int));
		for(i=0; i<annot->sizeG; i++)
			if(indexA[i] != END_INT) {
				sap->nameG[indexA[i]] = strdpl(annot->nameG[i]);
				for(j=0; j<sap->sizeA; j++)
					sap->mat[indexA[i]][j] = annot->mat[i][j];
			}
		for(i=0; i<annot->sizeA; i++)
			sap->nameA[i] = strdpl(annot->nameA[i]);
		sum = 0;
		for(i=0; i<sap->sizeG; i++)
			for(j=0; j<sap->sizeA; j++)
				if(sap->mat[i][j])
					sum++;
		sap->startG = (int*) malloc((sap->sizeG+1)*sizeof(int));
		sap->itemG = (int*) malloc(sum*sizeof(int));
		ind = 0;
		for(i=0; i<sap->sizeG; i++) {
			sap->startG[i] = ind;
			for(j=0; j<sap->sizeA; j++)
				if(sap->mat[i][j])
					sap->itemG[ind++] = j;
		}
		sap->startG[sap->sizeG] = ind;
		sap->startA = (int*) malloc((sap->sizeA+1)*sizeof(int));
		sap->itemA = (int*) malloc(sum*sizeof(int));
		ind = 0;
		for(j=0; j<sap->sizeA; j++) {
			sap->startA[j] = ind;
			for(i=0; i<sap->sizeG; i++)
				if(sap->mat[i][j])
					sap->itemA[ind++] = i;
		}
		sap->startA[sap->sizeA] = ind;
		freeLexiTree(dict);
		free((void*)indexA);
		return sap;
	}
}

void freeStatAnnotation(TypeStatAnnotation *sap) {
	if(sap != NULL) {
		int i;
		if(sap->nameA != NULL) {
			for(i=0; i<sap->sizeA; i++)
				if(sap->nameA[i] != NULL)
					free((void*)sap->nameA[i]);
			free((void*)sap->nameA);
		}
		if(sap->nameG != NULL) {
			for(i=0; i<sap->sizeG; i++)
				if(sap->nameG[i] != NULL)
					free((void*)sap->nameG[i]);
			free((void*)sap->nameG);
		}
		if(sap->mat != NULL) {
			for(i=0; i<sap->sizeG; i++)
				if(sap->mat[i] != NULL)
					free((void*)sap->mat[i]);
			free((void*)sap->mat);
		}
		if(sap->startA != NULL)
			free((void*)sap->startA);
		if(sap->startG != NULL)
			free((void*)sap->startG);
		if(sap->itemA != NULL)
			free((void*)sap->itemA);
		if(sap->itemG != NULL)
			free((void*)sap->itemG);
		free((void*)sap);
	}
}

void fprintStatAnnotationInfo(FILE *f, TypeStatAnnotation *sap) {
    fprintf(f, "%d genes %d annotations\n", sap->sizeG, sap->sizeA);
}

/*Fisher C(n_1, k) C(n-n_1, t - k) / C(n_1 + n_2, t) n_1 number of elt with ontologie o in the atom, n card of the atom, k card of o, t total*/

double getFisher(int k, int n1, int n2, int t) {
    double res = 0.;
    res = gsl_cdf_hypergeometric_Q(k-1, n1, n2, t);
    if(res == 0.)
        printf("Probable overflow : k %d n1 %d n2 %d t %d %lE\n", k, n1, n2, t,res);
    return res;
}

TypeSignificantTable *getSignificantTableList(TypeStatAnnotation *sap, int *list) {
    int *eff, **gene, j, k, c, inc_buffer = sap->startG[sap->sizeG]/10, size_buffer = inc_buffer, number = 0;
    TypeSignificantTable *res;
    double bonf = ((double) sap->sizeA);
    res = (TypeSignificantTable*) malloc(sizeof(TypeSignificantTable));
    res->table = (TypeSignificant*) malloc(size_buffer*sizeof(TypeSignificant));
    res->size = 0;
    res->total = sap->sizeG;
    eff = (int*) malloc(sap->sizeA*sizeof(int));
    gene = (int**) malloc(sap->sizeA*sizeof(int*));
	for(c=0; c<sap->sizeA; c++) {
		eff[c] = 0;
		gene[c] = NULL;
	}
	number = 0;
	for(j=0; list[j] != END_INT; j++)
		number++;
	if(number == 0)
		return NULL;
	for(j=0; list[j] != END_INT; j++) {
		for(k=sap->startG[list[j]]; k<sap->startG[list[j]+1]; k++) {
			if(gene[sap->itemG[k]] == NULL)
				gene[sap->itemG[k]] = (int*) malloc(number*sizeof(int));
			gene[sap->itemG[k]][eff[sap->itemG[k]]++] = list[j];
		}
	}
	for(c=0; c<sap->sizeA; c++) {
		if(eff[c]) {
			double tmp = getFisher(eff[c], sap->startA[c+1]-sap->startA[c], sap->sizeG-(sap->startA[c+1]-sap->startA[c]), number), corrected = tmp*bonf;
			if(1) {
				if(res->size>=size_buffer) {
					size_buffer += inc_buffer;
					res->table = (TypeSignificant*) realloc(res->table, size_buffer*sizeof(TypeSignificant));
				}
				res->table[res->size].effC = number;
				res->table[res->size].effO = sap->startA[c+1]-sap->startA[c];
				res->table[res->size].effCO = eff[c];
				res->table[res->size].ontology = c;
				res->table[res->size].fisher = tmp;
				res->table[res->size].corrected = corrected;
				res->table[res->size].gene = (int*) realloc((void*) gene[c], eff[c]*sizeof(int));
				res->size++;
			}
		}
	}
	qsort(&(res->table[0]), res->size, sizeof(TypeSignificant), compareSignificant);
    free((void*)eff);
    free((void*)gene);
    return res;
}

void freeSignificantTable(TypeSignificantTable *sig) {
	if(sig != NULL) {
		if(sig->table != NULL) {
			int i;
			for(i=0; i<sig->size; i++)
				if(sig->table[i].gene != NULL)
					free((void*)sig->table[i].gene);	
			free((void*)sig->table);
		}
		free((void*)sig);
	}
}

void fprintSignificantTable(FILE *f, TypeSignificantTable *sig, int totalGene, char **nameA, char **nameG, TypeOntologyInfo *info) {
    int i;
    TypeLexiTree *dict;

    if(sig->size == 0)
        return;
    if(info != NULL)
        dict = 	getDictFromTable(info->id, info->size);
    else
        dict = NULL;
    for(i=0; i<sig->size; i++) {
        int n;
        fprintf(f, "%s\t%.3lE\t(%d %d %d %d)", nameA[sig->table[i].ontology], sig->table[i].fisher, sig->table[i].effCO, sig->table[i].effC, sig->table[i].effO, totalGene);
        if(info != NULL) {
			if((n = findWordLexi(nameA[sig->table[i].ontology], dict)) != END_INT)
				fprintf(f, "\t%s", info->name[n]);
			else
				fprintf(f, "\t?");
		}
		fprintf(f, "\t%s", nameG[sig->table[i].gene[0]]);
		for(n=1; n<sig->table[i].effCO; n++)
			fprintf(f, ", %s", nameG[sig->table[i].gene[n]]);
       fprintf(f, "\n");
    }
    if(dict != NULL)
        freeLexiTree(dict);
}

void fprintSignificantLine(FILE *f, TypeSignificant line, int totalGene, char **nameA, char **nameG, TypeOntologyInfo *info) {
	int n;
	TypeLexiTree *dict;
	if(info != NULL)
		dict = 	getDictFromTable(info->id, info->size);
	else
		dict = NULL;
	fprintf(f, "%s\t%.3lE\t(%d %d %d %d)", nameA[line.ontology], line.fisher, line.effCO, line.effC, line.effO, totalGene);
	if(info != NULL) {
		if((n = findWordLexi(nameA[line.ontology], dict)) != END_INT)
			fprintf(f, "\t%s", info->name[n]);
		else
			fprintf(f, "\t?");
	}
	fprintf(f, "\t%s", nameG[line.gene[0]]);
	for(n=1; n<line.effCO; n++)
		fprintf(f, ", %s", nameG[line.gene[n]]);
	fprintf(f, "\n");
	if(dict != NULL)
		freeLexiTree(dict);
}


int compareSignificant(const void* a, const void* b) {
    return compareDouble(&(((TypeSignificant*)a)->corrected), &(((TypeSignificant*)b)->corrected));
}

void fillRandomList(const gsl_rng *rg, int *list, int n, int max) {
	int *table, i;
	table = (int*) malloc(max*sizeof(int));
	for(i=0; i<max; i++)
		table[i] = i;
	for(i=0; i<n; i++) {
		int tmp, ind = gsl_rng_uniform_int(rg, max-i);
		list[i] = table[ind];
		tmp = table[ind];
		table[ind] = table[max-i-1];
		table[max-i-1] = tmp;
	}
	list[n] = END_INT;
	free((void*)table);
}

void fillCorrectedEmpirical(TypeStatAnnotation *sap, TypeSignificantTable *res, int nList, int nSim) {
	int i, j, *list, *occ;
	gsl_rng *rg = gsl_rng_alloc(gsl_rng_taus);
	list = (int*) malloc((nList+1)*sizeof(int));
	occ = (int*) malloc(res->size*sizeof(int));
	for(j=0; j<res->size; j++)
		occ[j] = 0;
	for(i=0; i<nSim; i++) {
		TypeSignificantTable *sig;
		fillRandomList(rg, list, nList, res->total);
		sig = getSignificantTableList(sap, list);
		for(j=0; j<res->size && j<sig->size; j++)
			if(res->table[j].fisher >= sig->table[j].fisher)
				occ[j]++;
		for(; j<res->size; j++)
				occ[j]++;
		freeSignificantTable(sig);
	}
	for(j=0; j<res->size; j++)
		res->table[j].corrected = ((double)occ[j])/((double)nSim);
	free((void*)list);
	free((void*)occ);
	gsl_rng_free(rg);
}
