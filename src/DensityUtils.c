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
#include <math.h>
#include "DensityUtils.h"

void printDensity(FILE *f, double *data, int size, double ec) {
	double min, max;
	int i, *occ, tot = 0;

	min = data[0];
	max = data[0];
	for(i=1; i<size; i++) {
		if((data[i]+0.5*ec)<min)
			min = data[i]+0.5*ec;
		if((data[i]+0.5*ec)>max)
			max = data[i]+0.5*ec;
	}
	min = floor(min/ec)*ec;
	max = ceil(max/ec)*ec;
	tot = (int) ceil((max-min)/ec);
	occ = (int*) malloc(tot*sizeof(int));
	for(i=0; i<tot; i++)
		occ[i] = 0;
	for(i=0; i<size; i++)
		occ[(int)floor((data[i]-min+0.5*ec)/ec)]++;
	for(i=0; i<tot; i++)
		fprintf(f, "%le\t%d\t%lE\n", ((double)i)*ec+min, occ[i], ((double)occ[i])/((double)size));
	free((void*)occ);
}

void printDensityBis(FILE *f, double *data, int size) {
	double min, max, ec;
	int i, *occ, tot = 0;

	min = data[0];
	max = data[0];
	for(i=1; i<size; i++) {
		if(data[i]<min)
			min = data[i];
		if(data[i]>max)
			max = data[i];
	}
	ec = pow(10.,ceil(log(max-min)/log(10.))-2.);
	min += 0.5*ec;
	max += 0.5*ec;
	min = floor(min/ec)*ec;
	max = ceil(max/ec)*ec;
	tot = (int) ceil((max-min)/ec);
	occ = (int*) malloc(tot*sizeof(int));
	for(i=0; i<tot; i++)
		occ[i] = 0;
	for(i=0; i<size; i++)
		occ[(int)floor((data[i]+0.5*ec-min)/ec)]++;
	for(i=0; i<tot; i++)
		fprintf(f, "%le\t%d\t%lE\n", ((double)i)*ec+min, occ[i], ((double)occ[i])/((double)size));
	free((void*)occ);
}

TypeEmpiricalDensity *getDensityBis(double *data, int size) {
	double min, max, ec;
	int i, *occ, tot = 0;
	TypeEmpiricalDensity *res;
	res = (TypeEmpiricalDensity*) malloc(sizeof(TypeEmpiricalDensity));
	min = data[0];
	max = data[0];
	for(i=1; i<size; i++) {
		if(data[i]<min)
			min = data[i];
		if(data[i]>max)
			max = data[i];
	}
	ec = pow(10.,ceil(log(max-min)/log(10.))-2.);
	min = floor(min/ec)*ec;
	max = ceil((max+0.5*ec)/ec)*ec;
	tot = (int) ceil((max-min)/ec);
	occ = (int*) malloc(tot*sizeof(int));
	for(i=0; i<tot; i++)
		occ[i] = 0;
	for(i=0; i<size; i++) {
		occ[(int)floor((data[i]+0.5*ec-min)/ec)]++;
	}
	res->size = tot;
	res->x = (double*) malloc(res->size*sizeof(double));
	res->y = (double*) malloc(res->size*sizeof(double));
	for(i=0; i<tot; i++) {
		res->x[i] = ((double)i)*ec+min;
		res->y[i] = ((double)occ[i])/((double)size);
	}
	free((void*)occ);
	return res;
}
