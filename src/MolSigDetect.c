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
#include <pthread.h>
#include <signal.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "Tree.h"
#include "Alignment.h"
#include "Character.h"
#include "EvolutionModel.h"
#include "EvolutionModelProt.h"
#include "EvolutionModelProtStored.h"
#include "ConvergenceExpectation.h"
#include "ColumnLikelihood.h"
#include "AlignmentLikelihood.h"
#include "NLOpt.h"
#include "DensityUtils.h"
#include "AATree.h"

#define HELPMESSAGE "--------------------------\n\nNAME\n\tmsd - detection of molecular signatures of convergence events\n\t\nSYNOPSIS\n\tmsd [OPTIONS] [FILE TREE] [FILE CONV.] [ALIGN. FILE 1] [ALIGN. FILE 2] ...\n\nDESCRIPTION\n\treturn a table where each line displays the result of the detection of molecular signature of each alignment fils.\n\tOptions are\n\t\n\t-n [FILE]\n\t\tset the optimisation options to those contained in [FILE].\n\t-m [FILE]\n\t\tset the evolution model to that contained in [FILE].\n\t-o [NAME]\n\t\tset the name of the output file to [NAME]. By default this is 'output.txt'.\n\t-e [VALUE]\n\t\tsave an extra file named 'significant.txt', which contains all the names of alignments with a corrected p-value unde [VALUE].\n\t-t [NUMBER]\n\t\tset the maximal number of simultaneous threads to [NUMBER.\n\t-h\n\t\tdisplay help\n\n--------------------------"
#define EXTRACT_NAME "significant.txt"
#define OUTPUT_NAME "output.txt"
#define ERROR_NAME "warnings.txt"
#define REPORT_NAME "report.txt"
#define ITER_TP 1000000
#define MAX_SIMULATION 30000000
//#define MAX_SIMULATION 500000

//./msd -t 40 -o thomas_comp.txt -e 0.05 ../data/Thomas/species_tree.nwk ../data/Thomas/car.csv /home/olivier/these/data/matthew/alignmentProt/*.fasta
//./msd -t 40 -o thomas_comp_X.txt -e 0.05 -s 0.001 ../data/Thomas/species_tree.nwk ../data/Thomas/car.csv ~/OwnCloud/data/Thomas/alignments/*.fasta

typedef struct CONVERGENCE_DATA {
	int nSite, *site, unknown, status;
	double pvalue;
} TypeConvergenceData;

typedef struct ALIGNMENT_DATA {
	char *name;
	int size, nSite, *site, status;
	double pvalue, pvalueBH, unknown;
} TypeAlignmentData;

typedef struct INPUT_DATA {
	TypeTree *tree;
	int *character;
} TypeInput;

typedef struct THREAD_PARAMETER {
	int *number;
	pthread_mutex_t *mutex_result;
	pthread_mutex_t *mutex_number;
	pthread_mutex_t *mutex_error;
	pthread_mutex_t *mutex_report;
	pthread_cond_t *cond_number;
	FILE *fe, *fr;
	char *filename, *name;
	double threshold;
	TypeLexiTree *dict;
	TypeTree tree;
	TypeEvolutionModel *model;
	int *character, nGam;
	TypeAlignmentData *result;
} TypeThreadParameter;

static TypeConvergenceData getAlignmentConvergence(TypeTree *tree, int *character, TypeEvolutionModel *modelC, TypeAlignment *align, TypeThreadParameter *data, gsl_rng *rg);
static void reindexAlignments(TypeAlignment *alignment, TypeLexiTree *dict);
static int compareAlignmentData(void const *a, void const *b);
static void correctBenjaminiHochberg(TypeAlignmentData *alignData, int size);
static void fprintAlignmentData(FILE *f, TypeAlignmentData *alignData, int size);
static void *threadAlignment(void *data);
static char *getAlignmentName(char *filename);
static int fillUnknown(TypeAlignment *align, double *unknown);
static TypeInput fixInput(TypeTree *tree,  int *character, TypeAlignment *align, TypeLexiTree *dictTree);
static TypeEvolutionModel *cloneEvolutionModelProtStored(TypeEvolutionModel *model);
static void freeCloneEvolutionModelProtStored(TypeEvolutionModel *model);
static char getStatusChar(int status);

unsigned long int random_seed = 0;
int nb_sim = 10000;
double confidence = 0.9999, sig_thre = 0.005;

int main(int argc, char *argv[]) {
	char modelFileName[STRING_SIZE], outputFileName[STRING_SIZE], option[256];
	FILE *f, *fe, *fr;
	int i, j, nA, *character, nGam = 4, nT, maxT = 40, cont, size, writeIterF = 0;
	TypeTree *tree;
	TypeEvolutionModel *model;
	TypeLexiTree *dict;
	double threshold = 0.001, best, extractThreshold = END_DOUBLE;
	TypeAlignmentData *alignData;
	pthread_mutex_t mutexR = PTHREAD_MUTEX_INITIALIZER, mutexN = PTHREAD_MUTEX_INITIALIZER, mutexE = PTHREAD_MUTEX_INITIALIZER, mutexP = PTHREAD_MUTEX_INITIALIZER;
	pthread_cond_t condN = PTHREAD_COND_INITIALIZER;
	time_t t0, t1;

	time(&t0);
	strcpy(outputFileName, OUTPUT_NAME);
	modelFileName[0] = '\0';
	for(i=0; i<256; i++)
		option[i] = 0;	 
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['n']) {
			option['n'] = 0;
			FILE *fopt;
			if((i+1)<argc) {
				if((fopt = fopen(argv[++i], "r"))) {
					TypeNLOptOption nloptOption;
					fscanNLoptOptionTag(fopt, &nloptOption);
					fclose(fopt);
					setNLOption(nloptOption);
				} else {
					fprintf(stderr, "Can't open file %s\n", argv[++i]);
					exit(1);
				}
			} else {
				fprintf(stderr, "File name missing after option -o\n");
				exit(1);
			}
		}
		if(option['y']) {
			option['y'] = 0;
			writeIterF = 1;
		}
		if(option['m']) {
			option['m'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%s", modelFileName) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -m");
		}
		if(option['o']) {
			option['o'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%s", outputFileName) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a file name is required after option -o");
		}
		if(option['e']) {
			option['e'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &extractThreshold) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a value is required after option -e");
		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &threshold) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a value is required after option -s");
		}
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &maxT) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is required after option -t");
		}
		if(option['a']) {
			option['a'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lu", &random_seed) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is required after option -a");
		}
		if(option['i']) {
			option['i'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nb_sim) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is required after option -i");
		}
		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}
	if(i<argc && (f = fopen(argv[i], "r"))) {
		tree = readTree(f);
		reindexLeavesFirst(tree);
		fclose(f);
	} else {
		fprintf(stderr, "Error when reading %s\n", argv[i]);
		exit(1);
	}
	i++;
	dict = getLexiTreeLeaves(tree);
	if(i<argc && (f = fopen(argv[i], "r"))) {
		character  = getCharacters(f, dict);
		fclose(f);
	} else {
		fprintf(stderr, "Error when reading %s\n", argv[i]);
		exit(1);
	}
	i++;
	nA = argc-i;
	if(nA <= 0)
		exitProg(ErrorArgument, "at least one alignment file is required.");
	model = getEvolutionModelProt();
	best = calibRate(tree, character, model, NULL, CALIB_SYMBOL_PROT, CALIB_VALUE);
	model->freeModel(model);
	model = getEvolutionModelProtStored(tree, best);
	alignData = (TypeAlignmentData*) malloc(nA*sizeof(TypeAlignmentData));
	j=0;
	nT = 0;
	cont = 1;
	fe = fopen(ERROR_NAME, "w");
	if(writeIterF)
		fr = fopen(REPORT_NAME, "w");
	else
		fr = NULL;
	while(cont) {
		pthread_mutex_lock(&mutexN);
		while(i < argc && nT < maxT) {
			pthread_t thread;	
			int ret = 0;
			TypeThreadParameter *param;
			param = (TypeThreadParameter*) malloc(sizeof(TypeThreadParameter));
			param->filename = argv[i];
			param->name = getAlignmentName(param->filename);
			param->threshold = threshold;
			param->dict = dict;
			param->tree = *tree;
			param->model = model;
			param->character = character;
			param->nGam = nGam;
			param->result = &(alignData[j]);
			param->mutex_result = &mutexR;
			param->mutex_number = &mutexN;
			param->mutex_error = &mutexE;
			param->mutex_report = &mutexP;
			param->cond_number = &condN;
			param->fe = fe;
			param->fr = fr;
			param->number = &nT;
			if((ret = pthread_create(&thread, NULL, threadAlignment, (void*) param)) == 0) {
				int err;
				if((err = pthread_detach(thread)) == 0) {
					nT++; i++; j++;
				} else {
					fprintf (stderr, "Error %d while detaching thread: %s\n", err, (char*) strerror(err));
//					pthread_kill(thread, 0);
				}
			} else
				fprintf (stderr, "Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
fprintf(stdout, "alignment %s %d/%d, %d running thread(s)\n", param->name, i-argc+nA, nA, nT);
		}
		cont = (nT > 0);
		if(cont)
			pthread_cond_wait(&condN, &mutexN);
		pthread_mutex_unlock(&mutexN);
	}
	if(fe)
		fclose(fe);
	if(fr)
		fclose(fr);
	for(i=0, size=0; i<nA; i++)
		if(alignData[i].name != NULL)
			alignData[size++] = alignData[i];
	correctBenjaminiHochberg(alignData, size);
	if((f = fopen(outputFileName, "w"))) {
		fprintAlignmentData(f, alignData, size);
		fclose(f);
	} else {
		fprintf(stderr, "Error when writing %s\n", outputFileName);
		exit(1);
	}
	if(extractThreshold != END_DOUBLE) {
		if((f = fopen(EXTRACT_NAME, "w"))) {
			for(i=0; i<size && ((alignData[i].name && alignData[i].pvalueBH <= extractThreshold) || alignData[i].name == NULL); i++)
				if(alignData[i].name)
					fprintf(f, "%s\n", alignData[i].name);
			fclose(f);
		} else {
			fprintf(stderr, "Error when writing %s\n", EXTRACT_NAME);
			exit(1);
		}
	}
	for(i=0; i<size; i++) {
		if(alignData[i].name) {
			free((void*)alignData[i].name);
			if(alignData[i].site)
				free((void*)alignData[i].site);
		}
	}
	free((void*)alignData);
	model->freeModel(model);
	freeLexiTree(dict);
	free(character);
	freeTree(tree);
	time(&t1);
	fprintf(stderr, "Computation time: %.0lf hours %.0lf mins %.0lf secs\n", floor(difftime(t1, t0)/(3600.)), floor(difftime(t1, t0)/(60.))-60.*floor(difftime(t1, t0)/(3600.)), difftime(t1, t0)-60.*floor(difftime(t1, t0)/(60.)));
}

void *threadAlignment(void *data) {
	FILE *fa;
	if((fa = fopen(((TypeThreadParameter*)data)->filename, "r"))) {
		TypeAlignment *align;
		TypeConvergenceData conv;
		align = readAlignment(fa, TYPE_ALIGNMENT_PROTEIN);
		fclose(fa);
		if(align != NULL) {
//			gsl_rng *rg = gsl_rng_alloc(gsl_rng_taus);
			gsl_rng *rg = gsl_rng_alloc(gsl_rng_ranlux389);
			gsl_rng_set(rg, random_seed);
			TypeEvolutionModel *modelC;
			TypeInput input = fixInput(&(((TypeThreadParameter*)data)->tree), ((TypeThreadParameter*)data)->character, align, ((TypeThreadParameter*)data)->dict);
			if(input.tree != &(((TypeThreadParameter*)data)->tree)) {
				TypeEvolutionModel *model;
				double best;
				model = getEvolutionModelProt();
				best = calibRate(input.tree, input.character, model, NULL, CALIB_SYMBOL_PROT, CALIB_VALUE);
				model->freeModel(model);
				modelC = getEvolutionModelProtStored(input.tree, best);
			} else
				modelC = cloneEvolutionModelProtStored(((TypeThreadParameter*)data)->model);
			conv = getAlignmentConvergence(input.tree, input.character, modelC, align, (TypeThreadParameter*)data, rg);
printf("alignment %s\tsize %d\tconv %d\tpvalue %.2le\n", ((TypeThreadParameter*)data)->name, align->size, conv.nSite, conv.pvalue);
			gsl_rng_free(rg);
			if(input.tree != &(((TypeThreadParameter*)data)->tree)) {
				freeTree(input.tree);
				free((void*)input.character);
				modelC->freeModel(modelC);
			} else
				freeCloneEvolutionModelProtStored(modelC);
			pthread_mutex_lock(((TypeThreadParameter*)data)->mutex_result);
				((TypeThreadParameter*)data)->result->name = ((TypeThreadParameter*)data)->name;
				((TypeThreadParameter*)data)->result->size = align->size;		
				((TypeThreadParameter*)data)->result->unknown = ((double)conv.unknown)/((double)(align->size*align->number));
				((TypeThreadParameter*)data)->result->nSite = conv.nSite;
				((TypeThreadParameter*)data)->result->site = conv.site;
				((TypeThreadParameter*)data)->result->pvalue = conv.pvalue;
				((TypeThreadParameter*)data)->result->status = conv.status;
			pthread_mutex_unlock(((TypeThreadParameter*)data)->mutex_result);
			freeAlignment(align);
		} else {
			pthread_mutex_lock(((TypeThreadParameter*)data)->mutex_error);
			fprintf(((TypeThreadParameter*)data)->fe, "Error when reading %s\n", ((TypeThreadParameter*)data)->filename);
			pthread_mutex_unlock(((TypeThreadParameter*)data)->mutex_error);
			((TypeThreadParameter*)data)->result->name = NULL;
		}
	} else {
		pthread_mutex_lock(((TypeThreadParameter*)data)->mutex_error);
		fprintf(((TypeThreadParameter*)data)->fe, "Error when opening %s\n", ((TypeThreadParameter*)data)->filename);
		pthread_mutex_unlock(((TypeThreadParameter*)data)->mutex_error);
		((TypeThreadParameter*)data)->result->name = NULL;
	}
	pthread_mutex_lock(((TypeThreadParameter*)data)->mutex_number);
		(*((TypeThreadParameter*)data)->number)--;
		pthread_cond_signal(((TypeThreadParameter*)data)->cond_number);
	pthread_mutex_unlock(((TypeThreadParameter*)data)->mutex_number);
	free(data);
	return NULL;
}

char *getAlignmentName(char *filename) {
	char *tmp, *name;
	if((tmp=strrchr(filename, '/')) == NULL)
		tmp = filename;
	else
		tmp++;
	name = (char*) malloc((strlen(tmp)+1)*sizeof(char));
	strcpy(name, tmp);
	if((tmp = strrchr(name, '.')) != NULL)
		tmp[0] = '\0';
	return name;
}

TypeConvergenceData getAlignmentConvergence(TypeTree *tree, int *character, TypeEvolutionModel *modelC, TypeAlignment *align, TypeThreadParameter *param, gsl_rng *rg) {
	int i, j, g, a, s, pos, nSim, *column, *occ, *conv, nConv, *al, *N, *eff;
	double alpha, *gamma, *unknown, *ce, *v, val, max;
	TypeEvolutionModel **modelT, *model;
	TypeConvergenceData data;
	TypeAATree *t;
	model = getEvolutionModelProt();
	alpha = estimateGammaParameter(tree, model, rg, align, param->nGam);
	model->freeModel(model);
	gamma = getGamma(alpha, param->nGam);
	column = (int*) malloc(tree->size*sizeof(int));
	modelT = (TypeEvolutionModel**) malloc(param->nGam*sizeof(TypeEvolutionModel*));
	unknown = (double*) malloc(align->number*sizeof(double));
	data.unknown = fillUnknown(align, unknown);
	for(g=0; g<param->nGam; g++)
		modelT[g] = getEvolutionModelProtStored(tree, gamma[g]);
	nConv = 0;
	for(i=0; i<align->number; i++)
		if(character[i])
			nConv++;
	conv = (int*) malloc(nConv*sizeof(int));
	al = (int*) malloc((nConv+1)*sizeof(int));
	for(i=0, j=0; i<align->number; i++)
		if(character[i])
			conv[j++] = i;
	occ = (int*) malloc(modelC->cardinal*sizeof(int));
	for(i=0, j=0; i<modelC->cardinal; i++)
		occ[i] = 0;
	ce = (double*) malloc(align->size*sizeof(double));
	for(pos=0; pos<align->size; pos++) {
		int n;
		for(n=0; n<align->number; n++)
			column[n] = align->sequence[n][pos];
		for(i=0, a=0; i<nConv; i++)
			if(column[conv[i]] >= 0 && column[conv[i]] < modelC->cardinal) {
				occ[column[conv[i]]]++;
				if(occ[column[conv[i]]] == 2)
					al[a++] = column[conv[i]];
			}
		al[a] = END_INT;
		for(i=0; i<nConv; i++)
			if(column[conv[i]] >= 0 && column[conv[i]] < modelC->cardinal)
				occ[column[conv[i]]] = 0;
		if(a > 0)
			ce[pos] = getConvergenceIndexList(tree, character, modelC, column, al);
		else
			ce[pos] = 1.;
	}
	t = newAATree(align->size);
	for(pos=0; pos<align->size; pos++)
		insertAATree(t, ce[pos]);
	v = (double*) malloc(t->size*sizeof(double));
	setTableIndex(v, t);
	eff = (int*) malloc(t->size*sizeof(int));
	for(i=0; i<t->size; i++)
		eff[i] = 0;
	for(pos=0; pos<align->size; pos++) {
		int index = searchAATree(t, ce[pos]);
		if(index == NO_NODE_AATREE) {
			fprintf(stderr, "Error in binary tree of convergences indexes\n");
			exit(1);
		}
		eff[t->node[index].index]++;
	}
	for(i=1; i<t->size; i++)
		eff[i] += eff[i-1];
	N = (int*) malloc(t->size*sizeof(int));
	for(i=0; i<t->size; i++)
		N[i] = 0;
	nSim = 0;
	do {
		double sum;
		int start, end;
		for(s=0; s<nb_sim; s++) {
			drawColumnUnknown(tree, modelT[gsl_rng_uniform_int(rg, param->nGam)], unknown, (void*) rg, column);
			for(i=0, a=0; i<nConv; i++)
				if(column[conv[i]]>=0 && column[conv[i]]<modelC->cardinal) {
					occ[column[conv[i]]]++;
					if(occ[column[conv[i]]] == 2)
						al[a++] = column[conv[i]];
				}
			al[a] = END_INT;
			for(i=0; i<nConv; i++)
				if(column[conv[i]]>=0 && column[conv[i]]<modelC->cardinal)
					occ[column[conv[i]]] = 0;
			if(a > 0) {
				int ind = searchLowerAATree(t, getConvergenceIndexList(tree, character, modelC, column, al));
				if(ind != NO_NODE_AATREE)
					N[t->node[ind].index]++;
			}
		}
		nSim += nb_sim;
		max = gsl_cdf_binomial_P(N[0], param->threshold, nSim);
		val = HUGE_VAL;
		sum = max;
		start = N[0];
/*
end = start;
double thri = param->threshold*((double)nSim)- 4*sqrt(param->threshold*(1.-param->threshold)*((double)nSim));
double thrs = param->threshold*((double)nSim)+ 4*sqrt(param->threshold*(1.-param->threshold)*((double)nSim));
if(end>=thri)
	printf("%d\t%.2le (%d, %d, %d)\n", 0, max, 0, start, N[0]);
*/ 
		for(i=1; i<t->size && sum<=confidence; i++) {
			double tmp;
			end = start+N[i];
			tmp = gsl_cdf_binomial_P(end, param->threshold, nSim)-gsl_cdf_binomial_P(start, param->threshold, nSim);
/*
if(end>=thri && start<=thrs)
	printf("%d\t%.2le\t%.2le (%d, %d, %d)\n", i, v[i-1], tmp, start, end, N[i]);
*/ 
			if(tmp>max) {
				max = tmp;
				val = v[i-1];
			}
			sum += tmp;
			start = end;
		}
/*
printf("nSim %d max %.2lf val %.2lf %.1lf %.1lf [%.1lf, %.1lf] res %.2le (%d, %.2le)\n\n", nSim, max, val, param->threshold*((double)nSim), 4*sqrt(param->threshold*(1.-param->threshold)*((double)nSim)), param->threshold*((double)nSim)- 4*sqrt(param->threshold*(1.-param->threshold)*((double)nSim)), param->threshold*((double)nSim)+4*sqrt(param->threshold*(1.-param->threshold)*((double)nSim)), gsl_cdf_binomial_Q(eff[i-2]-1, param->threshold, align->size), i-2, 1-sum);
*/ 
	} while(max<confidence && gsl_cdf_binomial_Q(eff[i-2]-1, param->threshold, align->size)<sig_thre && nSim<MAX_SIMULATION);
	if(param->fr != NULL) {
		pthread_mutex_lock(param->mutex_report);
			fprintf(param->fr, "%s\t%d\n", param->name, nSim);
		pthread_mutex_unlock(param->mutex_report);
	}
	if(nSim>=MAX_SIMULATION) {
		double sum;
		int start, end;
		FILE *f;
		data.status = 1;
		if(param->fe == NULL)
			f = stderr;
		else
			f = param->fe;
		max = gsl_cdf_binomial_P(N[0], param->threshold, nSim);
		sum = max;
		start = N[0];
		double thri = param->threshold*((double)nSim)- 4*sqrt(param->threshold*(1.-param->threshold)*((double)nSim));
		double thrs = param->threshold*((double)nSim)+ 4*sqrt(param->threshold*(1.-param->threshold)*((double)nSim));
		pthread_mutex_lock(param->mutex_error);
		fprintf(f, "Warning: possible issue in estimating the p-value of %s\n", param->name);
		if(start>=thri) {
			fprintf(f, "with confidence %.4le (threshold %.2le) :\n", max, HUGE_VAL);
			fprintf(f, "%s\t%c\t%d\t%.2le\t%.4le\t%.4le\t%d\t\n", param->name, '*', align->size, 100*((double)data.unknown)/((double)(align->size*align->number)), 1., 0., 0);
		}
		for(i=1; i<t->size; i++) {
			double tmp;
			end = start+N[i];
			tmp = gsl_cdf_binomial_P(end, param->threshold, nSim)-gsl_cdf_binomial_P(start, param->threshold, nSim);
			if(end>=thri && start<=thrs) {
				fprintf(f, "with confidence %.4le (threshold %.2le) :\n", tmp, v[i-1]);
				fprintf(f, "%s\t%c\t%d\t%.2le\t%.4le\t%.4le\t%d", param->name, '*', align->size, 100*((double)data.unknown)/((double)(align->size*align->number)), gsl_cdf_binomial_Q(eff[i-1]-1, param->threshold, align->size), 0., eff[i-1]);
				if(eff[i-1] > 0) {
					for(pos=0; pos<align->size && ce[pos] < v[i-1]; pos++)
					;
					fprintf(f, "\t%d", pos);
					pos++;
					for(; pos<align->size; pos++)
						if(ce[pos] >= v[i-1])
							fprintf(f, ", %d", pos);
				}
				fprintf(f, "\n");
			}
			sum += tmp;
			start = end;
		}
		fprintf(f, "\n");
		pthread_mutex_unlock(param->mutex_error);
	} else
		data.status = 0;
	freeAATree(t);
	for(g=0; g<param->nGam; g++)
		modelT[g]->freeModel(modelT[g]);
	free((void*)modelT);
	free((void*)unknown);
	free((void*)gamma);
	data.site = (int*) malloc(align->size*sizeof(int));
	data.nSite = 0;
    for(pos=0; pos<align->size; pos++)
		if(ce[pos] >= val)
			data.site[data.nSite++] = pos;
	data.pvalue = (data.nSite>0)?gsl_cdf_binomial_Q(data.nSite-1, param->threshold, align->size):1;
	if(data.nSite>0)
		data.site = (int*) realloc((void*)data.site, data.nSite*sizeof(int));
	else {
		free((void*)data.site);
		data.site = NULL;
	}
	free((void*)eff);
	free((void*)N);
	free((void*)v);
	free((void*)ce);
	free((void*)conv);
	free((void*)al);
	free((void*)occ);
	free((void*)column);
	return data;
}

char getStatusChar(int status) {
	switch(status) {
		case 0:
			return ' ';
		case 1:
			return '*';
	}
	return '*';
}

void reindexAlignments(TypeAlignment *alignment, TypeLexiTree *dict) {
	char **newNamesOrder;
	TypeSymbol **newSequencesOrder;
	int i, n, size = 0, sizeDict;
	sizeDict = getSizeLexiTree(dict);
	newNamesOrder = (char **) malloc(sizeDict*sizeof(char *));
	newSequencesOrder = (TypeSymbol **) malloc(sizeDict*sizeof(TypeSymbol *));
	for(i=0; i<alignment->number; i++)
		if((n = findWordLexi(alignment->name[i], dict)) != END_INT) {
			newSequencesOrder[n] = alignment->sequence[i];
			newNamesOrder[n] = alignment->name[i];
			size++;
		} else {
			free(alignment->sequence[i]);
			free(alignment->name[i]);
		}
	if(size != sizeDict) {
		fprintf(stderr, "%d sequences missing in alignment\n", sizeDict-size);
//		exit(1);
	}
	alignment->number = size;
	free(alignment->sequence);
	free(alignment->name);
	alignment->sequence = newSequencesOrder;
	alignment->name = newNamesOrder;
}

int compareAlignmentData(void const *a, void const *b) {
	if(((TypeAlignmentData*)a)->pvalue > ((TypeAlignmentData*)b)->pvalue)
		return 1;
	if(((TypeAlignmentData*)a)->pvalue < ((TypeAlignmentData*)b)->pvalue)
		return -1;
	return 0;
}

void correctBenjaminiHochberg(TypeAlignmentData *alignData, int size) {
	int i;
	qsort(alignData, size, sizeof(TypeAlignmentData), compareAlignmentData);
	for(i=0; i<size; i++)
		alignData[i].pvalueBH = (alignData[i].pvalue*((double)size))/((double)i+1.);
}

void fprintAlignmentData(FILE *f, TypeAlignmentData *alignData, int size) {
	int i, j;
	fprintf(f, "#name\tstatus\tlength\tunknown percent\tpvalue brute\tpvalue corrected\tnumber of convergent sites\tconvergent sites\n");
	for(i=0; i<size; i++) {
		if(alignData[i].name) {
			fprintf(f, "%s\t%c\t%d\t%.2le\t%.4le\t%.4le\t%d", alignData[i].name, getStatusChar(alignData[i].status), alignData[i].size, 100*alignData[i].unknown, alignData[i].pvalue, alignData[i].pvalueBH, alignData[i].nSite);
			if(alignData[i].nSite > 0) {
				fprintf(f, "\t%d", alignData[i].site[0]+1);
				for(j=1; j<alignData[i].nSite; j++)
					fprintf(f, ", %d", alignData[i].site[j]+1);
			}
			fprintf(f, "\n");
		}
	}
}

int fillUnknown(TypeAlignment *align, double *unknown) {
	int n, tot=0;
	for(n=0; n<align->number; n++) {
		int pos;
		unknown[n] = 0.;
		for(pos=0; pos<align->size; pos++)
			if(align->sequence[n][pos] < 0 || align->sequence[n][pos] >= align->cardinal) {
				unknown[n]++;
				tot++;
			}
		unknown[n] /= (double) align->size;
	}
	return tot;
}

TypeInput fixInput(TypeTree *tree,  int *character, TypeAlignment *align, TypeLexiTree *dictTree) {
	int n, nL;
	TypeInput input;
	for(nL=0; tree->node[nL].child == NOSUCH; nL++)
		;
	reindexAlignments(align, dictTree);
//	purgeAlignment(align);
	if(nL == align->number) {
		input.tree = tree;
		input.character = character;
		return input;
	}
	TypeTree *treeTmp;
	TypeLexiTree *dictAlign;
	dictAlign = getDictFromTable(align->name, align->number);
	treeTmp = pruneLeavesFromDict(tree, dictAlign);
	input.tree = fixBinary(treeTmp);
	freeTree(treeTmp);
	input.character = (int*) malloc(align->number);
	for(n=0; n<align->number; n++)
		input.character[n] = character[findWordLexi(align->name[n], dictTree)];
	return input;
}

TypeEvolutionModel *cloneEvolutionModelProtStored(TypeEvolutionModel *model) {
	TypeEvolutionModel *result;
	result = (TypeEvolutionModel *) malloc(sizeof(TypeEvolutionModel));
	*result = *model;
	result->data = malloc(sizeof(TypeProtData));
	*((TypeProtData*)result->data) = *((TypeProtData*)model->data);
	((TypeProtData*)result->data)->user = malloc(sizeof(TypeProtStoredData));
	*((TypeProtStoredData*)((TypeProtData*)result->data)->user) = *((TypeProtStoredData*)((TypeProtData*)model->data)->user);
	return result;
}

void freeCloneEvolutionModelProtStored(TypeEvolutionModel *model) {
	free(((TypeProtData*)model->data)->user);
	free(model->data);
	free((void*) model);
}

