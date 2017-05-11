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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "Utils.h"
#include "Annotation.h"
#include "StatAnnotation.h"

//./enr -o ../data/GO_terms.txt -b ../data/background2.txt resBX.txt GOTestBX.txt
//./enr -o ../data/GO_terms.txt ../data/bis.csv ./bof.txt 
//./enr -o ../data/GO_terms.txt  -b ../data/background2.txt ../data/bis.csv ./bof.txt
#define HELPMESSAGE "\nusage: bonf [options] <annotation file> <partition file> [<output file>]\n\noptions are\n\t-o <ontology file>\t\tload ontology descriptors\n\t-t <real number>\tset the threshold (0.001 as default)\n\t-f c or f\t\t\tindicate the partition format\n\t-h\t\t\tdisplay help\n"

static char **readList(FILE *f);
static int *getIndexFromName(TypeLexiTree *dict, char **name);

int main(int argc, char **argv) {		
	char ontoFileName[STRING_SIZE], outputFileName[STRING_SIZE], option[256], format='s', **name, **background;
	FILE *f;
	int i, j, niter = 100, number = 50, total, *list, nList, nSim = 10000;
	TypeAnnotation *an;
	TypeStatAnnotation *sap;
	TypeSignificantTable *sig;
	TypeOntologyInfo *info;
	TypeLexiTree *dict;
	
	info = NULL;
	background = NULL;
	for(i=0; i<256; i++)
		option[i] = 0;
	   
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['i']) {
			option['i'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &niter) == 1)
				i++;
			else
				exitProg(ErrorArgument, "an integer number is required after option -i");
		}
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &format) == 1)
				i++;
			else
				exitProg(ErrorArgument, "'c' or 'f' are expected after option -f");
		}
		if(option['o']) {
			option['o'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%s", ontoFileName) == 1) {
				i++;
				if((f = fopen(ontoFileName, "r"))) {
					info = readOntologyInfo(f);
					fclose(f);
				} else
					exitProg(ErrorReading, ontoFileName);
			} else
				exitProg(ErrorArgument, "'c' or 'f' are expected after option -f");
		}
		if(option['n']) {
			option['n'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &number) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a float number is required after option -t");
		}
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &total) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a float number is required after option -t");
		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nSim) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a float number is required after option -s");
		}
		if(option['b']) {
			option['b'] = 0;
			if((i+1)<argc && (f = fopen(argv[i+1], "r"))) {
				background = readList(f);
				i++;
			} else
				exitProg(ErrorArgument, "a file name is required after option -b");
		}
		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}
	if(i>=argc)
		exitProg(ErrorArgument, "wrong or absent file name");
	if((f = fopen(argv[i], "r"))) {
		an = readAnnotation(f);
		fclose(f);
	} else {
		fprintf(stderr, "Error while opening %s\n", argv[i]);
		exit(1);
	}
//fprintAnnotation(stdout, an);
	i++;
	if(i>=argc)
		exitProg(ErrorArgument, "wrong or absent file name");
	if((f = fopen(argv[i], "r"))) {
		name = readList(f);
		fclose(f);
	} else {
		fprintf(stderr, "Error while opening %s\n", argv[i]);
		exit(1);
	}
	i++;
	sap = newStatAnnotation(an, background);
	dict = getDictFromTable(sap->nameG, sap->sizeG);
	list = getIndexFromName(dict, name);

	freeLexiTree(dict);
	sig = getSignificantTableList(sap, list);
	for(nList=0; list[nList]!=END_INT; nList++)
		;
	fillCorrectedEmpirical(sap, sig, nList, nSim);
	if(i<argc) {
		sprintf(outputFileName, "%s", argv[i]);
	} else {
		sprintf(outputFileName, "%s", "out.txt");
	}
	if((f = fopen(outputFileName, "w"))) {
		fprintSignificantTable(f, sig, sap->sizeG, sap->nameA, sap->nameG, info);
		fclose(f);
	} else {
		fprintf(stderr, "Error while opening %s\n", outputFileName);
		exit(1);
	}
	free((void*) list);
	for(i=0; name[i] != NULL; i++)
		free((void*) name[i]);
	free((void*) name);
	if(background != NULL) {
		for(i=0; background[i] != NULL; i++)
			free((void*) background[i]);
		free((void*) background);
	}
	freeAnnotation(an);
	freeOntologyInfo(info);
	freeStatAnnotation(sap);
	freeSignificantTable(sig);
	return 0;
}

#define MAX_SIZE_TMP 50
#define INC_BUFFER 50
#define IS_SEP(c) (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == ';' || c == ',')
#define IS_SEP_NO_LINE(c) (c == ' ' || c == '\t' || c == ';' || c == ',')

char **readList(FILE *f) {
	char c, tmp[MAX_SIZE_TMP+1], **list;
	int size, sizeBuffer;

	sizeBuffer = INC_BUFFER;
	list= (char**) malloc(sizeBuffer*sizeof(char*));
	size = 0;
	do {
		c = getc(f);
	} while(c!=EOF && IS_SEP(c)); 
	while(c != EOF) {
		int i;
		i = 0;
		while(i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c)) {
			tmp[i] = c;
			c = getc(f);
			i++;
		}
		tmp[i++] = '\0';
		if(i == MAX_SIZE_TMP) {
			fprintf(stderr, "Ident too long (%s) ...", tmp);
			exit(1);
		}
		if(i>1) {
			if(size >= sizeBuffer) {
				sizeBuffer += INC_BUFFER;
				list = (char**) realloc((void *) list, sizeBuffer*sizeof(char*));
			}
			list[size] = (char*) malloc((strlen(tmp)+1)*sizeof(char));
			strcpy(list[size], tmp);
			size++;
		}
		while(c!=EOF && IS_SEP(c))
			c=getc(f);
	}
	if(size >= sizeBuffer) {
		sizeBuffer += INC_BUFFER;
		list = (char**) realloc((void *) list, sizeBuffer*sizeof(char*));
	}
	list[size++] = NULL;
	return list;
}

int *getIndexFromName(TypeLexiTree *dict, char **name) {
	int i, ind, size, *list;
	TypeLexiTree *dl;
	dl = newLexiTree();
	if(name == NULL)
		return NULL;
	size = 0;
	for(i=0; name[i] != NULL; i++)
		size++;
	list = (int*) malloc((size+1)*sizeof(int));
	size = 0;
	for(i=0; name[i] != NULL; i++)
		if((ind = findWordLexi(name[i], dict)) != END_INT) {
			if(findWordLexi(name[i], dl) == END_INT) {
				list[size++] = ind;
				addWordLexi(name[i], i, dl);
			} else
				fprintf(stderr, "Warning - Duplicate item %s\n", name[i]);
		} else
			fprintf(stderr, "Warning - Unknown item %s\n", name[i]);
	freeLexiTree(dl);
	list[size++] = END_INT;
	list = (int*) realloc((void*)list, size*sizeof(int));
	return list;
}
