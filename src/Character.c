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
#include "Character.h"

#define MAX_NAME_SIZE 300
#define INC_BUFFER 100

int *getCharacters(FILE *characterFile, TypeLexiTree *dict) {
	TypeCharacters *charactersTab = NULL;
	int i, n, sizeDict;
	int *characters;
	sizeDict = getSizeLexiTree(dict);
	characters = (int *) malloc(sizeDict*sizeof(int));
	charactersTab = readCharacterTab(characterFile);
	for(i=0; i<sizeDict; i++)
		characters[i] = 0;
	for(i=0; i<charactersTab->size; i++) {
		if((n=findWordLexi(charactersTab->name[i], dict)) != END_INT)
			characters[n] = charactersTab->values[i]; 
		else {
			fprintf (stderr, "%s: non defini dans le fichier de character\n", charactersTab->name[i]);
			exit(1);
		}
	}
	freeCharacter(charactersTab);
	return characters;
}

void freeCharacter(TypeCharacters *caracter) {
	if (caracter == NULL)
		return;
	if (caracter->name != NULL) {
		int i;
		for (i = 0 ; i < caracter->size ; i++)
			free(caracter->name[i]);
		free(caracter->name);
	}
	if(caracter->values != NULL)
		free(caracter->values);
	free(caracter);
}

TypeCharacters *readCharacterTab(FILE *f) {
	int sizeBuf;
	char c;
	TypeCharacters *res;
	
	res = (TypeCharacters*) malloc(sizeof(TypeCharacters));
	sizeBuf = INC_BUFFER;
	res->name = (char**) malloc(sizeBuf*sizeof(char*));
	res->values = (int *) malloc(sizeBuf*sizeof(int));
	res->size = 0;
	do {
		char *tmp;
		int i;
		tmp = (char*) malloc((MAX_NAME_SIZE+1)*sizeof(char));
		for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f));
		if(c == '\'' || c == '"') {
			c = fgetc(f);
			for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
				tmp[i] = c;
				c = fgetc(f);
			}
			if(c == '\'' || c == '"')
				c = fgetc(f);
		} else {
			for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issep(c); i++) {
				tmp[i] = c;
				c = fgetc(f);
			}
		}
		if(i == MAX_NAME_SIZE){
			fprintf(stderr, "Name too much long while reading a character file...");
			exit(1);
		}
		tmp[i++] = '\0';
		if(i>1) {
			char bof[MAX_NAME_SIZE+1];
			if(res->size >= sizeBuf) {
				sizeBuf += INC_BUFFER;
				res->name   = (char**) realloc((void*) res->name  , sizeBuf * sizeof(char*));
				res->values = (int *)  realloc((void*) res->values, sizeBuf * sizeof(int));
			}
			res->name[res->size] = (char *) realloc((void*) tmp, i*sizeof(char));
			for(; c != EOF && issep(c); c = fgetc(f));
			while(c != EOF && !isline(c)) {
				for(i=0; c != EOF && !issepline(c) && i<MAX_NAME_SIZE; i++) {
					bof[i] = c;
					c = fgetc(f);
				}
				bof[i++] = '\0';
				res->values[res->size] = atoi(bof);
				for(; c != EOF && issep(c); c = fgetc(f));
			}
			res->size++;
		} else
			free((void*) tmp);
	} 
	while(c != EOF);
	if(res->size > 0)
		res->values = (int*) realloc((int*) res->values, (res->size)*sizeof(int));
	else {
		free((void*)res->values);
		res->values = NULL;
	}
	return res;
}
