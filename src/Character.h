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




#ifndef CharacterF
#define CharacterF

#include <stdio.h>
#include "Utils.h"

typedef struct NAME_CHARACTER_TAB {
	char **name;
	int size, *values;
} TypeCharacters;

TypeCharacters *readCharacterTab(FILE *f);
int *getCharacters(FILE *characterFile, TypeLexiTree *dict);
void freeCharacter(TypeCharacters *caracter);

#endif
