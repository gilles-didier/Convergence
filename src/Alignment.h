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




#ifndef AlignmentF
#define AlignmentF

#include <stdlib.h>
#include <stdio.h>

#define EMPTY 255
#define EMPTY_CAR '-'
#define WHATEVER_CAR_DNA 'N'
#define WHATEVER_CAR_PROT 'X'
#define TYPE_ALIGNMENT_PROTEIN  'p'
#define TYPE_ALIGNMENT_DNA  'd'
#define TYPE_ALIGNMENT_RNA 'r'
#define DNA "ACGT" /*"ACGT" + "YRMKWSBDHVN"*/
#define RNA "ACGUYRMKWSBDHVN" /*"ACGU" + "YRMKWSBDHVN"*/
#define PRO "ARDNCEQGHILKMFPSTWYV"

typedef int TypePosition;
typedef int TypeSymbol;
typedef int TypeNumber;

typedef enum TAF {
	fasta=0,
	clustal,
	msf,
	markx,
	srs,
	unknown
} TypeAlignmentFile;

typedef struct ALIGNMENT {
	TypeNumber number;
	TypePosition size;
	TypeSymbol **sequence, empty, cardinal, whatever;
	char **name, *table;
} TypeAlignment;

TypeAlignment *readAlignment(FILE *f, char typeAlphabet);
TypeAlignment *readAlignmentFasta(FILE *f, char *table, int canInc);
TypeAlignment *readAlignmentMsf(FILE *f, char *table, int canInc);
TypeAlignment *readAlignmentClustal(FILE *f, char *table, int canInc);
void printAlignmentFasta(FILE *f, TypeAlignment *al, int sizeLine);
void printAlignmentMsf(FILE *f, TypeAlignment *al, int sizeLine);
void printAlignmentSrs(FILE *f, TypeAlignment *al, int sizeLine);
void printAlignmentMarkX(FILE *f, TypeAlignment *al, int sizeLine);
void printAlignmentTex(FILE *f, TypeAlignment *al, int sizeLine);
void printAlignmentTexBis(FILE *f, TypeAlignment *al, int sizeLine);
void printHeadPair(FILE *f, TypeAlignment *al);
void printHeadMulti(FILE *f, TypeAlignment *al);
void freeAlignment(TypeAlignment *alignment);
void purgeAlignment(TypeAlignment *al);

#endif
