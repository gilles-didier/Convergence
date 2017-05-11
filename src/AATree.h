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




#ifndef AATreeF
#define AATreeF

#include <stdio.h>
#include <limits.h>

#define NO_NODE_AATREE INT_MAX

typedef struct AA_NODE {
   double val;
   int lv, l, r, index;
} TypeAANode;

typedef struct {
	int root, nullnode, deleted, last, size, buffer;
	TypeAANode *node;
} TypeAATree;

// new empty AA tree
TypeAATree *newAATree(int initSize);
// delete the AA tree(release all resources)
void freeAATree(TypeAATree *t);
// search an item in the tree by the key, and return the value
int searchAATree(TypeAATree *t, double val);
int searchLowerAATree(TypeAATree *t, double val);
int searchUpperAATree(TypeAATree *t, double val);
// insert an item by key-value into the tree
void insertAATree(TypeAATree *t, const double val);
void fprintAATree(FILE *f, TypeAATree *t);
void setTableIndex(double *tab, TypeAATree *t);

#endif
