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




#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "AATree.h"

#define SIZE_BUFFER 100

static int skewAATree(TypeAATree *t, int n);
static int splitAATree(TypeAATree *t, int n);
static int searchAATreeRec(TypeAATree *t,  double val, int n);
static int searchLowerAATreeRec(TypeAATree *t, double val, int last, int n);
static int searchUpperAATreeRec(TypeAATree *t, double val, int last, int n);
static int insertAATreeRec(TypeAATree *t, double val, int n);
static void setTableIndexRec(double *tab, TypeAATree *t, int *cur, int n);


TypeAATree *newAATree(int initSize) {
	TypeAATree *t = (TypeAATree*)malloc(sizeof(TypeAATree));
	if(t == NULL)
		return NULL;
	if(initSize <= 0)
		initSize = SIZE_BUFFER;
	t->buffer = initSize;
	t->node = (TypeAANode*) malloc((t->buffer+1)*sizeof(TypeAANode));
	if(t->node == NULL) {
		free(t);
		return NULL;
	}
	t->node[0].lv = 0;
	t->node[0].l = -1;
	t->node[0].r = -1;
	t->node[0].index = NO_NODE_AATREE;
	t->node++;
	t->nullnode = -1;
	t->root = -1;
	t->deleted = NO_NODE_AATREE;
	t->size = 0;
	return t;
}

void freeAATree(TypeAATree *t) {
	assert(t && t->node);
	t->node--;
	free((void*) t->node);
	free((void*) t);
}

int skewAATree(TypeAATree *t, int n) {
	assert(t && t->node && n<t->size);
	if(t->node[n].lv != t->node[t->node[n].l].lv)
		return n;
	int left = t->node[n].l;
	t->node[n].l = t->node[left].r;
	t->node[left].r = n;
	return left;
}

int splitAATree(TypeAATree *t, int n) {
	assert(t && t->node && n<t->size);
	if(t->node[t->node[t->node[n].r].r].lv != t->node[n].lv)
		return n;
	int right = t->node[n].r;
	t->node[n].r = t->node[right].l;
	t->node[right].l = n;
	t->node[right].lv++;
	return right;
}

int searchAATreeRec(TypeAATree *t, double val, int n) {
	if(n == t->nullnode)
		return NO_NODE_AATREE;
	else {
		if(t->node[n].val > val)
			return searchAATreeRec(t, val, t->node[n].l);
		else {
			if(t->node[n].val < val)
				return searchAATreeRec(t, val, t->node[n].r);
			else
				return n;
		}
	}
}

int searchAATree(TypeAATree *t, double val) {
   return searchAATreeRec(t, val, t->root);
}

int searchLowerAATreeRec(TypeAATree *t, double val, int last, int n) {
	if(n == t->nullnode)
		return last;
	else {
		if(t->node[n].val > val)
			return searchLowerAATreeRec(t, val, last, t->node[n].l);
		else {
			if(t->node[n].val < val)
				return searchLowerAATreeRec(t, val, n, t->node[n].r);
			else
				return n;
		}
	}
}

int searchLowerAATree(TypeAATree *t, double val) {
   return searchLowerAATreeRec(t, val, NO_NODE_AATREE, t->root);
}

int searchUpperAATreeRec(TypeAATree *t, double val, int last, int n) {
	if(n == t->nullnode)
		return last;
	else {
		if(t->node[n].val > val)
			return searchUpperAATreeRec(t, val, n, t->node[n].l);
		else {
			if(t->node[n].val < val)
				return searchUpperAATreeRec(t, val, last, t->node[n].r);
			else
				return n;
		}
	}
}

int searchUpperAATree(TypeAATree *t, double val) {
   return searchUpperAATreeRec(t, val, NO_NODE_AATREE, t->root);
}

int insertAATreeRec(TypeAATree *t, const double val, int n) {
	if(n == t->nullnode) {
		if(t->size >= t->buffer) {
			t->buffer += SIZE_BUFFER;
			t->node--;
			t->node = (TypeAANode*) realloc((void*)t->node, t->buffer*sizeof(TypeAANode));
			t->node++;
		}
		t->node[t->size].val = val;
		t->node[t->size].lv = 1;
		t->node[t->size].index = NO_NODE_AATREE;
		t->node[t->size].l = t->nullnode;
		t->node[t->size].r = t->nullnode;
		t->size++;
		return t->size-1;
	}
	if(val < t->node[n].val)
		t->node[n].l = insertAATreeRec(t, val, t->node[n].l);
	else {
		if(val > t->node[n].val)
			t->node[n].r = insertAATreeRec(t, val, t->node[n].r);
		else
			return n;
	}
	n = skewAATree(t, n);
	n = splitAATree(t, n);
	return n;
}

void insertAATree(TypeAATree *t, const double val) {
	t->root = insertAATreeRec(t, val, t->root);
}


void fprintAATree(FILE *f, TypeAATree *t) {
	int i;
	fprintf(f, "size %d\troot %d\n", t->size, t->root);
	for(i=0; i<t->size; i++)
		fprintf(f, "node %d/%.2le l %d r %d (%d)\n", i, t->node[i].val, t->node[i].l, t->node[i].r, t->node[i].index);
}

void setTableIndexRec(double *tab, TypeAATree *t, int *index, int n) {
	if(n == t->nullnode)
		return;
	setTableIndexRec(tab, t, index, t->node[n].r);
	t->node[n].index = *index;
	tab[(*index)++] = t->node[n].val;
	setTableIndexRec(tab, t, index, t->node[n].l);
}

void setTableIndex(double *tab, TypeAATree *t) {
	int index = 0;
	setTableIndexRec(tab, t, &index, t->root);
}
/*
int removeAATreeRec(TypeAATree *t, const double val, int n) {
	if(n == t->nullnode)
		return n;
	t->last = n;
	if(val < t->node[n].val)
		t->node[n].l = removeAATreeRec(t, val, t->node[n].l);
	else {
		t->deleted = n;
		t->node[n].r = removeAATreeRec(t, val, t->node[n].r);
	}
	if(n == t->last && t->deleted != t->nullnode && val == t->node[t->deleted]val) {
		t->node[t->deleted].val = t->node[n].val;
		t->deleted = t->nullnode;
		n = t->node[n].r;
		free(t->last->val);
		free(t->last);
		t->num_entries--;
	} else if(t->node[n].l->lv < t->node[n].lv-1 || t->node[n].r->lv < t->node[n].lv-1) {
	t->node[n].lv--;
	if(t->node[n].r->lv > t->node[n].lv) t->node[n].r->lv = t->node[n].lv;
	n = aa_skew(n);
	t->node[n].r = aa_skew(t->node[n].r);
	t->node[n].r->r = aa_skew(t->node[n].r->r);
	n = aa_split(n);
	t->node[n].r = aa_split(t->node[n].r);
	}
	return n;
}

void removeAATree(TypeAATree *t, const double val) {
   t->root = removeAATreeRec(t, val, t->root);
}

static void aa_foreach_(const aanode *nn, const aanode *n, AAForeach cb, void *arg) {
   if(n == nn) return;
   (*cb)(n, arg);
   aa_foreach_(nn, t->node[n].l, cb, arg);
   aa_foreach_(nn, t->node[n].r, cb, arg);
}
void aa_foreach(const aat *t, AAForeach callback, void *arg) {
   aa_foreach_(t->nullnode, t->root, callback, arg);
}

static void aa_map_(const aanode *nn, aanode *n, AAMap cb) {
   if(n == nn) return;
   void *tmp = (*cb)(n);
   free(t->node[n].val);
   t->node[n].val = tmp;
   aa_map_(nn, t->node[n].l, cb);
   aa_map_(nn, t->node[n].r, cb);
}
void aa_map(const aat *t, AAMap callback) {
   aa_map_(t->nullnode, t->root, callback);
}

static aanode* aa_first_(aanode *nn, aanode *n) {
   if(n == nn) return nn;
   else if(t->node[n].l == nn) return n;
   else return aa_first_(nn, t->node[n].l);
}
aanode* aa_first(const aat *t) {
   return aa_first_(t->nullnode, t->root);
}

static aanode* aa_last_(aanode *nn, aanode *n) {
   if(n == nn) return nn;
   else if(t->node[n].r == nn) return n;
   else return aa_last_(nn, t->node[n].r);
}
aanode* aa_last(const aat *t) {
   return aa_last_(t->nullnode, t->root);
}
*/
