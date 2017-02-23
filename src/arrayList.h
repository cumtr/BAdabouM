#ifndef _H_ArrayList 
#define _H_ArrayList

#include <stdio.h>
#include "utils.h"

#define ARRAYLIST_INCRSIZE 10

#define AL_compareFElemElem int (*AL_compareFuncEE)(const PtrVoid,\
                                                    const PtrVoid)


typedef struct s_arrayList {
	AL_compareFElemElem;
    int nbElem;
    int nbMaxElem;
    PtrVoid* data;
} ArrayList, *PtrArrayList;

PtrArrayList           createArrayList(AL_compareFElemElem);
void                   freeArrayList(PtrArrayList p_al);

void    AL_addLast  (PtrArrayList p_al, PtrVoid elem);
void    AL_addFirst (PtrArrayList p_al, PtrVoid elem);
void    AL_addSorted(PtrArrayList p_al, PtrVoid elem);

void    AL_insert   (PtrArrayList p_al, PtrVoid elem, int index);
PtrVoid AL_popFirst (PtrArrayList p_al);
PtrVoid AL_popLast  (PtrArrayList p_al);
PtrVoid AL_pop      (PtrArrayList p_al, int index);
int     AL_isEmpty  (PtrArrayList p_al);
void    AL_printf   (FILE* f, PtrArrayList p_al);

#endif
