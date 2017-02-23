#include "utils.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "circlePit.h"


PtrCirclePit createCirclePit(int maxSize) {
  PtrCirclePit p_cp;

  p_cp           = (PtrCirclePit) malloc(sizeof(CirclePit));
  p_cp->sizeData = CIRCLEPIT_BEGINSIZE;
  p_cp->data     = (void**) malloc(p_cp->sizeData*sizeof(void*));
  p_cp->begin    = 0;
  p_cp->nbElem   = 0;

  return p_cp;
}

void freeCirclePit(PtrCirclePit p_cp) {
  free(p_cp->data);
  free(p_cp);
}

void  push(PtrCirclePit p_cp, void* elem) {

  if (p_cp->nbElem==p_cp->sizeData) {

    /*
      Il faut augmenter la taille
     */
    /* 
       1. allouer de la memoire : on double l'espace
     */
    p_cp->data = (void**) realloc(p_cp->data, (p_cp->sizeData*2)*sizeof(void*));
    /* 
       2. recopier les elements avant begin
     */
    DEBUGP("deplacement mem %p <- %p, %lu octets\n", (void*) p_cp->data+(p_cp->sizeData*sizeof(void*)), (const void *) p_cp->data, (size_t) p_cp->begin*sizeof(void*));
    memcpy((void*) p_cp->data+(p_cp->sizeData*sizeof(void*)), (const void *) p_cp->data, (size_t) p_cp->begin*sizeof(void*));

    /* 
       3. mettre a jour la structure
     */
    p_cp->sizeData = p_cp->sizeData *2;
    DEBUGP("Taille augmentee a : %d\n", p_cp->sizeData);
  }

  p_cp->data[(p_cp->begin+p_cp->nbElem)%p_cp->sizeData] = elem;
  p_cp->nbElem++;

}

void* first (PtrCirclePit p_cp){

  if (p_cp->nbElem==0) {
    ERROR(2, "circlePit is empty\n");
  }

  return p_cp->data[p_cp->begin];
}

void* get (PtrCirclePit p_cp, int i) {
  if (p_cp->nbElem<=i) {
	ERROR(2,"circlePit has not enough element\n");
  }
  return p_cp->data[(p_cp->begin+i)%p_cp->sizeData];
}


void* last (PtrCirclePit p_cp){

  return get(p_cp, p_cp->nbElem-1);
/*
  if (p_cp->nbElem==0) {
    ERROR("circlePit is empty\n",2);
  }

  return p_cp->data[(p_cp->begin+p_cp->nbElem-1)%p_cp->sizeData];
*/
}

void* pop (PtrCirclePit p_cp) {
  
  int oldBegin;

  if (p_cp->nbElem==0) {
    ERROR(2, "circlePit is empty\n");
  }

  p_cp->nbElem--;
  oldBegin = p_cp->begin;
  p_cp->begin = (p_cp->begin+1)%p_cp->sizeData;

  return p_cp->data[oldBegin];
}

int   isEmptyCirclePit(PtrCirclePit p_cp) {
  return (p_cp->nbElem==0);
}

void  printfCirclePit(FILE* f, PtrCirclePit p_cp) {
  fprintf(f, "circlePit at %p\n", p_cp);
  fprintf(f, "  sizeData %d\n", p_cp->sizeData);
  fprintf(f, "  begin    %d\n", p_cp->begin);
  fprintf(f, "  nbElem   %d\n", p_cp->nbElem);
}
