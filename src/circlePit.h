#ifndef _H_CirclePit 
#define _H_CirclePit

#include <stdio.h>

#define CIRCLEPIT_BEGINSIZE 100

typedef struct s_CirclePit {
  
  int    sizeData;
  int    nbElem;
  int    begin;
  void** data;
} CirclePit, *PtrCirclePit;


PtrCirclePit createCirclePit();
void         freeCirclePit(PtrCirclePit p_cp);

void  push   (PtrCirclePit p_cp, void* elem);
void* first  (PtrCirclePit p_cp);
void* get    (PtrCirclePit p_cp, int i);
void* last   (PtrCirclePit p_cp);
void* pop    (PtrCirclePit p_cp);
int   isEmptyCirclePit(PtrCirclePit p_cp);
void  printfCirclePit(FILE* f, PtrCirclePit p_cp);

#endif
