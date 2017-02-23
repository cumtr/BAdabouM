#include <stdlib.h>
#include <string.h>
#include "arrayList.h"


PtrArrayList createArrayList(AL_compareFElemElem){
	PtrArrayList p_al;

	p_al = (PtrArrayList) malloc(sizeof(ArrayList));

	p_al->AL_compareFuncEE = AL_compareFuncEE;
	p_al->nbElem = 0;
	p_al->nbMaxElem = 0;
	p_al->data = NULL;

	return p_al;
}

void freeArrayList(PtrArrayList p_al) {

	free(p_al->data);
	free(p_al);
}

void __AL_checkSize(PtrArrayList p_al) {
	if (p_al->nbElem==p_al->nbMaxElem) {
		p_al->data = realloc(p_al->data, (ARRAYLIST_INCRSIZE + p_al->nbElem) * sizeof(PtrVoid));
		p_al->nbMaxElem+=ARRAYLIST_INCRSIZE;
	}
}

void __AL_swap(PtrVoid* p1, PtrVoid* p2) {
    PtrVoid tmp;
	tmp = *p1;
	*p1  = *p2;
	*p2  = tmp;
}

void AL_addSorted(PtrArrayList p_al, PtrVoid elem) {

	AL_addLast(p_al, elem);

	int i = p_al->nbElem-1;
	while(i>0 && p_al->AL_compareFuncEE(p_al->data[i-1], p_al->data[i])>0) {
		__AL_swap(&(p_al->data[i-1]), &(p_al->data[i]));
		i--;
	}
}

void AL_insert(PtrArrayList p_al, PtrVoid elem, int index) {
	__AL_checkSize(p_al);
	if (p_al->nbElem>0)
		memmove(p_al->data + (index+1),
				p_al->data + index,
				(p_al->nbElem-index)*sizeof(PtrVoid));
	p_al->data[index] = elem;
	p_al->nbElem++;
}

void AL_addLast(PtrArrayList p_al, PtrVoid elem) {
	AL_insert(p_al, elem, p_al->nbElem);
}

void AL_addFirst(PtrArrayList p_al, PtrVoid elem) {
	AL_insert(p_al, elem, 0);
}


PtrVoid AL_pop(PtrArrayList p_al, int index) {
    PtrVoid ret;
	ret = p_al->data[index];
	if (index<(p_al->nbElem-1))
		memmove(p_al->data + index,
				p_al->data + (index+1),
				(p_al->nbElem-index-1)*sizeof(PtrVoid));
	p_al->nbElem--;
	return ret;
}

PtrVoid AL_popLast (PtrArrayList p_al) {
	return AL_pop(p_al, p_al->nbElem-1);
}

PtrVoid AL_popFirst(PtrArrayList p_al) {
	return AL_pop(p_al, 0);
}

int AL_isEmpty(PtrArrayList p_al) {
	return p_al->nbElem==0;

}
void AL_printf(FILE* f, PtrArrayList p_al) {
	fprintf(f, "ArrayList at %p with %d elem (max %d)\n", p_al, p_al->nbElem, p_al->nbMaxElem);
}


/*
typedef struct {
	int nb;
	char name[50];
} Toto, *PtrToto;

void fprintf_toto(FILE* f, Toto p) {
	fprintf(f, "%s %d\n",p.name, p.nb);
}

int compar_toto(const void * p1, const void *p2) {
	PtrToto p_t1 = (PtrToto) p1;
	PtrToto p_t2 = (PtrToto) p2;

	return p_t2->nb - p_t1->nb;
}

int AL__test() {

	PtrArrayList p_al = createArrayList(NULL);
	int i;

	Toto fredo;
	fredo.nb=10;
	strcpy(fredo.name, "Fredo");

	Toto marta;
	marta.nb=256;
	strcpy(marta.name, "Marta");

	Toto choupinou;
	choupinou.nb=25;
	strcpy(choupinou.name, "Choupinou");

	Toto aurelie;
	aurelie.nb=30;
	strcpy(aurelie.name, "Aurelie");

	PtrToto p_toto;

	AL_addLast(p_al, (PtrVoid) &fredo);
	AL_printf(stdout, p_al);
	for(i=0;i<p_al->nbElem;i++) {
		p_toto =(PtrToto)(p_al->data[i]);
		fprintf(stdout, "%p\n", p_toto);
		fprintf_toto(stdout, *p_toto);
	}

	AL_addFirst(p_al, (PtrVoid) &marta);
	AL_printf(stdout, p_al);
	for(i=0;i<p_al->nbElem;i++) {
		p_toto =(PtrToto)(p_al->data[i]);
		fprintf(stdout, "%p\n", p_toto);

		fprintf_toto(stdout, *p_toto);
	}

	AL_insert(p_al, (PtrVoid) &choupinou, 1);
	AL_printf(stdout, p_al);
	for(i=0;i<p_al->nbElem;i++) {
		p_toto =(PtrToto)(p_al->data[i]);
		fprintf_toto(stdout, *p_toto);
	}

	AL_insert(p_al, (PtrVoid) &aurelie, 1);
	AL_printf(stdout, p_al);
	for(i=0;i<p_al->nbElem;i++) {
		p_toto = (PtrToto)(p_al->data[i]);
		fprintf_toto(stdout, *p_toto);
	}

	fprintf(stdout,"\n");
	fredo.nb=100;
	fprintf_toto(stdout, fredo);

	fprintf(stdout,"\n");
	p_toto = ((PtrToto) AL_popFirst(p_al));
	fprintf_toto(stdout, *p_toto);

	fprintf(stdout,"\n");
	AL_printf(stdout, p_al);
	for(i=0;i<p_al->nbElem;i++) {
		p_toto =(PtrToto)(p_al->data[i]);
		fprintf_toto(stdout, *p_toto);
	}

	while (! AL_isEmpty(p_al)) {
		p_toto = (PtrToto)AL_popLast(p_al);
		fprintf(stdout, "Pop\n");
		fprintf_toto(stdout, *p_toto);
	}
	freeArrayList(p_al);

	fprintf(stdout,"\n***\n");

	p_al = createArrayList(compar_toto);

	AL_addSorted(p_al, (PtrVoid) &choupinou);
	AL_addSorted(p_al, (PtrVoid) &aurelie);
	AL_addSorted(p_al, (PtrVoid) &marta);

	AL_printf(stdout, p_al);
	for(i=0;i<p_al->nbElem;i++) {
		p_toto =(PtrToto)(p_al->data[i]);
		fprintf_toto(stdout, *p_toto);
	}

	fprintf(stdout,"\n***\n");

	while (! AL_isEmpty(p_al)) {
		p_toto = (PtrToto)AL_popLast(p_al);
		fprintf_toto(stdout, *p_toto);
	}
	freeArrayList(p_al);

	return 1;
}

*/
