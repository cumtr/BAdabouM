#include "hashTable.h"
#include <stdlib.h>
#include <string.h>

typedef struct s_HashTableCell {
	char* key;
    PtrVoid elem;
} HashTableCell, *PtrHashTableCell;


unsigned long
djb2(unsigned char *str)
{
    unsigned long hash = 5381;
    int c;

    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}

unsigned long __HT_hash(char *name)
{
	return djb2(name);
}

int HT_compareFuncCellCell (const PtrVoid p_cell1, const PtrVoid p_cell2) {
	/*
	 * fprintf(stderr, "Comparing %s with %s = %d\n", ((PtrHashTableCell) p_cell1)->key, ((PtrHashTableCell)p_cell2)->key, strcmp(((PtrHashTableCell) p_cell1)->key, ((PtrHashTableCell)p_cell2)->key));
	 */
	return strcmp(((PtrHashTableCell) p_cell1)->key, ((PtrHashTableCell)p_cell2)->key);
}


PtrHashTable createHashTable(HT_FreeElem) {
    PtrHashTable p_ht;
    int i;
    p_ht = (PtrHashTable) malloc(sizeof(HashTable));

    p_ht->HT_FreeElemFunc  = HT_FreeElemFunc;
    
    p_ht->nbElem=0;
    p_ht->nbOps =0;

    p_ht->size=HASHTABLE_INITSIZE;
    
    p_ht->data = (ArrayList*) malloc(HASHTABLE_INITSIZE*sizeof(ArrayList));
    
    for(i=0;i<p_ht->size;i++) {
    	p_ht->data[i].AL_compareFuncEE = HT_compareFuncCellCell;
        p_ht->data[i].nbElem=0;
        p_ht->data[i].nbMaxElem=0;

        p_ht->data[i].data=NULL;
    }
    
    return p_ht;
}

void freeHashTable(PtrHashTable p_ht) {
	int i, j;
    
    for(i=0;i<p_ht->size;i++) {
    	ArrayList al = p_ht->data[i];
    	for(j=0;j<al.nbElem;j++) {

    		/* the key and data stored */
    		PtrHashTableCell p_cell = al.data[j];
    		free(p_cell->key);
    		if (p_ht->HT_FreeElemFunc)
    			p_ht->HT_FreeElemFunc(p_cell->elem);
    		/* and the cell */
    		free(p_cell);
     	}
    	/* freeing the tab of cells */
    	free(al.data);
    }
    free(p_ht->data);
    free(p_ht);
}


void emptyingHashTable(PtrHashTable p_ht) {
	int i, j;

	for(i=0;i<p_ht->size;i++) {
    	ArrayList al = p_ht->data[i];
    	for(j=0;j<al.nbElem;j++) {

    		/* the key and data stored */
    		PtrHashTableCell p_cell = al.data[j];
    		free(p_cell->key);
    		if (p_ht->HT_FreeElemFunc)
    			p_ht->HT_FreeElemFunc(p_cell->elem);
    		/* and the cell */
    		free(p_cell);
     	}
		(&p_ht->data[i])->nbElem = 0;
    }
    p_ht->nbElem = 0;
    p_ht->nbOps++;
}

void  HT_add    (PtrHashTable p_ht, char* key, PtrVoid elem) {
    PtrHashTableCell p_cell=NULL;
    if (((p_ht->nbElem+1)/p_ht->size)>HASHTABLE_LIMIT) {
        /* TODO: should rehash here */
    	HT_rehash (p_ht, p_ht->size*HASHTABLE_INCRSIZE);
    }
    p_cell = (PtrHashTableCell) malloc(sizeof(HashTableCell));
    /* insert */
    p_cell->key  = (char*) malloc((strlen(key)+1)*sizeof(char));
    strcpy(p_cell->key, key);

    p_cell->elem = elem;
    unsigned long hash = __HT_hash(p_cell->key);
    unsigned long indexAL = hash%p_ht->size;
    AL_addSorted(&(p_ht->data[indexAL]), p_cell);
    p_ht->nbElem++;
    p_ht->nbOps++;

}

int HT_compareFuncKeyCell (const PtrVoid key, const PtrVoid* p_p_cell) {
	/*fprintf(stderr,"compare %s et %s\n", (char *) key, ((PtrHashTableCell) (*p_p_cell))->key);*/
	return strcmp((char *) key, ((PtrHashTableCell) (*p_p_cell))->key);
}

int __HT_find(PtrHashTable p_ht, int al_index, char* key) {
    /* Given the index of the ArrayList to be searched and
     * the key that is searched
     * search the element key in the *sorted* ArrayList
     * with the bsearch method
     */
	PtrHashTableCell p_p_cell;
	ArrayList al = p_ht->data[al_index];
	/*
	 * The compar routine is expected to have two arguments which point to the key object and to an array member, in that order.
	 */
	p_p_cell = al.data[0];

	/*
	 * fprintf(stderr, "bsearch(key=%s, p=%p, nb=%d)\n", key, al.data, al.nbElem);
	 */

	p_p_cell = (PtrHashTableCell) bsearch(key, al.data, al.nbElem, sizeof(PtrVoid), (int (*) (const PtrVoid, const PtrVoid)) HT_compareFuncKeyCell);
    
    if (p_p_cell)
    	return ((PtrVoid*) p_p_cell - (PtrVoid*) al.data);
    else
    	return -1;
}


PtrVoid HT_get(PtrHashTable p_ht, char* key, int pop) {
    PtrHashTableCell p_cell;
    PtrVoid elem = NULL;
    unsigned long index, indexAL;
    unsigned long hash = __HT_hash(key);
    indexAL = hash%p_ht->size;

    //fprintf(stderr, "***\n");

    index    = __HT_find(p_ht, indexAL, key);

    if (index!=-1) {
    	p_cell = p_ht->data[indexAL].data[index];
    	elem = p_cell->elem;
    	if (pop) {
    		/* Have to:
    		 * 1. Remove the cell in the ArrayList
    		 * 2. Move the other cells to ensure ArrayList is sorted
    		 * 3. free cell
    		 */
    		 if ((p_ht->data[indexAL].nbElem-index)>1) {
    			 fprintf(stderr,"HT_GET index = %d (nbElem = %d) memmove(%p, %p, %d);\n", index,p_ht->data[indexAL].nbElem, p_ht->data[indexAL].data + index, p_ht->data[indexAL].data + (index+1), (p_ht->data[indexAL].nbElem-index-1)*sizeof(PtrVoid));
    			 memmove(p_ht->data[indexAL].data + index, p_ht->data[indexAL].data + (index+1), (p_ht->data[indexAL].nbElem-index-1)*sizeof(PtrVoid));
    		 }
    	     p_ht->data[indexAL].nbElem--;

    	     free(p_cell->key);
    	     free(p_cell);
    	     p_ht->nbOps++;
    	}
    }
    return elem;
}


void  HT_printf (FILE* f, PtrHashTable p_ht) {
#define PP fprintf(f,
	PP "HashTable %d/%d (%.2f)\n",p_ht->nbElem, p_ht->size, p_ht->nbElem/(float)p_ht->size);

#undef PP
}

void  HT_rehash (PtrHashTable p_ht, int newSize) {
//WARNING("Not yet implemented\n");
}



PtrHashTableIterator createHashTableIterator(PtrHashTable p_ht) {
	PtrHashTableIterator p_iter = (PtrHashTableIterator) malloc(sizeof(HashTableIterator));
	p_iter->nbOps=p_ht->nbOps;
	p_iter->nbArray=0;
	p_iter->nbElemInA=0;
	p_iter->p_ht = p_ht;

	return p_iter;
}
void                   freeHashTableIterator(PtrHashTableIterator p_iter) {
	free(p_iter);
}
PtrVoid HTI_next(PtrHashTableIterator p_iter) {
/* one should test nbOps before returning ANY VALUE */

	while(p_iter->nbArray<p_iter->p_ht->size) {
    	ArrayList al = p_iter->p_ht->data[p_iter->nbArray];
    	if (p_iter->nbElemInA<al.nbElem) {

    		/* the key and data stored */
    		PtrHashTableCell p_cell = al.data[p_iter->nbElemInA];

    		p_iter->nbElemInA++;
    		/* return the elem */
    		return p_cell->elem;

    	}
		p_iter->nbArray++;
		p_iter->nbElemInA=0;

	}

	return NULL;
}


/*
	p_ht->size;

	ArrayList* newData = (ArrayList*) malloc(newSize*sizeof(ArrayList));

	for(i=0;i<p_ht->size;i++) {
		newData[i].AL_compareFuncEE = data[0].AL_compareFuncEE;
		newData[i].nbElem=0;
		newData[i].nbMaxElem=0;
		newData[i].data=NULL;
	}

	//ICI: copier les elements


	free(p_ht->data);
	p_ht->data = newData;
	p_ht->size = newSize;
	return p_ht;

}
*/

/* test section */

/*
typedef struct {
	int nb;
	char name[50];
} Toto, *PtrToto;

void fprintf_toto(FILE* f, Toto p) {
	fprintf(f, "%s %d\n",p.name, p.nb);
}

PtrVoid freeToto (PtrVoid elem) {
	return NULL;
}



int main() {
	Toto fredo;
	fredo.nb=10;
	strcpy(fredo.name, "Fredo");

	Toto marta;
	marta.nb=256;
	strcpy(marta.name, "Marta");

	Toto choupinou;
	choupinou.nb=25;
	strcpy(choupinou.name, "Choupinou is a motherfucker !");

	Toto aurelie;
	aurelie.nb=30;
	strcpy(aurelie.name, "Aurelie");

	PtrToto p_who_is_that;

	PtrHashTable p_ht = createHashTable(freeToto);

	HT_add(p_ht, fredo.name, &fredo);
	HT_add(p_ht, marta.name, &marta);
	HT_add(p_ht, "Choupinou", &choupinou);
	HT_add(p_ht, "Choupinou", &choupinou);
	HT_add(p_ht, aurelie.name, &aurelie);

	fprintf(stderr, "Get choupinou !\n");
	p_who_is_that = HT_get(p_ht, "Choupinou", 1);
	fprintf_toto(stdout, *p_who_is_that);
	fprintf(stderr, "Get choupinou again !\n");
	p_who_is_that = HT_get(p_ht, "Choupinou", 1);
	fprintf_toto(stdout, *p_who_is_that);

	fprintf(stderr, "Get choupinou again, we love Choupinou so much !\n");
	p_who_is_that = HT_get(p_ht, "Choupinou", 1);
	if (p_who_is_that) {
		fprintf_toto(stdout, *p_who_is_that);
	}
	else {
		fprintf(stderr,"Choupinou disapeared ! Fucking bastard !\n");
	}
	freeHashTable(p_ht);
	return 1;
}

*/

