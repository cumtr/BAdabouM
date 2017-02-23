#ifndef _H_HashTable 
#define _H_HashTable

#include <stdio.h>
#include "arrayList.h"
#include "utils.h"

#define HASHTABLE_INITSIZE 1000000
#define HASHTABLE_INCRSIZE 2
#define HASHTABLE_LIMIT    0.8


#define HT_FreeElem void (*HT_FreeElemFunc)(PtrVoid elem)



typedef struct s_HashTable {
    HT_FreeElem;
    unsigned int nbOps; /* to be used to prevent a modification while iterating over the elements */
    int nbElem;
    int size;
    ArrayList* data;
} HashTable, *PtrHashTable;

typedef struct s_HashTableIter {
	unsigned int nbOps; /* to be used to prevent a modification while iterating over the elements */
    int nbArray;
    int nbElemInA;
    PtrHashTable p_ht;
} HashTableIterator, *PtrHashTableIterator;


PtrHashTable           createHashTable(HT_FreeElem);
void                   freeHashTable(PtrHashTable p_ht);
void                   emptyingHashTable(PtrHashTable p_ht);

void          HT_add    (PtrHashTable p_ht, char* key, PtrVoid elem);
PtrVoid       HT_get    (PtrHashTable p_ht, char* key, int pop);
void          HT_printf (FILE* f, PtrHashTable p_ht);

void          HT_rehash (PtrHashTable p_ht, int newSize);

PtrHashTableIterator   createHashTableIterator(PtrHashTable p_ht);
void                   freeHashTableIterator(PtrHashTableIterator p_iter);

PtrVoid HTI_next(PtrHashTableIterator p_iter);

#endif
