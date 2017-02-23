#include "utils.h"
#include "circlePit.h"
#include "hashTable.h"

#include "sam.h"

//#include "gsl/gsl_cdf.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include <getopt.h>
#include <ctype.h>


#define REVERSE  0
#define FORWARD  2
#define DANGLING 0
#define OK       1
#define ADD      1
#define SUB      -1

#define MIN(x,y) (((x)<(y))?(x):(y))

int comp_int(const int * a,const int * b) {
    return (*a)-(*b);
}

int comp_bam1_t(const bam1_t* * a, const bam1_t* * b) {
    int diffChr = ((*a)->core.tid)-((*b)->core.tid);
    if (diffChr==0) return ((*a)->core.pos)-((*b)->core.pos);
    return diffChr;
}

static inline void bam_destroy(bam1_t *b) {
    DEBUGP("bam_destroy %p\n", b);
    bam_destroy1(b);
}

static inline bam1_t* bam_init() {
    bam1_t *b = NULL;
    b = bam_init1();
    DEBUGP("bam_init %p\n", b);
    return b;
}

static inline bam1_t* bam_copy(bam1_t *b1,bam1_t *b2) {
    DEBUGP("bam_copy %p %p\n", b1,b2);
    return bam_copy1(b1,b2);
}


void printfDangling(FILE* stream, int matrixCount[]) {
    fprintf(stream, "   %3s %3s\n","F","R");
    fprintf(stream, "OK %3d %3d\nD  %3d %3d\n", matrixCount[FORWARD+OK],
                                                matrixCount[REVERSE+OK],
                                                matrixCount[FORWARD+DANGLING],
                                                matrixCount[REVERSE+DANGLING]);
}

void printfDanglingBis(FILE* stream, int matrixCount1[], int matrixCount2[], int matrixCount3[]) {
    fprintf(stream, "   %3s %3s       %3s %3s       %3s %3s\n","F","R","F","R","F","R");
    fprintf(stream, "OK %3d %3d       %3d %3d       %3d %3d\n", matrixCount1[FORWARD+OK],
            matrixCount1[REVERSE+OK],
            matrixCount2[FORWARD+OK],
            matrixCount2[REVERSE+OK],
            matrixCount3[FORWARD+OK],
            matrixCount3[REVERSE+OK]);
    fprintf(stream, "D  %3d %3d       %3d %3d       %3d %3d\n", matrixCount1[FORWARD+DANGLING],
            matrixCount1[REVERSE+DANGLING],
            matrixCount2[FORWARD+DANGLING],
            matrixCount2[REVERSE+DANGLING],
            matrixCount3[FORWARD+DANGLING],
            matrixCount3[REVERSE+DANGLING]);
}

void printfDanglingf(FILE* stream, float matrixCount[]) {
    fprintf(stream, "   %5s %5s\n","F","R");
    fprintf(stream, "OK %3.2f %3.2f\nD  %3.2f %3.2f\n", matrixCount[FORWARD+OK],
                                                        matrixCount[REVERSE+OK],
                                                        matrixCount[FORWARD+DANGLING],
                                                        matrixCount[REVERSE+DANGLING]);
}


typedef struct _bam_fetch_MQ_data{ //
      
	char* mate_name; 
    int MQ; 
	
} bam_fetch_MQ_data, *Prt_bam_fetch_MQ_data;

Prt_bam_fetch_MQ_data create_bam_fetch_MQ_data(){
	Prt_bam_fetch_MQ_data DATA = (Prt_bam_fetch_MQ_data) malloc(sizeof(bam_fetch_MQ_data));
	
	//DATA->mate_name = NULL;
	DATA->MQ = 0;
	
	return DATA;
}

void destroy_bam_fetch_MQ_data(Prt_bam_fetch_MQ_data p_mq){
	free(p_mq);
}

static int getMQ(bam1_t *b, Prt_bam_fetch_MQ_data userData){
		
	if (strcmp(bam1_qname(b), userData->mate_name) == 0) {
		fprintf (stderr, "MateQuality : %d\n",  b->core.qual);
		userData->MQ = b->core.qual;
	}
	return -1;
}

int isDangling(bam1_t *b, char* fileName, int min_qual){
	
	samfile_t *sf;
	sf = samopen(fileName, "rb", 0);			
	
	bam_index_t *idx;
	idx = bam_index_load(fileName);

	Prt_bam_fetch_MQ_data p_mq_data = NULL;
	p_mq_data = create_bam_fetch_MQ_data();

	p_mq_data->mate_name = bam1_qname(b);

	
	
	bam_fetch(sf->x.bam, idx, b->core.mtid, b->core.mpos, (b->core.mpos + 1), p_mq_data, (bam_fetch_f) getMQ);
	fprintf (stderr, "2 MateQual : %d\n",  p_mq_data->MQ);
	
	bam_index_destroy(idx);
	samclose(sf);
				
	if (p_mq_data->MQ >= min_qual) return 1;
	else return 0;
	
	destroy_bam_fetch_MQ_data(p_mq_data);
}

/**/

/*
void countDangling(bam1_t *b, int matrixCount[], int toAdd) {
    if (b->core.flag&BAM_FMUNMAP) {
        DEBUGP("Dangling %s tid  %d pos  %d : mtid %d mpos %d - %d\n", bam1_qname(b),
               b->core.tid,  b->core.pos,
               b->core.mtid, b->core.mpos, b->core.mpos-b->core.pos);
    }
    matrixCount[((b->core.flag&BAM_FREVERSE)?REVERSE:FORWARD)+ ((b->core.flag&BAM_FPAIRED) && (b->core.flag&BAM_FMUNMAP)?DANGLING:OK)]+=toAdd;
}
*/

void countDangling(bam1_t *b, int matrixCount[], int min_qual, int toAdd) {

	uint32_t ptr = NULL;
	ptr = (uint32_t) bam_aux2i(bam_aux_get(b, "MQ"));

	if (ptr){
		matrixCount[((b->core.flag&BAM_FREVERSE)?REVERSE:FORWARD) + ((b->core.flag&BAM_FPAIRED) && ((b->core.flag&BAM_FMUNMAP) || (ptr < min_qual))?DANGLING:OK)]+=toAdd;
	}
//	if (ptr == NULL){
//		matrixCount[((b->core.flag&BAM_FREVERSE)?REVERSE:FORWARD) + ((b->core.flag&BAM_FPAIRED) && (b->core.flag&BAM_FMUNMAP)?DANGLING:OK)]+=toAdd;
//	}
	if ((! ptr)){
		matrixCount[((b->core.flag&BAM_FREVERSE)?REVERSE:FORWARD) + ((b->core.flag&BAM_FPAIRED) && (b->core.flag&BAM_FMUNMAP)?DANGLING:OK)]+=toAdd;
	}
}

int filterRead(const bam1_t *b, int min_qual) {
    return (b->core.flag&BAM_FUNMAP) || (b->core.flag&BAM_FDUP)|| (b->core.flag&BAM_FSECONDARY) || (b->core.flag&BAM_FQCFAIL) || b->core.qual<min_qual;
}


/*    
    typedef struct { 
        int32_t tid; 
        int32_t pos; 
        uint32_t bin:16, qual:8, l_qname:8; 
        uint32_t flag:16, n_cigar:16; 
        int32_t l_qseq; 
        int32_t mtid; 
        int32_t mpos; 
        int32_t isize; 
    } bam1_core_t;
*/

void FileExistOrIndex(const char* filename, const char* IndexFile) {
  FILE* fptr = fopen(IndexFile, "r");
  if (fptr != NULL) {
    fclose(fptr);
    //fprintf(stderr, "Index file exists\n");
  }
  else {
  	WARNING("Index File doesnt exist,\n\t\t Indexing\n");
  	//fprintf(stderr, "Index File doesnt exist, Indexing\n");
  	bam_index_build(filename);
  }
}

typedef struct _composite_cp {
    PtrCirclePit p_5P;
    PtrCirclePit p_middle;
    PtrCirclePit p_3P;

    int size_5P;
    int size_middle;
    int size_3P;

    int matrixCount_5P[4];
    int matrixCount_middle[4];
    int matrixCount_3P[4];

	// BreakPoint

	int nbSoftCliped_5p;
	int nbSoftCliped_3p;
	int BreakPoint_min;
	int BreakPoint_max;

	// DEL

	int nbReads_DEL_5p;         
	float meanInsert_DEL_5p;   
	float delta_DEL_5p;

	int nbReads_DEL_3p;

	int prev_DEL_pos;
	int prev_DEL_chr;
	
	int Config_DEL;
	int Start_Config_DEL;
	float Start_meanInser_DEL;
	
	int breakPoint_5p_Start_DEL;
	int breakPoint_3p_Start_DEL;
	
	//INV
	
	int prev_INV_pos;
	int prev_INV_chr;	    

	int nbReads_INV_5p;
	float meanDist_INV_5p;
	float varianceDist_INV_5p;
	float delta_INV_5p;
	float m2_INV_5p;
	
	int nbReads_INV_3p;
	float meanDist_INV_3p;
	float varianceDist_INV_3p;
	float delta_INV_3p;
	float m2_INV_3p;
	
	int prev_INV;
	int Config_INV;
	int Start_Config_INV;
	
	// CTX
	
	int nbReads_CTX_5p;
	int chr_start_CTX_5p;
	int	pos_start_CTX_5p;
	int chr_end_CTX_5p;
	int pos_end_CTX_5p;
	
	int nbReads_CTX_3p;
	int chr_start_CTX_3p;
	int pos_start_CTX_3p;
	int chr_end_CTX_3p;
	int pos_end_CTX_3p;
	
	int prev_CTX_pos;
	int prev_CTX_chr;	
	
	// ITX
	
	int nbReads_ITX_5p;
	int chr_start_ITX_5p;
	int	pos_start_ITX_5p;
	int chr_end_ITX_5p;
	int pos_end_ITX_5p;
	
	int nbReads_ITX_3p;
	int chr_start_ITX_3p;
	int pos_start_ITX_3p;
	int chr_end_ITX_3p;
	int pos_end_ITX_3p;
	
	int prev_ITX_pos;
	int prev_ITX_chr;	
	
} CompositeWindow, *PtrCompositeWindow;

PtrCompositeWindow create_compositeWindow(int readSize, int libraryInsertSize) {
    PtrCompositeWindow p_window = (PtrCompositeWindow) malloc(sizeof(CompositeWindow));
    p_window->p_5P     = createCirclePit();
    p_window->p_middle = createCirclePit();
    p_window->p_3P     = createCirclePit();

    p_window->size_5P = libraryInsertSize;
    p_window->size_middle = readSize;
    p_window->size_3P = libraryInsertSize;

    int i;
    for(i=0;i<4;i++){
        p_window->matrixCount_5P[i]=0;
        p_window->matrixCount_middle[i]=0;
        p_window->matrixCount_3P[i]=0;
    }
    
    /* SoftCliped */
    
    p_window->nbSoftCliped_5p = 0;
    p_window->nbSoftCliped_3p = 0;
    p_window->BreakPoint_min = 0;
    p_window->BreakPoint_max = 0;
    
    /* DEL */
	p_window->nbReads_DEL_5p = 0;        
	p_window->meanInsert_DEL_5p = 0.0; 
    p_window->delta_DEL_5p = 0.0;
    
    p_window->nbReads_DEL_3p = 0; 
    
    p_window->prev_DEL_pos = 0;
	p_window->prev_DEL_chr = 0;
    
    p_window->Config_DEL = 0;
	p_window->Start_Config_DEL = 0;
    p_window->Start_meanInser_DEL = 0;
    
    p_window->breakPoint_5p_Start_DEL = 0;
	p_window->breakPoint_3p_Start_DEL = 0;
    
    /* INV */
    p_window->prev_INV_pos = 0;
	p_window->prev_INV_chr = 0;	    
	
    p_window->nbReads_INV_5p = 0;
    p_window->nbReads_INV_3p = 0;
    
    p_window->meanDist_INV_5p = 0.0;
    p_window->varianceDist_INV_5p = 0.0;
    p_window->delta_INV_5p = 0.0;
    p_window->m2_INV_5p = 0.0;
    
    p_window->meanDist_INV_3p = 0.0;
    p_window->varianceDist_INV_3p = 0.0;
    p_window->delta_INV_3p = 0.0;
    p_window->m2_INV_3p = 0.0;
    
    p_window->prev_INV = 0;
    p_window->Config_INV = 0;
    p_window->Start_Config_INV = 0;
    
    /* CTX */
  	p_window->nbReads_CTX_5p = 0;
	p_window->chr_start_CTX_5p = 0;
	p_window->pos_start_CTX_5p = 0;
	p_window->chr_end_CTX_5p = 0;
	p_window->pos_end_CTX_5p = 0;
	
	p_window->nbReads_CTX_3p = 0;
  	p_window->chr_start_CTX_3p = 0;
	p_window->pos_start_CTX_3p = 0;
	p_window->chr_end_CTX_3p = 0;
	p_window->pos_end_CTX_3p = 0; 
  		
	p_window->prev_CTX_pos = 0;
	p_window->prev_CTX_chr = 0;

    /* ITX */
  	p_window->nbReads_ITX_5p = 0;
	p_window->chr_start_ITX_5p = 0;
	p_window->pos_start_ITX_5p = 0;
	p_window->chr_end_ITX_5p = 0;
	p_window->pos_end_ITX_5p = 0;
	
	p_window->nbReads_ITX_3p = 0;
  	p_window->chr_start_ITX_3p = 0;
	p_window->pos_start_ITX_3p = 0;
	p_window->chr_end_ITX_3p = 0;
	p_window->pos_end_ITX_3p = 0; 
  		
	p_window->prev_ITX_pos = 0;
	p_window->prev_ITX_chr = 0;
    
    
    return p_window;
}

void emptyThisWindow(PtrCompositeWindow p_w, void (*destroy) (PtrVoid p_v)) {
    while (!isEmptyCirclePit(p_w->p_5P)) {
        destroy(pop(p_w->p_5P));
    }
    while (!isEmptyCirclePit(p_w->p_middle)) {
        destroy(pop(p_w->p_middle));
    }
    while (!isEmptyCirclePit(p_w->p_3P)) {
        destroy(pop(p_w->p_3P));
    }
    int i;
    for(i=0;i<4;i++){
        p_w->matrixCount_5P[i]=0;
        p_w->matrixCount_middle[i]=0;
        p_w->matrixCount_3P[i]=0;
    }
    
    /* SoftCliped */
    
    p_w->nbSoftCliped_5p = 0;
    p_w->nbSoftCliped_3p = 0;
    p_w->BreakPoint_min = 0;
    p_w->BreakPoint_max = 0;
    
    /* DEL */
    
	p_w->nbReads_DEL_5p = 0;        
	p_w->meanInsert_DEL_5p = 0.0; 
    p_w->delta_DEL_5p = 0.0;
    
    p_w->nbReads_DEL_3p = 0; 
    
    p_w->prev_DEL_pos = 0;
	p_w->prev_DEL_chr = 0;
    
    p_w->Config_DEL = 0;
	p_w->Start_Config_DEL = 0;
    p_w->Start_meanInser_DEL = 0;
    
    p_w->breakPoint_5p_Start_DEL = 0;
	p_w->breakPoint_3p_Start_DEL = 0;
    
    /* INV */
    p_w->prev_INV_pos = 0;
	p_w->prev_INV_chr = 0;	    
	
    p_w->nbReads_INV_5p = 0;
    p_w->nbReads_INV_3p = 0;
    
    p_w->meanDist_INV_5p = 0.0;
    p_w->varianceDist_INV_5p = 0.0;
    p_w->delta_INV_5p = 0.0;
    p_w->m2_INV_5p = 0.0;
    
    p_w->meanDist_INV_3p = 0.0;
    p_w->varianceDist_INV_3p = 0.0;
    p_w->delta_INV_3p = 0.0;
    p_w->m2_INV_3p = 0.0;
    
    p_w->prev_INV = 0;
    p_w->Config_INV = 0;
    p_w->Start_Config_INV = 0;
    
    /* CTX */
  	p_w->nbReads_CTX_5p = 0;
	p_w->chr_start_CTX_5p = 0;
	p_w->pos_start_CTX_5p = 0;
	p_w->chr_end_CTX_5p = 0;
	p_w->pos_end_CTX_5p = 0;
	
	p_w->nbReads_CTX_3p = 0;
  	p_w->chr_start_CTX_3p = 0;
	p_w->pos_start_CTX_3p = 0;
	p_w->chr_end_CTX_3p = 0;
	p_w->pos_end_CTX_3p = 0; 
  		
	p_w->prev_CTX_pos = 0;
	p_w->prev_CTX_chr = 0;

    /* ITX */
  	p_w->nbReads_ITX_5p = 0;
	p_w->chr_start_ITX_5p = 0;
	p_w->pos_start_ITX_5p = 0;
	p_w->chr_end_ITX_5p = 0;
	p_w->pos_end_ITX_5p = 0;
	
	p_w->nbReads_ITX_3p = 0;
  	p_w->chr_start_ITX_3p = 0;
	p_w->pos_start_ITX_3p = 0;
	p_w->chr_end_ITX_3p = 0;
	p_w->pos_end_ITX_3p = 0; 
  		
	p_w->prev_ITX_pos = 0;
	p_w->prev_ITX_chr = 0;

}

void free_compositeWindow(PtrCompositeWindow p_w, void (*destroy) (PtrVoid p_v)) {
    emptyThisWindow(p_w, destroy);
    freeCirclePit(p_w->p_5P);
    freeCirclePit(p_w->p_middle);
    freeCirclePit(p_w->p_3P);
}

int isEmptyWindow(PtrCompositeWindow p_w) {
    return isEmptyCirclePit(p_w->p_5P) && isEmptyCirclePit(p_w->p_middle) && isEmptyCirclePit(p_w->p_3P) ;
}

PtrVoid firstInWindow(PtrCompositeWindow p_w) {
    PtrVoid p_first = NULL;
    //DEBUGP("firstInWindow\n");

    if (! isEmptyCirclePit(p_w->p_5P))
        p_first = first(p_w->p_5P);
    else if (! isEmptyCirclePit(p_w->p_middle))
        p_first = first(p_w->p_middle);
    else if (! isEmptyCirclePit(p_w->p_3P))
        p_first = first(p_w->p_3P);

    //DEBUGP("firstInWindow (out)\n");

    return p_first;
}

PtrVoid lastInWindow(PtrCompositeWindow p_w) {
    PtrVoid p_last = NULL;

    //DEBUGP("lastInWindow\n");
    if (! isEmptyCirclePit(p_w->p_3P))
        p_last = last(p_w->p_3P );
    else if (! isEmptyCirclePit(p_w->p_middle))
        p_last = last(p_w->p_middle );
    else if (! isEmptyCirclePit(p_w->p_5P))
        p_last = last(p_w->p_5P );

    //DEBUGP("lastInWindow (out)\n");

    return p_last;
}

PtrVoid lastIn5P(PtrCompositeWindow p_w) {
    PtrVoid p_last = NULL;

    //DEBUGP("lastInWindow\n");
    if (! isEmptyCirclePit(p_w->p_5P))
        p_last = last(p_w->p_5P );

    //DEBUGP("lastInWindow (out)\n");

    return p_last;
}

inline int sizeWindow(PtrCompositeWindow p_w) {
    return p_w->size_5P + p_w->size_middle + p_w->size_3P;
}

inline int nbElemInWindow(PtrCompositeWindow p_w) {
    return p_w->p_3P->nbElem + p_w->p_middle->nbElem + p_w->p_5P->nbElem;
}

void CountSoftCliped(bam1_t *b, PtrCompositeWindow p_w, int toAdd){
		uint32_t *cigar = bam1_cigar(b);
		int len = (int) bam_cigar2qlen((bam1_core_t*) b, cigar);
	
		int BreakPoint;
		
		int op_5p = cigar[0]&BAM_CIGAR_MASK;
		if (op_5p == BAM_CSOFT_CLIP){
			int l_5p = cigar[0]>>BAM_CIGAR_SHIFT;
			if (l_5p > len/5){ 
				p_w->nbSoftCliped_5p = p_w->nbSoftCliped_5p +toAdd;
				BreakPoint = b->core.pos + 1;	// +1 necessary because alignament start at 0
		 		if (p_w->BreakPoint_min > BreakPoint || p_w->BreakPoint_min < BreakPoint - 100){ 
 					p_w->BreakPoint_min = BreakPoint;
 				}
 				if (p_w->BreakPoint_max < BreakPoint || p_w->BreakPoint_max > BreakPoint + 100){
 					p_w->BreakPoint_max = BreakPoint;
 				}
			}
		}
		
		int nbCigar = b->core.n_cigar-1;
		int op_3p = cigar[nbCigar] & BAM_CIGAR_MASK;
		if (op_3p == BAM_CSOFT_CLIP){
			int l_3p = cigar[nbCigar]>>BAM_CIGAR_SHIFT;
			if (l_3p > len/5){
				p_w->nbSoftCliped_3p = p_w->nbSoftCliped_3p +toAdd;
				int k;
				for (k = 0, BreakPoint = b->core.pos + 1; k < b->core.n_cigar; ++k) {
					int op = cigar[k]&BAM_CIGAR_MASK;
					if(op == BAM_CMATCH || op == BAM_CINS || op == BAM_CHARD_CLIP){
						int Mapped = cigar[k]>>BAM_CIGAR_SHIFT;
						BreakPoint = BreakPoint + Mapped;
					}
				}
			
	 		if (p_w->BreakPoint_min > BreakPoint || p_w->BreakPoint_min < BreakPoint - 100){ 
				p_w->BreakPoint_min = BreakPoint;
			}
 			if (p_w->BreakPoint_max < BreakPoint || p_w->BreakPoint_max > BreakPoint + 100){
 				p_w->BreakPoint_max = BreakPoint;
 			}
		}
	}
}

typedef struct _FenConfirm{ 

        int nbReads_CTX_5p; 
		int nbReads_CTX_3p;
		int pos_CTX_5p;
		int pos_CTX_3p;
		
		int nbReads_ITX_5p; 
		int nbReads_ITX_3p;
		int pos_ITX_5p;
		int pos_ITX_3p;
		
		int nbReads_INV_5p;
		int nbReads_INV_3p;
		int pos_INV_5p;
		int pos_INV_3p;
		
		int nbReads_soft_5p;
		int nbReads_soft_3p;
		int BreakPoint_5p;
		int BreakPoint_3p;
		
    } FenConfirm, *PtrFenConfirm;

PtrFenConfirm create_FenConfirm(){
	PtrFenConfirm Fen = (PtrFenConfirm) malloc(sizeof(FenConfirm));
	
	Fen->nbReads_CTX_5p = 0;
	Fen->nbReads_CTX_3p = 0;
	Fen->pos_CTX_5p = 0;
	Fen->pos_CTX_3p = 0;
	
	Fen->nbReads_ITX_5p = 0;
	Fen->nbReads_ITX_3p = 0;
	Fen->pos_ITX_5p = 0;
	Fen->pos_ITX_3p = 0;
	
	Fen->nbReads_INV_5p = 0;
	Fen->nbReads_INV_3p = 0;
	Fen->pos_INV_5p = 0;
	Fen->pos_INV_3p = 0;	
	
	Fen->nbReads_soft_5p = 0;
	Fen->nbReads_soft_3p = 0;
	Fen->BreakPoint_5p = 0;
	Fen->BreakPoint_3p = 0;
	
	return Fen;
}

void destroy_FenConfirm(PtrFenConfirm Fen){
	free(Fen);
}

static int fetch_func(const bam1_t *b, PtrFenConfirm data) {
	
		/* CTX */
		
	if(! filterRead(b, 60) &&
	(b->core.flag&BAM_FPAIRED) &&   		// paired
	  (! (b->core.flag&BAM_FMUNMAP)) &&		// the mate is NOT unmapped
	  ((b->core.flag&BAM_FREVERSE)) &&		//the read is mapped to the reverse strand 
	  (b->core.mtid != b->core.tid) ){
		data->nbReads_CTX_5p = data->nbReads_CTX_5p + 1;
		if(data->pos_CTX_5p == 0) {data->pos_CTX_5p = b->core.pos;}
	}
	if(! filterRead(b, 60) &&
	  (b->core.flag&BAM_FPAIRED) &&   		// paired
	  (! (b->core.flag&BAM_FMUNMAP)) &&		// the mate is NOT unmapped
	  (! (b->core.flag&BAM_FREVERSE)) &&	//the read is NOT mapped to the reverse strand 
	  (b->core.mtid != b->core.tid)){
		data->nbReads_CTX_3p = data->nbReads_CTX_3p + 1;
		if(data->pos_CTX_3p == 0 || data->pos_CTX_3p < b->core.pos) {data->pos_CTX_3p = b->core.pos;}
	}

		/* ITX */
		
	if(! filterRead(b, 60) &&
	(b->core.flag&BAM_FPAIRED) &&   		// paired
	  (! (b->core.flag&BAM_FMUNMAP)) &&		// the mate is NOT unmapped
	  (! (b->core.flag&BAM_FREVERSE)) &&	//the read is mapped to the reverse strand 
	  (b->core.tid==b->core.mtid) &&
	  ( (b->core.mpos-b->core.pos)>= 2000 ||
	  (b->core.pos-b->core.mpos)>= 2000) ) {
		data->nbReads_ITX_5p = data->nbReads_ITX_5p + 1;
		if(data->pos_ITX_5p == 0){data->pos_ITX_5p = b->core.pos;}
	}
	
	if(! filterRead(b, 60) &&
	  (b->core.flag&BAM_FPAIRED) &&   		// paired
	  (! (b->core.flag&BAM_FMUNMAP)) &&		// the mate is NOT unmapped
	  ((b->core.flag&BAM_FREVERSE)) &&		//the read is NOT mapped to the reverse strand 
	  (b->core.tid==b->core.mtid) &&
	  ((b->core.mpos-b->core.pos)>= 2000 ||
	  ((b->core.pos-b->core.mpos)>= 2000))){
		data->nbReads_ITX_3p = data->nbReads_ITX_3p + 1;
		if(data->pos_ITX_3p == 0 || data->pos_ITX_3p < b->core.pos){data->pos_ITX_3p = b->core.pos;}
	}
	
		/* INV */
	
	if(! filterRead(b, 60) &&
	  (b->core.flag&BAM_FPAIRED) &&   		// paired
	  (! (b->core.flag&BAM_FMUNMAP)) &&		// the mate is NOT unmapped
	  (! (b->core.flag&BAM_FREVERSE)) &&	//the read is NOT mapped to the reverse strand 
	  (! (b->core.flag&BAM_FMREVERSE)) &&	//the mate NOT is mapped to the reverse strand
	  (b->core.tid==b->core.mtid)){			// Reads on the same chr
		data->nbReads_INV_5p = data->nbReads_INV_5p + 1;
		if(data->pos_INV_5p == 0){data->pos_INV_5p = b->core.pos;}
	}
	
	if(! filterRead(b, 60) &&
	  (b->core.flag&BAM_FPAIRED) &&			// paired
	  (! (b->core.flag&BAM_FMUNMAP)) &&		// the mate is NOT unmapped
	  ((b->core.flag&BAM_FREVERSE)) &&		// the read is mapped to the reverse strand
	  ((b->core.flag&BAM_FMREVERSE)) && 	// the mate is mapped to the reverse strand
	  (b->core.tid==b->core.mtid)){			// Reads on the same chr
		data->nbReads_INV_3p = data->nbReads_INV_3p + 1;
		if(data->pos_INV_3p == 0 || data->pos_INV_3p < b->core.pos){data->pos_INV_3p = b->core.pos;}
	}
	
	/* BreakPoint */
			
	uint32_t *cigar = bam1_cigar(b);
	int32_t len = bam_cigar2qlen((bam1_core_t*) b, cigar);

	int BreakPoint;

	int op_5p = cigar[0]&BAM_CIGAR_MASK;
	if (op_5p == BAM_CSOFT_CLIP){
		int l_5p = cigar[0]>>BAM_CIGAR_SHIFT;
		if (l_5p > len/5){ 
			data->nbReads_soft_5p = data->nbReads_soft_5p + 1;
			BreakPoint = b->core.pos+1;	// +1 necessary because alignament start at 0
	 		if (data->BreakPoint_5p > BreakPoint || data->BreakPoint_5p < BreakPoint - 100){ 
				data->BreakPoint_5p = BreakPoint;
			}
			if (data->BreakPoint_3p < BreakPoint || data->BreakPoint_3p > BreakPoint + 100){
				data->BreakPoint_3p = BreakPoint;
			}
		}
	}
	
	int nbCigar = b->core.n_cigar-1;
	int op_3p = cigar[nbCigar] & BAM_CIGAR_MASK;
	if (op_3p == BAM_CSOFT_CLIP){
		int l_3p = cigar[nbCigar]>>BAM_CIGAR_SHIFT;
		if (l_3p > len/5){
			data->nbReads_soft_3p = data->nbReads_soft_3p + 1;
	
			int k; 
			for (k = 0, BreakPoint = b->core.pos+1; k < b->core.n_cigar; ++k) {
				int op = cigar[k]&BAM_CIGAR_MASK;
				if(op == BAM_CMATCH || op == BAM_CINS || op == BAM_CHARD_CLIP){
					int Mapped = cigar[k]>>BAM_CIGAR_SHIFT;
					BreakPoint = BreakPoint + Mapped;
				}
			}
 			if (data->BreakPoint_5p > BreakPoint || data->BreakPoint_5p < BreakPoint - 100){ 
				data->BreakPoint_5p = BreakPoint;
			}
			if (data->BreakPoint_3p < BreakPoint || data->BreakPoint_3p > BreakPoint + 100){
				data->BreakPoint_3p = BreakPoint;
			}
		}
	}
		
	return 0;
}

/* DEL */

void update_5p_DEL(bam1_t *b, PtrCompositeWindow p_w, int maxdist, int toAdd) {
    if ((b->core.flag&BAM_FPAIRED) &&		// paired
       (! (b->core.flag&BAM_FMUNMAP)) &&	// mate is NOT unmapped
       (! (b->core.flag&BAM_FREVERSE)) &&	// Read NOT mapped to the reverse strand
       (b->core.flag&BAM_FMREVERSE) &&		// mate mapped to the reverse strand
       (b->core.tid==b->core.mtid) &&		// Reads on the same chromosome
       ((b->core.mpos-b->core.pos) >= maxdist) &&	// Reads far one from the other
       b->core.mpos > b->core.pos) {		
                
		if((b->core.pos < p_w->prev_DEL_pos + maxdist) && (p_w->prev_DEL_chr == b->core.tid)){ // Same DEL
			p_w->nbReads_DEL_5p = (p_w->nbReads_DEL_5p) + toAdd;
			p_w->delta_DEL_5p = (b->core.mpos-b->core.pos) - p_w->meanInsert_DEL_5p;
			if ((int) (p_w->meanInsert_DEL_5p + (p_w->delta_DEL_5p / p_w->nbReads_DEL_5p)) > 0) {
				p_w->meanInsert_DEL_5p = p_w->meanInsert_DEL_5p + (p_w->delta_DEL_5p / p_w->nbReads_DEL_5p);
			}
		}
		
		else{
			p_w->prev_DEL_pos = b->core.pos;
			p_w->prev_DEL_chr = b->core.tid;
			p_w->nbReads_DEL_5p = 1;
			p_w->meanInsert_DEL_5p = (b->core.mpos-b->core.pos);
		}
	}
}

void update_3p_DEL(bam1_t *b, PtrCompositeWindow p_w, int maxdist, int toAdd) {
    if ((b->core.flag&BAM_FPAIRED) &&		// paired
       (! (b->core.flag&BAM_FMUNMAP)) &&	// mate is NOT unmapped
       (b->core.tid==b->core.mtid) &&		// Reads on the same chromosome
       ((b->core.pos-b->core.mpos)>=maxdist)) {
			p_w->nbReads_DEL_3p = p_w->nbReads_DEL_3p + toAdd;
	}
}


/* INV */

void update_nb5p_INV(bam1_t *b, PtrCompositeWindow p_w, int maxdist, int toAdd){
	if((b->core.flag&BAM_FPAIRED) &&   		// paired
	  (! (b->core.flag&BAM_FMUNMAP)) &&		// the mate is NOT unmapped
	  (! (b->core.flag&BAM_FREVERSE)) &&	//the read is NOT mapped to the reverse strand 
	  (! (b->core.flag&BAM_FMREVERSE)) &&	//the mate NOT is mapped to the reverse strand
	  (b->core.tid==b->core.mtid) &&		//reads on the same chr
	  ((b->core.mpos-b->core.pos)>=1) ){
		if((b->core.pos < p_w->prev_INV_pos + maxdist) && (p_w->prev_INV_chr == b->core.tid)){ // if same INV
		
			p_w->nbReads_INV_5p = (p_w->nbReads_INV_5p) + toAdd;
			p_w->meanDist_INV_5p = b->core.mpos;
			
		}
		else{ // si nouvelle INV
			p_w->prev_INV_pos = b->core.pos;
			p_w->prev_INV_chr = b->core.tid;
			p_w->nbReads_INV_5p = 1;

			p_w->meanDist_INV_5p = b->core.mpos;
			
			p_w->varianceDist_INV_5p = 0.0;
			p_w->m2_INV_5p = 0.0;					
		}
	}
}

void update_nb3p_INV(bam1_t *b, PtrCompositeWindow p_w, int maxdist, int toAdd){
	if((b->core.flag&BAM_FPAIRED) &&		// paired
	  (! (b->core.flag&BAM_FMUNMAP)) &&		// the mate is NOT unmapped
	  ((b->core.flag&BAM_FREVERSE)) &&		// the read is mapped to the reverse strand
	  ((b->core.flag&BAM_FMREVERSE)) && 	// the mate is mapped to the reverse strand
	  (b->core.tid==b->core.mtid) &&		//reads on the same chr
	  ((b->core.mpos-b->core.pos)>=1) ){	//first in pair
		if((b->core.pos < p_w->prev_INV_pos + maxdist) && (p_w->prev_INV_chr == b->core.tid)){ // same INV
		
			p_w->nbReads_INV_3p = (p_w->nbReads_INV_3p) + toAdd;
			p_w->meanDist_INV_3p = b->core.mpos;
			
		}
		else{  // si nouvelle INV
			p_w->prev_INV_pos = b->core.pos;
			p_w->prev_INV_chr = b->core.tid;
			p_w->nbReads_INV_3p = 1;

			p_w->meanDist_INV_3p = b->core.mpos;
			
			p_w->varianceDist_INV_3p = 0.0;
			p_w->m2_INV_3p = 0.0;					
		}
	}
}

/* CTX */

void update_nb5p_CTX(bam1_t *b, PtrCompositeWindow p_w, int maxdist, int toAdd){
	if((b->core.flag&BAM_FPAIRED) &&   		// paired
	  (! (b->core.flag&BAM_FMUNMAP)) &&		// the mate is NOT unmapped
	  (! (b->core.flag&BAM_FREVERSE)) &&
	  (b->core.mtid != b->core.tid)){  	// Reads are not on the same chr

	  	if (b->core.pos < (p_w->pos_start_CTX_5p + maxdist) && b->core.tid == p_w->chr_start_CTX_5p){ // same variant
	  	
	  		if (b->core.mpos < p_w->pos_end_CTX_5p + maxdist && b->core.mtid == p_w->chr_end_CTX_5p){ /* Position du mate semblable aux precedentes */
		 		p_w->nbReads_CTX_5p = p_w->nbReads_CTX_5p + toAdd; 
	  		}
	  	
	  	}	
	  	else { // New variant
	  		p_w->nbReads_CTX_5p = 1;
	  		p_w->chr_start_CTX_5p = b->core.tid;
	  		p_w->pos_start_CTX_5p = b->core.pos;
	  		p_w->chr_end_CTX_5p = b->core.mtid;
	  		p_w->pos_end_CTX_5p = b->core.mpos;
	  		//fprintf(stdout, "J'en ai un !\n");
	  		
		}
	}
}

void update_nb3p_CTX(bam1_t *b, PtrCompositeWindow p_w, int maxdist, int toAdd){
	if((b->core.flag&BAM_FPAIRED) &&   		// paired
	  (! (b->core.flag&BAM_FMUNMAP)) &&		// the mate is NOT unmapped
	  ((b->core.flag&BAM_FREVERSE)) &&		//the read is mapped to the reverse strand 
	  (b->core.mtid != b->core.tid)){  	

	  	if (b->core.pos < (p_w->pos_start_CTX_3p + maxdist) && b->core.tid == p_w->chr_start_CTX_3p){ // same variant
	  	
	  		if (b->core.mpos < (p_w->pos_end_CTX_3p + maxdist) && b->core.mtid==p_w->chr_end_CTX_3p) { /* Position du mate semblable aux precedentes */
		 		p_w->nbReads_CTX_3p = p_w->nbReads_CTX_3p + toAdd;
	  		}
	  		
	  		else { /* Position du mate differente */
	  			p_w->nbReads_CTX_3p = p_w->nbReads_CTX_3p + toAdd;
	  			p_w->pos_end_CTX_3p = b->core.mpos;
	  			p_w->chr_end_CTX_3p = b->core.mtid;
	  		}
	  	}	
	  	else { // New variant
	  		p_w->nbReads_CTX_3p = 1;
	  		p_w->chr_start_CTX_3p = b->core.tid;
	  		p_w->pos_start_CTX_3p = b->core.pos;
	  		p_w->chr_end_CTX_3p = b->core.mtid;
	  		p_w->pos_end_CTX_3p = b->core.mpos;
	  	}
	}
}


/*ITX*/

void update_nb5p_ITX(bam1_t *b, PtrCompositeWindow p_w, int maxdist, int toAdd){
	if((b->core.flag&BAM_FPAIRED) &&   		// paired
	  (! (b->core.flag&BAM_FMUNMAP)) &&		// the mate is NOT unmapped
	  (! (b->core.flag&BAM_FREVERSE)) &&	// the read is NOT mapped to the reverse strand 
	  ((b->core.flag&BAM_FMREVERSE)) &&		// the mate is mapped to the reverse strand
	  (b->core.mtid == b->core.tid) &&      // reads on the same chromosome
	  ((b->core.mpos-b->core.pos) >= (maxdist * 2) || ((b->core.pos-b->core.mpos) >= (maxdist * 2)) )){  	// Reads are suspectly distant

	  	if (b->core.pos < p_w->pos_start_ITX_5p + maxdist && b->core.tid == p_w->chr_start_ITX_5p){ // same variant
	  		if (b->core.mpos < p_w->pos_end_ITX_5p + maxdist && b->core.mtid==p_w->chr_end_ITX_5p){ // Position du mate semblable aux precedentes
		 		p_w->nbReads_ITX_5p = p_w->nbReads_ITX_5p + toAdd;
		 		if (b->core.mpos < p_w->pos_end_ITX_5p){ // Determinaiton point cassure mate
		 			p_w->pos_end_ITX_5p = b->core.mpos;
		 		} 
	  		}
	  		else { // Position du mate differente
	  			p_w->nbReads_ITX_5p = p_w->nbReads_ITX_5p + toAdd;
	  			p_w->pos_end_ITX_5p = b->core.mpos;
	  			p_w->chr_end_ITX_5p = b->core.mtid;
	  		}
	  	}	
	  	else { // New variant
	  		p_w->nbReads_ITX_5p = 1;
	  		p_w->chr_start_ITX_5p = b->core.tid;
	  		p_w->pos_start_ITX_5p = b->core.pos;
	  		p_w->chr_end_ITX_5p = b->core.mtid;
	  		p_w->pos_end_ITX_5p = b->core.mpos;
	  		
		}
	}
}

void update_nb3p_ITX(bam1_t *b, PtrCompositeWindow p_w, int maxdist, int toAdd){
	if((b->core.flag&BAM_FPAIRED) &&   		// paired
	  (! (b->core.flag&BAM_FMUNMAP)) &&		// the mate is NOT unmapped
	  ((b->core.flag&BAM_FREVERSE)) &&		//the read is mapped to the reverse strand 
	  (! (b->core.flag&BAM_FMREVERSE)) &&	//the mate NOT is mapped to the reverse strand
	  (b->core.mtid == b->core.tid) &&      // reads on the same chromosome
	  ((b->core.mpos-b->core.pos)>= (maxdist * 2) || ((b->core.pos-b->core.mpos)>= (maxdist * 2)) )){  	

	  	if (b->core.pos < p_w->pos_start_ITX_3p + maxdist && b->core.tid == p_w->chr_start_ITX_3p){ // same variant
	  		if (b->core.mpos < (p_w->pos_end_ITX_3p + maxdist) && b->core.mtid==p_w->chr_end_ITX_3p) { // Position du mate semblable aux precedentes
		 		p_w->nbReads_ITX_3p = p_w->nbReads_ITX_3p + toAdd;
	  		}
	  		else { // Position du mate differente
	  			p_w->nbReads_ITX_3p = p_w->nbReads_ITX_3p + toAdd;
	  			p_w->pos_end_ITX_3p = b->core.mpos;
	  			p_w->chr_end_ITX_3p = b->core.mtid;
	  		}
	  	}	
	  	else { // New variant
	  		p_w->nbReads_ITX_3p = 1;
	  		p_w->chr_start_ITX_3p = b->core.tid;
	  		p_w->pos_start_ITX_3p = b->core.pos;
	  		p_w->chr_end_ITX_3p = b->core.mtid;
	  		p_w->pos_end_ITX_3p = b->core.mpos;
	  	}
	}
}



void addWindowInit(PtrCompositeWindow p_w, bam1_t * b, int min_qual, int DistValue){
	
    if (b->core.pos <= p_w->size_5P) {
        push(p_w->p_5P, b);
        countDangling(b, p_w->matrixCount_5P, min_qual, ADD);
        update_5p_DEL(b, p_w, DistValue, ADD);
        update_nb5p_INV(b, p_w, DistValue, ADD);
        update_nb5p_CTX(b, p_w, DistValue, ADD);
        update_nb5p_ITX(b, p_w, DistValue, ADD);
    } 
    else if (b->core.pos <= (p_w->size_5P + p_w->size_middle)) {
        push(p_w->p_middle, b);
        countDangling(b, p_w->matrixCount_middle, min_qual, ADD);
        CountSoftCliped(b, p_w, ADD);
    }
    else {
        push(p_w->p_3P, b);
        countDangling(b, p_w->matrixCount_3P, min_qual, ADD);
        update_nb3p_INV(b, p_w, DistValue, ADD);
        update_nb3p_CTX(b, p_w, DistValue, ADD);
        update_nb3p_ITX(b, p_w, DistValue, ADD);
        update_3p_DEL(b, p_w, DistValue, ADD);
    }
}

int shiftWindow(PtrCompositeWindow p_w, bam1_t * b, void (*destroy) (PtrVoid elem), uint32_t * lenChr, int min_qual, int DistValue) {

    DEBUGP("shiftWindow\n");

    int begin_window=-1;
    int pos_read;
    int pos_first;

    int min_decal = 0;
    int decal_read = 0;
    int decal_window[3];

    int added = 0;
    bam1_t * b_inside;

    if (isEmptyWindow(p_w)) {
        if (b == NULL) {
            added = 0;
        }
        if (b != NULL) {
            addWindowInit(p_w, b, min_qual, DistValue);
            added = 1;
        }
    }
    else {
        //DEBUGP("Determining the minimum shift to do\n");

        /* computing the minimum shift to do */
        //DEBUGP("Determing shift for 5P\n");
        if (! isEmptyCirclePit(p_w->p_5P)) {
            pos_first = ((bam1_t *)first(p_w->p_5P ))->core.pos;
            begin_window = pos_first;
            decal_window[0] = pos_first - begin_window;
            if (decal_window[0]<min_decal) min_decal = decal_window[0];
        }
        //DEBUGP("Determing shift for middle\n");
        if (! isEmptyCirclePit(p_w->p_middle)) {
            pos_first = ((bam1_t *)first(p_w->p_middle ))->core.pos;
            if (begin_window==-1) {
                begin_window = pos_first - p_w->size_5P;
            }
            decal_window[1] = pos_first - (begin_window + p_w->size_5P);
            if (decal_window[1]<min_decal) min_decal = decal_window[1];

        }
        //DEBUGP("Determing shift for 3P\n");
        if (! isEmptyCirclePit(p_w->p_3P)) {
            pos_first = ((bam1_t *)first(p_w->p_3P ))->core.pos;
            if (begin_window==-1) {
                begin_window = pos_first - p_w->size_5P - p_w->size_middle;
            }
            decal_window[2] = pos_first - (begin_window + p_w->size_5P + p_w->size_middle);
            if (decal_window[2]<min_decal) min_decal = decal_window[2];
        }
        //DEBUGP("Determing shift for b\n");
        if (b != NULL && ((bam1_t *) firstInWindow(p_w))->core.tid == b->core.tid) {
            pos_read = b->core.pos;
            decal_read = pos_read - (begin_window + p_w->size_5P + p_w->size_middle + p_w->size_3P);
            if (decal_read<min_decal) min_decal = decal_read;
        }
        else {
            //DEBUGP("b is not on the same chromosome\n");
        }
        /*
        DEBUGP("MINIMUM SHIFT ARE (%d) 5P: %d(%s) middle: %d(%s) 3P: %d(%s) b: %d(%s)\n",min_decal, decal_window[0],
                                                                                    ! isEmptyCirclePit(p_w->p_5P)?"OK":"NO",
                                                                                    decal_window[1],
                                                                                    !isEmptyCirclePit(p_w->p_middle)?"OK":"NO",
                                                                                    decal_window[2],
                                                                                    ! isEmptyCirclePit(p_w->p_3P)?"OK":"NO",
                                                                                    decal_read,
                                                                                    b!=NULL?((((bam1_t *) firstInWindow(p_w))->core.tid == b->core.tid)?"OK":"OCHR"):"NO");
         */
        /*
         * we have to empty the window iff :
         * - the read is not on the same chromosome AND
         * - we can not shift the window
         */
         
         /* TODO : check if isEmptyWindow(p_w) is necessary here as we are in the else part of the first test 'if (isEmptyWindow(p_w)) {' 
                   check if b != NULL &&       is necessary as it's after '(b==NULL)'
         */
         
        if (((b==NULL) || 
             isEmptyWindow(p_w) || 
             (b != NULL && ((bam1_t *) firstInWindow(p_w))->core.tid != b->core.tid)) /*&&
                 ( (! isEmptyCirclePit(p_w->p_5P) && decal_window[0]==min_decal &&     (((bam1_t *)first(p_w->p_5P ))->core.pos + p_w->size_5P + p_w->size_middle + p_w->size_3P) > lenChr[((bam1_t *)firstInWindow(p_w))->core.tid]) ||
                   (! isEmptyCirclePit(p_w->p_middle) && decal_window[1]==min_decal && (((bam1_t *)first(p_w->p_middle ))->core.pos + p_w->size_middle + p_w->size_3P)            > lenChr[((bam1_t *)firstInWindow(p_w))->core.tid]) ||
                   (! isEmptyCirclePit(p_w->p_3P) && decal_window[2]==min_decal &&     (((bam1_t *)first(p_w->p_3P ))->core.pos + p_w->size_3P)                                   > lenChr[((bam1_t *)firstInWindow(p_w))->core.tid])
                 )*/
                ) {
            DEBUGP("Case where we are at the end of the chromosome\n");
			
			//fprintf(stderr, "%d\t%d\n", ((bam1_t *) firstInWindow(p_w))->core.tid, b->core.tid);
			
            emptyThisWindow(p_w, destroy);
            if (b!=NULL) {
                addWindowInit(p_w, b, min_qual, DistValue);
                added = 1;
            }
            else {
                added=0;
            }
        }
        /* determining the minimum shift for b */
        else {

            /* general case where
             * we shift and/or add the read
             */

            /* shifting here */

            /*
             * we have to shift the read before the window because the 3P window can become empty
             * after shifting
             */
            int wasEmpty_3P = isEmptyCirclePit(p_w->p_3P);

            //DEBUGP("Testing shifting for read\n");
            if (b!=NULL && ((bam1_t *) firstInWindow(p_w))->core.tid == b->core.tid && decal_read==min_decal) {
                push(p_w->p_3P, (PtrVoid) b);
                countDangling(b, p_w->matrixCount_3P, min_qual, ADD);
                update_nb3p_INV(b ,p_w, DistValue, ADD);
                update_nb3p_CTX(b ,p_w, DistValue, ADD);
                update_nb3p_ITX(b ,p_w, DistValue, ADD);
                update_3p_DEL(b, p_w, DistValue, ADD);
                added= 1;
                //DEBUGP("Inserting read into 3P\n");
            }
            else {
                added= 0;
            }
            /*
             * shifting the window now
             */
            //DEBUGP("Testing shifting for 5P\n");
            
            if ((! isEmptyCirclePit(p_w->p_5P)) && decal_window[0]==min_decal) {
                pos_first = ((bam1_t *)first(p_w->p_5P ))->core.pos;
                while((! isEmptyCirclePit(p_w->p_5P)) && (pos_first==((bam1_t *)first(p_w->p_5P))->core.pos)) {
                    /* remove first there */
                    b_inside = (bam1_t *) pop(p_w->p_5P);
                    countDangling(b_inside, p_w->matrixCount_5P, min_qual, SUB);
                    update_5p_DEL(b_inside, p_w, DistValue, SUB);
                    update_nb5p_INV(b_inside,p_w, DistValue, SUB);
                    update_nb5p_CTX(b_inside,p_w, DistValue, SUB);
                    update_nb5p_ITX(b_inside,p_w, DistValue, SUB);

                    destroy(b_inside);
                    DEBUGP("Destroy one read from 5P\n");
                }
            }
            //DEBUGP("Testing shifting for middle\n");
            if ((! isEmptyCirclePit(p_w->p_middle)) && decal_window[1]==min_decal) {
                pos_first = ((bam1_t *)first(p_w->p_middle ))->core.pos;
                while((! isEmptyCirclePit(p_w->p_middle)) && (pos_first==((bam1_t *)first(p_w->p_middle))->core.pos)) {
                    /* move first there from middle to 5P */
                    b_inside = (bam1_t *) pop(p_w->p_middle);
                    countDangling(b_inside, p_w->matrixCount_5P, min_qual, ADD);
                    CountSoftCliped(b_inside, p_w, SUB);
                    
                    update_5p_DEL(b_inside, p_w, DistValue, ADD);
                    update_nb5p_INV(b_inside, p_w, DistValue, ADD);
                    update_nb5p_CTX(b_inside, p_w, DistValue, ADD);
                    update_nb5p_ITX(b_inside, p_w, DistValue, ADD);


                    countDangling(b_inside, p_w->matrixCount_middle, min_qual, SUB);
                    push(p_w->p_5P, b_inside);
                    DEBUGP("Move one read from middle to 5P\n");
                }
            }
            //DEBUGP("Testing shifting for 3P\n");
            if ((! wasEmpty_3P) && decal_window[2]==min_decal) {
                pos_first = ((bam1_t *)first(p_w->p_3P ))->core.pos;
                while((! isEmptyCirclePit(p_w->p_3P)) && (pos_first==((bam1_t *)first(p_w->p_3P))->core.pos)) {
                    /* move first there from 3P to middle */
                    b_inside = (bam1_t *) pop(p_w->p_3P);
                    countDangling(b_inside, p_w->matrixCount_middle, min_qual, ADD);
                    countDangling(b_inside, p_w->matrixCount_3P, min_qual, SUB);
                    CountSoftCliped(b_inside, p_w, ADD);
                    update_nb3p_INV(b_inside, p_w, DistValue, SUB);
                    update_nb3p_CTX(b_inside, p_w, DistValue, SUB);
                    update_nb3p_ITX(b_inside, p_w, DistValue, SUB);
                    update_3p_DEL(b_inside, p_w, DistValue, SUB);
                    push(p_w->p_middle, b_inside);
                    DEBUGP("Move one read from 3P to middle\n");
                }
            }
        }
    }
    return added;
}

void printHelp() {
#define PP fprintf(stdout,
    PP "\n");
    PP "\tBAdabouM - Compiled on %s\n", __DATE__);
    PP "\n");
    PP  "    Detects putative genomic Structural Varaitions from BAM files\n");
    PP "\n");
    PP "\tGeneral options :\n");
	PP "\n");
    PP "  -h  (Default: Off)\n");
    PP "       Print this help\n");
    PP "  -r  INTEGER (Default: 100)\n");
    PP "       Mean read length\n");
    PP "  -q  INTEGER (Default: 60)\n");
    PP "       min mapping quality threshold\n");
    PP "  -l  INTEGER (Default: Estimated during init phase)\n");
    PP "       Mean library size\n");
    PP "  -I  INTEGER (Default: 100000)\n");
    PP "       Number of read for initiation phase\n");
    //PP "  -p  (Default: Off)\n");
    //PP "       Print reads in window\n");
	PP "\n");
	PP "\n");
	PP "\tSVs reporting options :\n");
	PP "\n");
	PP "  -m   INTEGER (Default: Auto-evaluated) \n");
    PP "       Minimum number of soft-clipped reads for reporting a SV\n");
    PP "  -s  INTEGER (Default: Auto-evaluated)\n");
    PP "       Required number of abnormally mapped read on each side of the breakpoint for report a SV\n");
	PP "\n");
    PP "  -M  INTEGER (Default: mean library size * 1.5)\n");
    PP "       Min DELetion size to report\n");
	PP "\n");
    PP "  -C  INTEGER (Default: mean library size)\n");
    PP "       Min CNV size to report\n");
	PP "\n");
    PP "  -i  INTEGER (Default: mean library size * 1.5)\n");
    PP "       Min INVersion size to report\n");
    PP "\n");
#undef PP
}

int main(int argc, char *argv[]) {
    int status=EXIT_SUCCESS;

    int i;
    samfile_t *fp_in = NULL;
    bam1_t *b = NULL;
    bam_header_t *header = NULL;
    int someData;
    int backgroundMatrixCount[4] = {0,0,0,0};

    PtrCompositeWindow p_window = NULL;

    PtrHashTable p_ht;
    PtrHashTable p_htPrintF;
    PtrHashTable p_htPrintR;

    PtrHashTable p_htPrint;

    i=0;
	
	int Col_Names = 0;

    /*********************/
    /*   Default values  */
    /*********************/
    
    int read_length = 100;
    int library_length = 0;
    float coverage = 0.0;
    int maxPosInit = -1, minPosInit = -1;
    int firstChr = -1;
    int first_pass = 100000;
    int min_qual = 60;
    
    int min_side = 0;
    int min_softCliped = 0;
    int lastInChr;
    int hflag = 0;
    int printreadflag = 0;
    
    int min_size_CNV = 0;    
    int min_mean_DEL = 0;	 
    int min_size_INV = 0;
    
    /*********************/
    /*   Initialisation  */
    /*********************/
    
    /* CNV */
    int is_in_CNV = 0;
    int32_t chr_in_CNV = 0;
    int32_t pos_in_CNV = 0;
    
    int breakpoint_in_CNV_3p = 0;
    int breakpoint_in_CNV_5p = 0;
    
    /* DEL */
    bam1_t *p_last = NULL;
    int32_t prev_chr_DEL = 0;
    
    /* INS */
    int32_t prev_out_INS = 0;
    int32_t prev_chr_INS = 0; 
    
    
    /***************************************/
    /*									   */
    /*		 Argument processing		   */
    /*								  	   */
    /***************************************/
     
    char *wvalue = NULL;
    int c;
    char const * endptr;
    char* fileName;
    opterr = 0;

    while ((c = getopt (argc, argv, "hpr:I:l:s:q:M:C:i:m:")) != -1)
      switch (c)
        {
        case 'h':
            hflag = 1;
            break;
            
        case 'p':
            printreadflag = 1;
            break;
            
         case 'r':
            wvalue = optarg;
            read_length = strtol(wvalue, (char**) &endptr, 10);
            if (read_length==0 || wvalue==endptr) {
                ERROR(3, "Invalid value for Read length parameter\n");
            }
            break;
            
        case 'I':
            wvalue = optarg;
            first_pass = strtol(wvalue, (char**) &endptr, 10);
            if (first_pass==0 || wvalue==endptr) {
                ERROR(3, "Invalid value for Init phase parameter\n");
            }
            break;
            
        case 'l':
            wvalue = optarg;
            library_length = strtol(wvalue, (char**) &endptr, 10);
            if (library_length==0 || wvalue==endptr) {
                WARNING("Invalid value for Library parameter, value will be computed\n");
            }
            break;
            
        case 's':
            wvalue = optarg;
            min_side = strtol(wvalue, (char**) &endptr, 10);
            if (min_side==0 || wvalue==endptr) {
                ERROR(3, "Invalid value for min_side parameter\n");
            }
            break;

        case 'q':
            wvalue = optarg;
            min_qual = strtol(wvalue, (char**) &endptr, 10);
            if (min_qual==0 || wvalue==endptr) {
                ERROR(3, "Invalid value for Minimum quality parameter\n");
            }
            break;

        case 'm':
            wvalue = optarg;
             min_softCliped = strtol(wvalue, (char**) &endptr, 10);
            if ( min_softCliped==0 || wvalue==endptr) {
                ERROR(3, "Invalid value for min_softCliped parameter\n");
            }            
            break; 

            
            /* DEL */
            
        case 'M':
            wvalue = optarg;
             min_mean_DEL = strtol(wvalue, (char**) &endptr, 10);
            if ( min_mean_DEL==0 || wvalue==endptr) {
                ERROR(3, "Invalid value for minimum DELetion size parameter\n");
            }
            break;
            
            /* CNV */
            
        case 'C':
            wvalue = optarg;
             min_size_CNV = strtol(wvalue, (char**) &endptr, 10);
            if ( min_size_CNV==0 || wvalue==endptr) {
                ERROR(3, "Invalid value for minimum CNV size parameter\n");
            }
            break;  
              
            /* INV */
                      
        case 'i':
            wvalue = optarg;
             min_size_INV = strtol(wvalue, (char**) &endptr, 10);
            if ( min_size_INV==0 || wvalue==endptr) {
                ERROR(3, "Invalid value for minimum INVersion size parameter\n");
            }
            break;
               
                        
        case '?':
          if (optopt == 'c')
            fprintf (stderr, "Option -%c requires an argument.\n", optopt);
          else if (isprint (optopt))
            fprintf (stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf (stderr,
                     "Unknown option character `\\x%x'.\n",
                     optopt);
          return 1;
        default:
          abort ();
        }

    if (hflag || optind >= argc){
        printHelp();
        exit(1);
    }

    fileName = argv[optind];
    
    
    /*********************************************/
    /*                                           */
    /*     END of argument processing            */
    /*                                           */
    /*********************************************/
	
    fp_in = samopen(fileName, "rb", 0);

	/******************************/
	/*  Index file if necessary   */
	/******************************/

	char IndexFile[1000];
	strcpy(IndexFile,fileName);
	strcat(IndexFile,".bai");
	
	FileExistOrIndex(fileName, IndexFile);
	
	
	 /****************************************/
	 /*                                      */
	 /*         Initialization phase         */
	 /*                                      */
	 /****************************************/
	
    if(NULL == fp_in) {
      ERROR(3, "Could not open file.\n");
    }

    int* tabL = NULL;
    if (library_length==0) {
        tabL = (int*) malloc(sizeof(int)*first_pass);
        memset(tabL, 0, sizeof(int)*first_pass);
    }

    b = bam_init();                                      /* allocation de memoire pour un alignement */
    while(samread(fp_in, b) > 0 && i<first_pass) {       /* loop over the records */
                        
        if (filterRead(b, min_qual)) { 					 /* unmapped ou pourrave */
            bam_destroy(b);
            b = bam_init();
            continue;
        }
        
        /* coverage */
        if (minPosInit==-1) {
            minPosInit=b->core.pos;
            firstChr= b->core.tid;
        }
         
        if (b->core.tid==firstChr) {
            coverage+=read_length;
            maxPosInit = b->core.pos;
        }
        
        /* dangling */
        countDangling(b, backgroundMatrixCount, min_qual, ADD);

        /* library length */
        if (library_length==0 && !(b->core.flag&BAM_FREVERSE) &&
           !(b->core.flag&BAM_FMUNMAP) && 
           (b->core.tid==b->core.mtid) &&
           (b->core.mpos-b->core.pos) > 0) {
            tabL[i] = (b->core.mpos-b->core.pos);
            i++;
        }
        bam_destroy(b);
        b = bam_init();/* allocation de memoire pour un alignement */
    }
    bam_destroy(b);

    coverage = coverage/(maxPosInit-minPosInit+1);
    WARNING("Coverage estimated to %f\n", coverage);
    
    
    if (library_length==0) {
        qsort(tabL, first_pass, sizeof(int), (const void *) comp_int);
        library_length = tabL[(int) (i/2)]+read_length;
        free(tabL);
        //WARNING("Library length estimated to %d\n", library_length);
    }
    
    if(min_side == 0){
    	min_side = (int) ceil(((coverage / read_length ) * library_length) / 8);
    	WARNING("min_side estimated to %d\n", min_side);
    }
    
    if(min_softCliped == 0){
    	min_softCliped = (int) ceil(((coverage / read_length) * read_length) / 8);
    	//WARNING("min_softCliped estimated to %d\n", min_softCliped);
    }
    
    if (min_mean_DEL==0) {
		min_mean_DEL = library_length * 1.5;  
    }

	if (min_size_CNV==0) {
		min_size_CNV = library_length;
	}
	
	if (min_size_INV==0) {
		min_size_INV = library_length * 1.5;
	}
    //WARNING("Normal number of reads : %f \n",(library_length*coverage/read_length));
    
	 /****************************************/
	 /*                                      */
	 /*         Looping over records         */
	 /*                                      */
	 /****************************************/    

    fp_in = samopen(fileName, "rb", 0);
    header = fp_in->header;

    p_window = create_compositeWindow(read_length, library_length);
    p_ht = createHashTable((void *) bam_destroy);
    p_htPrintF = createHashTable((void *) bam_destroy);
    p_htPrintR = createHashTable((void *) bam_destroy);

    b = bam_init();/* allocation de memoire pour un alignement */

    DEBUGP("BEGIN\n");

    do {                                                          /* loop over the records */ 
        DEBUGP("TRY GETTING A READ\n");
        someData = (samread(fp_in, b)>0);

#ifdef DEBUG
        if (someData) {
            char* s = bam_format1(header, b);
            DEBUGP("%s\n", s);
            free(s);
        }
#endif

        /* ne pas considerer les reads pourris :
            - non mappes 
            - qualite pas bonne
        */
        
        if (someData && filterRead(b, min_qual)) { /* unmapped ou pourrave */
            DEBUGP("I'm unmapped ou tout pourris %s\n", bam1_qname(b));

            /* Sauvegarde du read pour une utilisation future */
            
            if (printreadflag) {
                HT_add(p_ht, bam1_qname(b), (PtrVoid) b);
            }
            else bam_destroy(b);
        }
        
        /**************************************
         * Cas particulier debut du chromosome
         * on ne va pas calculer avec la fenetre
         * sans ce read
         **************************************/
                
        else if (someData &&
                 ( isEmptyWindow(p_window) || ((bam1_t *) lastInWindow(p_window))->core.tid==b->core.tid) &&
                 b->core.pos <= sizeWindow(p_window)
                ) {
           // on ajoute le read, c'est tout
            DEBUGP("BEGINNING OF CHR\n");
            addWindowInit(p_window, b, min_qual, (library_length * 2));
            DEBUGP("added INIT %p %s %d %d\n",b, bam1_qname(b), b->core.tid, b->core.pos);
        } 
        
        /***********************************************
         * Cas particulier meme position que le dernier
         * on ne va pas calculer avec la fenetre
         * sans ce read
         ***********************************************/
        
        else if (  someData &&
                   ( ! isEmptyWindow(p_window)) &&
                   ((bam1_t *) lastInWindow(p_window))->core.tid==b->core.tid &&
                   ((bam1_t *) lastInWindow(p_window))->core.pos==b->core.pos
                ) {
           /*
            * on ajoute le read, c'est tout
            * en plus c'est facile car ca ne change a la fenetre
            * on peut le mettre direct en 3P
            */
            DEBUGP("SAME POSITION AS LAST ONE\n");
            push(p_window->p_3P, b);
            countDangling(b, p_window->matrixCount_3P, min_qual, ADD);
            update_nb3p_INV(b ,p_window, (library_length * 2), ADD);
            update_nb3p_CTX(b ,p_window, (library_length * 2), ADD);
            update_nb3p_ITX(b ,p_window, (library_length * 2), ADD);

            update_3p_DEL(b, p_window, (library_length * 2), ADD);
            DEBUGP("added in 3P %p %s %d %d\n",b, bam1_qname(b), b->core.tid, b->core.pos);
        }
        
        /*****************************************************************************
         * cas ou on ne peut pas l'inserer directement :
         * - il faut calculer avec la fenetre actuelle (CNV, INS, DEL, ...)
         * - shifter
         * - inserer
         *****************************************************************************/
        
        else {
            DEBUGP("WANT TO ADD ONE\n");
            
            if (someData) {
                DEBUGP("WANT TO ADD %s %d %d\n",bam1_qname(b), b->core.tid, b->core.pos);
            }
            else {
                //WARNING("NO MORE DATA\n", 1);
                /*
                 * we just have to empty the window
                 */
            }

            int added=0;
            
            do {
                 
                /***********************/
                /*      Col_Names      */
                /***********************/ 
                 
                 if (Col_Names == 0){
                 	fprintf(stdout, "Chr_Start\tPos_Start_5p\tPos_Start_3p\tChr_End\tPos_End_5p\tPos_End_3p\tType\tSize\n");
                 	Col_Names = 1;
                 }
                 
                /********************************************************************/
                /*                          SVs  DETECTION                          */ 
                /********************************************************************/   
                                  
                                       
                /*****************/
                /*      INS      */
                /*****************/              
                 
                 
                if (!isEmptyWindow(p_window) &&
                
                   ((library_length*coverage/read_length) * 0.5 <= p_window->p_5P->nbElem) && ((library_length*coverage/read_length) * 2 >= p_window->p_5P->nbElem) &&
                   ((library_length*coverage/read_length) * 0.5 <= p_window->p_3P->nbElem) && ((library_length*coverage/read_length) * 2 >= p_window->p_3P->nbElem) &&
                
                    (p_window->matrixCount_5P[FORWARD+DANGLING] >= min_side &&
                    p_window->matrixCount_3P[REVERSE+DANGLING] >= min_side) ) {
                    
					if ((((bam1_t *) first(p_window->p_5P))->core.pos > prev_out_INS ||
					   (((bam1_t *) first(p_window->p_5P))->core.tid) > prev_chr_INS) &&
					   p_window->nbSoftCliped_5p >= min_softCliped && p_window->nbSoftCliped_3p >= min_softCliped) {
						fprintf(stdout,"%s\t%d\t%d\t%s\t%d\t%d\tINS\t%d\n",
								header->target_name[((bam1_t *)lastInWindow(p_window))->core.tid],
								p_window->BreakPoint_min,
								p_window->BreakPoint_max,
								header->target_name[((bam1_t *)firstInWindow(p_window))->core.tid],
								p_window->BreakPoint_min,
								p_window->BreakPoint_max,
								p_window->BreakPoint_max - p_window->BreakPoint_min);
															
						prev_out_INS = ((bam1_t *) lastInWindow(p_window))->core.pos;                        
						prev_chr_INS = ((bam1_t *) lastInWindow(p_window))->core.tid;
					}
                                        
				if ( DEBUG )
                    printfDanglingBis(stderr, p_window->matrixCount_5P, p_window->matrixCount_middle, p_window->matrixCount_3P);

                    /*
                     * use bam1_t *bam_dup1(const bam1_t *src) to make copies of alignment
                     */


                    if (printreadflag) {
                        bam1_t *b1, *b2;
                        DEBUGP("WINDOW 5P\n");

                        for(i=0;i<p_window->p_5P->nbElem;i++){
                            b2 = (bam1_t *) get(p_window->p_5P, i);
                            if (b2->core.flag&BAM_FREVERSE) {
                                p_htPrint=p_htPrintF;
                            }
                            else {
                                p_htPrint=p_htPrintR;
                            }

                            if (HT_get(p_htPrint, bam1_qname(b2), 0)==NULL) {
                                /* b2 n'est deja stocke dans la table de hash */
                                /* il faut en faire une copie */
                                b1 = bam_init();
                                bam_copy(b1, b2);
                                HT_add(p_htPrint, bam1_qname(b1), (PtrVoid) b1);
                            }
                        }
                        DEBUGP("WINDOW 3P\n");

                        for(i=0;i<p_window->p_3P->nbElem;i++){
                            b2 = (bam1_t *) get(p_window->p_3P, i);

                            if (b2->core.flag&BAM_FREVERSE) {
                                p_htPrint=p_htPrintF;
                            }
                            else {
                                p_htPrint=p_htPrintR;
                            }

                            if (HT_get(p_htPrint, bam1_qname(b2), 0)==NULL) {
                                /* b2 n'est deja stocke dans la table de hash */
                                /* il faut en faire une copie */
                                b1 = bam_init();
                                bam_copy(b1, b2);
                                HT_add(p_htPrint, bam1_qname(b1), (PtrVoid) b1);
                            }
                		}
                	}
                }
                        
                      
                /*****************/
                /*      CNV      */
                /*****************/
                                

                if (!isEmptyWindow(p_window) &&
                	( p_window->p_3P->nbElem >= ((library_length * coverage / read_length) * 2.5)) &&
                    ( p_window->p_5P->nbElem >= (library_length*coverage/read_length) * 0.5) &&
                    ( p_window->p_5P->nbElem <= (library_length*coverage/read_length) * 1.5) &&
                    is_in_CNV==0) {
                
					if ((p_window->nbSoftCliped_5p >= (min_softCliped)) ) { //&& (p_window->nbSoftCliped_3p >= (min_softCliped))
												
						is_in_CNV = 1;
						
						chr_in_CNV = ((bam1_t *)lastInWindow(p_window))->core.tid;
						pos_in_CNV = ((bam1_t *)lastInWindow(p_window))->core.pos;
				
						breakpoint_in_CNV_5p = p_window->BreakPoint_min;
						breakpoint_in_CNV_3p = p_window->BreakPoint_max;
					}          
                }

                if (is_in_CNV==1 &&
                   !isEmptyWindow(p_window) &&
                   (p_window->p_5P->nbElem >= ((library_length*coverage/read_length) * 2.5)) &&
                   ((library_length*coverage/read_length) * 1.5 >= p_window->p_3P->nbElem) &&
                   ((library_length*coverage/read_length) * 0.5 <= p_window->p_3P->nbElem) ) {
                    
                    if ((p_window->nbSoftCliped_3p)>=(min_softCliped) &&
                       (p_window->BreakPoint_max - breakpoint_in_CNV_5p) > min_size_CNV &&
                       (p_window->BreakPoint_max - p_window->BreakPoint_min) < (read_length/2)) {
                       
                    	fprintf(stdout, "%s\t%d\t%d\t%s\t%d\t%d\tCNV\t%d\n",
                            header->target_name[((bam1_t *)firstInWindow(p_window))->core.tid],
                           	breakpoint_in_CNV_5p,
                           	breakpoint_in_CNV_3p,
                           	header->target_name[((bam1_t *)firstInWindow(p_window))->core.tid],
 							p_window->BreakPoint_min,
 							p_window->BreakPoint_max,
                           	(p_window->BreakPoint_max - breakpoint_in_CNV_5p));
                    		
                    		is_in_CNV = 0;
                    	
                    }
                 }
                
                if (is_in_CNV==1 &&
                   (p_window->p_5P->nbElem <= (library_length*coverage/read_length)) &&
                   (p_window->p_3P->nbElem <= (library_length*coverage/read_length))) {
                                      
                	is_in_CNV = 0;
                }

                
                /*****************/
                /*       DEL     */
                /*****************/
                
				if (!isEmptyWindow(p_window) &&
					p_window->nbReads_DEL_5p >= min_side &&
					p_window->meanInsert_DEL_5p > min_mean_DEL &&
					((p_window->Start_Config_DEL + min_mean_DEL < ((bam1_t *)firstInWindow(p_window))->core.pos) ||
					(p_window->Start_Config_DEL > ((bam1_t *)firstInWindow(p_window))->core.pos))) { // New DEL
												
						p_window->Config_DEL = 1;
						p_window->Start_Config_DEL = ((bam1_t *)firstInWindow(p_window))->core.pos;
						p_window->Start_meanInser_DEL = p_window->meanInsert_DEL_5p;
						prev_chr_DEL = ((bam1_t *)firstInWindow(p_window))->core.tid;
				
				}
				
				if (p_window->Config_DEL == 1 && p_window->Config_DEL != 2 &&	// DEL ouverte
					p_window->nbSoftCliped_3p >= min_softCliped &&	//Assez de reads softClipped
					//(p_window->Start_Config_DEL + read_length) < b->core.pos && 	//rea
					(p_window->BreakPoint_min + (read_length/2)) > p_window->BreakPoint_max) {	// Same DEL : determining BreakPoint
																
						p_window->breakPoint_5p_Start_DEL = p_window->BreakPoint_min;
						p_window->breakPoint_3p_Start_DEL = p_window->BreakPoint_max;
						p_window->Config_DEL = 2;
								
				}
				
				if (p_window->Config_DEL == 2 &&	// DEL Open
				   ((p_window->breakPoint_5p_Start_DEL + read_length) < b->core.pos) &&	// not the same BreakPoint
					p_window->nbSoftCliped_5p >= min_softCliped &&						 // assez de reads 
					(p_window->BreakPoint_max - p_window->breakPoint_5p_Start_DEL) < (int) (p_window->Start_meanInser_DEL + min_mean_DEL) && 	// Verif que la sortie est coherente
					(p_window->BreakPoint_max - p_window->breakPoint_5p_Start_DEL) > (int) (p_window->Start_meanInser_DEL - min_mean_DEL) ){	// Verif que la sortie est coherente
					   			   
 					   fprintf(stdout, "%s\t%d\t%d\t%s\t%d\t%d\tDEL\t%d\n",
 							header->target_name[((bam1_t *)firstInWindow(p_window))->core.tid],
 							p_window->breakPoint_5p_Start_DEL,
 							p_window->breakPoint_3p_Start_DEL,
 							header->target_name[((bam1_t *)firstInWindow(p_window))->core.tid],
 							p_window->BreakPoint_min,
 							p_window->BreakPoint_max,
 							p_window->BreakPoint_max - p_window->breakPoint_5p_Start_DEL);
 						p_window->Config_DEL = 0;
				}				

				if ( (p_window->Config_DEL != 0) && 
				   (p_window->Start_Config_DEL + (int) p_window->Start_meanInser_DEL + (p_window->Start_meanInser_DEL + (library_length * 2 + read_length))) <= b->core.pos) {
					p_window->Config_DEL = 0;
				}
				
				
                /*****************/
                /*       INV     */
                /*****************/
				
				
				if (!isEmptyWindow(p_window) &&	
				   p_window->nbReads_INV_5p >= min_side &&	// enought reads
				   p_window->nbReads_INV_3p >= min_side &&	// enought reads
				   ((p_window->meanDist_INV_5p - p_window->meanDist_INV_3p) < (library_length * 2)) &&
				   ((p_window->meanDist_INV_5p - p_window->meanDist_INV_3p) > - (library_length * 2))) {
				   
					if (p_window->prev_INV + min_size_INV < ((bam1_t *)firstInWindow(p_window))->core.pos) { // New INV
						p_window->Config_INV = 1;
						p_window->Start_Config_INV = ((bam1_t *)firstInWindow(p_window))->core.pos;
					}
				}
				
				if ( ((p_window->nbSoftCliped_5p + p_window->nbSoftCliped_3p) >= min_softCliped) &&
					(p_window->meanDist_INV_5p - (((bam1_t *)lastInWindow(p_window))->core.pos) > (library_length * 2) ) &&  // Taille min de l'inv !
				    p_window->Config_INV == 1) {
					
					samfile_t *sf;
					sf = samopen(fileName, "rb", 0);
				
					bam_index_t *idx;
					idx = bam_index_load(fileName);
				
					PtrFenConfirm Fenetre = NULL;
					Fenetre = create_FenConfirm();									
					
					bam_fetch(sf->x.bam, idx, ((bam1_t *)firstInWindow(p_window))->core.tid, 
							 ((int) p_window->meanDist_INV_3p) - (read_length * 2), 
							 ((int) p_window->meanDist_INV_3p) + (read_length * 2),
							 Fenetre,
							 (bam_fetch_f) fetch_func);
					
					bam_index_destroy(idx);
					samclose(sf);					
				
					if ((Fenetre->nbReads_INV_5p >= min_side) && (Fenetre->nbReads_INV_3p >= min_side) &&
					((Fenetre->nbReads_soft_5p + Fenetre->nbReads_soft_3p) >= min_softCliped)) {

					   fprintf(stdout, "%s\t%d\t%d\t%s\t%d\t%d\tINV\t%d\n",
							header->target_name[((bam1_t *)firstInWindow(p_window))->core.tid],
							p_window->BreakPoint_min,
							p_window->BreakPoint_max,
							header->target_name[((bam1_t *)firstInWindow(p_window))->core.tid],
							Fenetre->BreakPoint_5p,
							Fenetre->BreakPoint_3p,
							Fenetre->BreakPoint_3p - p_window->BreakPoint_min);
					}
				
					destroy_FenConfirm(Fenetre);	
					
						p_window->prev_INV = ((bam1_t *)firstInWindow(p_window))->core.pos;
						p_window->prev_INV_pos = ((bam1_t *)firstInWindow(p_window))->core.pos;
						p_window->prev_INV_chr = ((bam1_t *)firstInWindow(p_window))->core.tid;
						p_window->Config_INV = 0;
					
				}
				
				
			if (p_window->Start_Config_INV + min_size_INV < b->core.pos) {
				p_window->Config_INV = 0;	
			}
			
			
                /*****************/
                /*       CTX     */
                /*****************/
				
				if(!isEmptyWindow(p_window) &&	
				   p_window->nbReads_CTX_5p >= min_side && // enought reads
				   p_window->nbReads_CTX_3p >= min_side && // enought reads
				   ((p_window->pos_end_CTX_5p - p_window->pos_end_CTX_3p) < (library_length * 2)) && 	// destination windows close one from the other
				   ((p_window->pos_end_CTX_5p - p_window->pos_end_CTX_3p) > 0) &&
				   (p_window->chr_end_CTX_5p == p_window->chr_end_CTX_3p) /*&&			// destimation windows on same chr
				   (p_window->chr_start_CTX_5p != p_window->chr_end_CTX_5p)*/) { 			// destination not on the chr than start
					if (p_window->prev_CTX_pos + (library_length * 2) < ((bam1_t *)firstInWindow(p_window))->core.pos ||
						p_window->prev_CTX_chr != ((bam1_t *)firstInWindow(p_window))->core.tid) {
						
						if (((p_window->nbSoftCliped_5p + p_window->nbSoftCliped_3p) >= min_softCliped * 1.5) ||
						   (p_window->nbSoftCliped_5p  >= min_softCliped && p_window->nbSoftCliped_3p >= min_softCliped)) {
						   						   
							samfile_t *sf;
							sf = samopen(fileName, "rb", 0);
							
							bam_index_t *idx;
							idx = bam_index_load(fileName);
							
							PtrFenConfirm Fenetre = NULL;
							Fenetre = create_FenConfirm();
							
							bam_fetch(sf->x.bam, idx, p_window->chr_end_CTX_3p, 
									  p_window->pos_end_CTX_3p - (read_length * 2),
									  p_window->pos_end_CTX_5p + (read_length * 2) ,
									  Fenetre,
									  (bam_fetch_f) fetch_func);
						
							bam_index_destroy(idx);
							samclose(sf);
							
							if (Fenetre->nbReads_CTX_5p >= min_side && Fenetre->nbReads_CTX_3p >= min_side &&
							   (Fenetre->nbReads_soft_5p >= min_softCliped || Fenetre->nbReads_soft_3p >= min_softCliped) /*&&
							   (Fenetre->pos_CTX_5p) <= (Fenetre->BreakPoint_5p) && (Fenetre->BreakPoint_5p) <= (Fenetre->pos_CTX_3p)*/){
								
								if (p_last!= NULL) bam_destroy(p_last);
								p_last = bam_init();
								bam_copy(p_last, b);
								
								fprintf(stdout, "%s\t%d\t%d\t%s\t%d\t%d\tCTX\t0\n",
									header->target_name[p_last->core.tid],
									p_window->BreakPoint_min,
									p_window->BreakPoint_max,
									header->target_name[p_window->chr_end_CTX_5p],
									Fenetre->BreakPoint_5p,
									Fenetre->BreakPoint_3p);
							
							}
						
						destroy_FenConfirm(Fenetre);
						
						p_window->prev_CTX_pos = ((bam1_t *)firstInWindow(p_window))->core.pos;
						p_window->prev_CTX_chr = ((bam1_t *)firstInWindow(p_window))->core.tid;	
						
						}
					}
				}
				
				
                /*****************/
                /*       ITX     */
                /*****************/
				
				if(!isEmptyWindow(p_window) &&	
				   p_window->nbReads_ITX_5p >= min_side && // enought reads
				   p_window->nbReads_ITX_3p >= min_side && // enought reads
				   ((p_window->pos_end_ITX_5p - p_window->pos_end_ITX_3p) < (library_length * 2)) && 	// destination windows close one from the other
				   ((p_window->pos_end_ITX_5p - p_window->pos_end_ITX_3p) > 0) &&
				   (p_window->chr_end_ITX_5p == p_window->chr_end_ITX_3p)) { 			// destination not on the chr than start
				   
					if (p_window->prev_ITX_pos + (library_length * 2) < ((bam1_t *)firstInWindow(p_window))->core.pos ||
						p_window->prev_ITX_chr != ((bam1_t *)firstInWindow(p_window))->core.tid) {
						
						if (((p_window->nbSoftCliped_5p + p_window->nbSoftCliped_3p) >= min_softCliped * 1.5) ||
						   (p_window->nbSoftCliped_5p  >= min_softCliped && p_window->nbSoftCliped_3p >= min_softCliped)) {
						   		   
							samfile_t *sf;
							sf = samopen(fileName, "rb", 0);
							
							bam_index_t *idx;
							idx = bam_index_load(fileName);
							
							PtrFenConfirm Fenetre = NULL;
							Fenetre = create_FenConfirm();
							
							bam_fetch(sf->x.bam, idx, p_window->chr_end_ITX_3p, 
									  p_window->pos_end_ITX_3p - library_length * 2, 
									  p_window->pos_end_ITX_5p + library_length * 2,
									  Fenetre,
									  (bam_fetch_f) fetch_func);
						
							bam_index_destroy(idx);
							samclose(sf);
							
 							fprintf(stderr, "ITX fetched ! chr : %d  3p : %d 5p : %d - %d %d %d %d \n",
 							                p_window->chr_end_ITX_3p,
 							                p_window->pos_end_ITX_5p,
 							                p_window->pos_end_ITX_3p,
 							                Fenetre->nbReads_ITX_5p,
 							                Fenetre->nbReads_ITX_3p,
 							                Fenetre->nbReads_soft_5p,
 							                Fenetre->nbReads_soft_3p);
 							
							
							if (Fenetre->nbReads_ITX_5p >= min_side && Fenetre->nbReads_ITX_3p >= min_side &&
							   (Fenetre->nbReads_soft_5p >= min_softCliped || Fenetre->nbReads_soft_3p >= min_softCliped) ){
								
								if (p_last!= NULL) bam_destroy(p_last);
								p_last = bam_init();
								bam_copy(p_last, b);
								
								fprintf(stdout, "%s\t%d\t%d\t%s\t%d\t%d\tITX\t0\n",
									header->target_name[p_last->core.tid],
									p_window->BreakPoint_min,
									p_window->BreakPoint_max,
									header->target_name[p_window->chr_end_ITX_5p],
									Fenetre->BreakPoint_5p,
									Fenetre->BreakPoint_3p);
							
							}
						
						destroy_FenConfirm(Fenetre);
						
						p_window->prev_ITX_pos = ((bam1_t *)firstInWindow(p_window))->core.pos;
						p_window->prev_ITX_chr = ((bam1_t *)firstInWindow(p_window))->core.tid;	
						
						}
					}
				}				
				
				
				
				

                if (someData) {
                    DEBUGP("%s %d %d WAITING\n",bam1_qname(b), b->core.tid, b->core.pos);
                }

                /* shift d'une position en fonction de b */
                lastInChr = someData && (!isEmptyWindow(p_window)) && b->core.tid!=((bam1_t*)(lastInWindow(p_window)))->core.tid;

                added = shiftWindow(p_window, someData?b:NULL, (void *) bam_destroy, header->target_len, min_qual, (library_length * 2));
                if (added) {
                    DEBUGP("JUST ADDED %s %d %d\n",bam1_qname(b), b->core.tid, b->core.pos);
                }

                DEBUGP("%d %d %d\n", (someData!= 0), (added==1), !isEmptyWindow(p_window));             
                
            DEBUGP("Fin du while\n");

            
            } while( (someData && ! added) || (! someData && !isEmptyWindow(p_window)));


            if (printreadflag && (lastInChr || isEmptyWindow(p_window))){
                //print
                /* pour tout les elements b2 de p_htPrint */
                bam1_t** p_tab = (bam1_t**) malloc((p_htPrintF->nbElem+p_htPrintR->nbElem)*sizeof(bam1_t*));

                bam1_t *b1, *b2;


                PtrHashTableIterator p_iter = createHashTableIterator(p_htPrintF);
                for(i=0;i<p_htPrintF->nbElem;i++) {
                    p_tab[i]=HTI_next(p_iter);
                }
                freeHashTableIterator(p_iter);
                p_iter = createHashTableIterator(p_htPrintR);
                for(i=0;i<p_htPrintR->nbElem;i++) {
                    p_tab[p_htPrintF->nbElem+i]=HTI_next(p_iter);
                }
                freeHashTableIterator(p_iter);

                qsort(p_tab, p_htPrintF->nbElem+p_htPrintR->nbElem, sizeof(bam1_t*), (const void *) comp_bam1_t);

                for(i=0;i<(p_htPrintF->nbElem+p_htPrintR->nbElem);i++) {
                    b2 = p_tab[i];
                    if (b2->core.flag&BAM_FMUNMAP){
                        b1 = (bam1_t *) HT_get(p_ht, bam1_qname(b2), 0);
                        if (b1!=NULL) {
                            DEBUGP("GOT %s %d %d AS MATE OF %s %d %d FROM HT\n",bam1_qname(b1), b1->core.tid, b1->core.pos,
                                                                                bam1_qname(b2), b2->core.tid, b2->core.pos);

                            char* sdangling = bam_format1(header, b2);
                            char* smate = bam_format1(header, b1);

                            fprintf(stdout, "ALIGNED %s%s\n", b2->core.flag&BAM_FREVERSE?"reverse ":"direct ", sdangling);
                            fprintf(stdout, "DANGLING %s%s\n", b2->core.flag&BAM_FREVERSE?"direct ":"reverse ", smate);

                            free(sdangling);
                            free(smate);

                        }
                        else{
                            ERROR("UNABLE TO GET MATE OF %s %d %d FROM HT\n",bam1_qname(b2), (int) b2->core.tid, b2->core.pos);
                        }
                    }
                    else {
                        char* sdangling = bam_format1(header, b2);
                        fprintf(stdout, "ALIGNEDPAIR %s%s\n", b2->core.flag&BAM_FREVERSE?"reverse ":"direct ", sdangling);
                        free(sdangling);
                    }
                }
                //free
                free(p_tab);
                emptyingHashTable(p_htPrintF);
                emptyingHashTable(p_htPrintR);
                emptyingHashTable(p_ht);
            }

            if (isEmptyWindow(p_window)) {
                /*
                 * stop iterating over records as none has been added
                 */
                break;
            }
        }

        b = bam_init();
    } while(1);
    
    bam_destroy(b);
    //WARNING("free_compositeWindow(p_window, bam_destroy);\n");
    free_compositeWindow(p_window, (void *) bam_destroy);
    //WARNING("freeHashTable(p_ht);\n");
    freeHashTable(p_ht);
    freeHashTable(p_htPrintF);
    freeHashTable(p_htPrintR);
    //WARNING("bam_header_destroy(header);\n");
    bam_header_destroy(header);

    return status;
}
