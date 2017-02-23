#include <stdio.h>
#include <samtools/sam.h>


int len_query(const bam1_t *b) {
  // This one includes soft-clipped nucleotides
  return  bam_cigar2qlen(&b->core, bam1_cigar(b));
}

int len_genome(const bam1_t *b) {
  // length in genome-cooredinates
  return  bam_calend(&b->core, bam1_cigar(b)) - b->core.pos;
}

int len_align(const bam1_t *b)  {
  // pedestrian way, fill in whatever choice you like, currently caculates number of (non-insert, non-delete) 
  // aligned (match & missmatch) nt's of the query
  int i;
  uint32_t *cigar = bam1_cigar(b);
  int l  = 0;
  for(i=0; i < b->core.n_cigar; ++i) {
    const int opl  = cigar[i] >> BAM_CIGAR_SHIFT;
    const int op   = cigar[i] & BAM_CIGAR_MASK;
    if(op == BAM_CMATCH)   l+= opl;
  }
  return l;
}


int main(int argc, char **argv) {
  int (*len_func)(const bam1_t *);
  if(argc < 4) {
    fprintf(stderr, "Insuficient arguments, Usage: %s len mode in.bam out.bam\n", argv[0]);
    exit(-1);
  } 
  switch(*argv[2]) {
  case 'Q': 
    len_func = len_query;
    break;

  case 'A': 
    len_func = len_align;
    break;

  case 'G': 
    len_func = len_genome;
    break;

  default:
    fprintf(stderr, "unsupported mode, known are _Q_uery-length, _A_ligned length and _G_enomic length\n");
    exit(-1);
    break;
  }

  int len=atoi(argv[1]);

  samfile_t *I = samopen(argv[3], "rb", 0);
  if(!I) {
    printf("Can not open sam input %s\n", argv[3]);
    exit(-1);
  }
  samfile_t *O = samopen(argv[4], "wb", I->header);
  if(!O) {
    printf("Can not open sam output %s\n", argv[4]);
    exit(-1);
  }

  bam1_t *b = bam_init1();

  while (samread(I, b) >= 0) {
    int l = len_func(b);
    if(len != len_func(b)) continue;
    samwrite(O, b);
  }
  samclose(I);
  samclose(O);
  bam_destroy1(b);
  exit(0);
}
