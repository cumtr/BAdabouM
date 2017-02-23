#include "utils.h"

#include <stdlib.h>
#include "circlePit.h"


#include <stdio.h>


int main() {
  PtrCirclePit p_cp ;
  int i;
  int j;
  int val;
  p_cp = createCirclePit();

  for (i=0;i<1000000;i++) {
    push(p_cp, (void*) i);

    if (i%50000==0) {
      val = (int) pop(p_cp);
      fprintf(stderr, "je sors avec %d\n", val);
    }
    if (i==400000) {
      for (j=0;j<(400000-10);j++)
	val = (int) pop(p_cp);
    }
  }

  circlepit_printf(stderr, p_cp);

  freeCirclePit(p_cp);
  fprintf(stdout, "Salut\n");
  exit(0);
}
