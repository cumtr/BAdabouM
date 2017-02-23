ifeq ($(strip $(DEBUG)),)
        DEB := -DDEBUG=0
else
        DEB := -DDEBUG=$(DEBUG)
endif

CC           = gcc
CFLAGS       = -g -Wall -O2 $(DEB)
CFLAGS       = -Wall -O3 $(DEB)
INCLUDES     = -I../includes/samtools-0.1.19/ 
INCLUDES_LIB = -L../includes/samtools-0.1.19/ 
CFLAGS   = -g -Wall -O2 $(shell echo $$CFLAGS)
LDFLAGS  = $(shell echo $$LDFLAGS)


BAdabouM: arrayList.o hashTable.o circlePit.o BAdabouM.o
	$(CC) $+ -o $@  $(INCLUDES_LIB) -lbam -lz -lpthread -lm

%.o: %.c utils.h
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

clean:
	\rm *.o BAdabouM
