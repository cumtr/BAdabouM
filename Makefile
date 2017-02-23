all:
	cd includes/samtools-0.1.19/; make
	cd src/; make


clean:
	cd includes/samtools-0.1.19/; make clean
	cd src/; make clean
