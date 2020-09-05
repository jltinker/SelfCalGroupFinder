# for my laptop
hd = $(HOME)/cosmo/lib
LIB = -lm -L${hd} -lcutil
CC = gcc

# -- for sirocco
hd = $(HOME)/lib
LIB = -lm -fopenmp -L${hd} -lcutil 
CC = gcc
CFLAGS = -O2 -fopenmp


OBJS2 = kdGroupFinder_omp.o qromo.o midpnt.o polint.o sham.o spline.o splint.o \
	zbrent.o sort2.o kdtree.o fit_clustering_omp.o gasdev.o ran1.o search.o \
	group_center.o fof.o
kdGroupFinder_omp:	$(OBJS2)
	$(CC) -o $@ $(OBJS2) $(LIB)
	cp -f $@ $(HOME)/exec/$@

clean:
	rm -f *.o
