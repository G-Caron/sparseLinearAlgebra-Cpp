CC=g++
CFLAGS=-W -O3 -Wall -ansi -pedantic -std=c++11
LIBS=
LDFLAGS=
EXEC=csr

all: $(EXEC)

csr: gnuplot.o mathTools.o CSRmatrix.o CSRmatrixSym.o CSRmain.o
	$(CC) -o csr gnuplot.o mathTools.o CSRmatrix.o CSRmatrixSym.o CSRmain.o
	
gnuplot.o: gnuplot.cc
	$(CC) -o gnuplot.o -c gnuplot.cc $(CFLAGS)
	
mathTools.o: mathTools.cc
	$(CC) -o mathTools.o -c mathTools.cc $(CFLAGS)
	
CSRmatrixSym.o: CSRmatrixSym.cc mathTools.h gnuplot.h	
	$(CC) -o CSRmatrixSym.o -c CSRmatrixSym.cc $(CFLAGS)
	
CSRmatrix.o: CSRmatrix.cc mathTools.h gnuplot.h
	$(CC) -o CSRmatrix.o -c CSRmatrix.cc $(CFLAGS)
	
CSRmain.o: CSRmain.cc CSRmatrix.h mathTools.h
	$(CC) -o CSRmain.o -c CSRmain.cc $(CFLAGS)
	
clean:
	rm -rf *.o
	
mrproper: clean
	rm -rf csr
