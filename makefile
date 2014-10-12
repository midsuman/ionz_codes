LINKLIB= -L${TACC_FFTW2_LIB} -lfftw -lsrfftw -lsfftw 
INCLUDE= -I${TACC_FFTW2_INC}
CFLAGS=-g
CC=mpicc
CFLAGS+= -DPARALLEL
POSTFLAGS= -lm

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/c/cs/cs390/local/fftw-2.1.5/install/lib

all: ionz_main


read_param.o: read_param.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(LINKLIB) read_param.c $(POSTFLAGS)

ionz_mpi.o: ionz_mpi.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(LINKLIB) ionz_mpi.c $(POSTFLAGS)

ionz_misc.o: ionz_misc.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(LINKLIB) ionz_misc.c $(POSTFLAGS)


ionz_main.o: 	ionz_main.c 
	$(CC) -c $(CFLAGS) $(INCLUDE) $(LINKLIB) ionz_main.c $(POSTFLAGS)


ionz_funcs.o:	ionz_funcs.c
	$(CC) -c $(CFLAGS) $(INCLUDE)  $(LINKLIB) ionz_funcs.c $(POSTFLAGS)


ionz_main: ionz_main.o allotarrays.o  ionz_funcs.o  ionz_mpi.o read_param.o
	$(CC) $(CFLAGS) -o ionz_main  ionz_main.o allotarrays.o ionz_funcs.o ionz_mpi.o read_param.o $(POSTFLAGS)


clean:
	rm -rf *.o
	rm -rf *~









