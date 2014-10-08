LINKLIB= -L/home/c/cs/cs390/local/fftw-2.1.5/install/lib -lsrfftw -lsfftw -lm
INCLUDE= -I/home/c/cs/cs390/local/fftw-2.1.5/install/include/
CFLAGS=-g
CC=mpicc
CFLAGS+= -DPARALLEL

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/c/cs/cs390/local/fftw-2.1.5/install/lib

all: ionz_main


read_param.o: read_param.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(LINKLIB) read_param.c

ionz_mpi.o: ionz_mpi.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(LINKLIB) ionz_mpi.c

allotarrays.o:	allotarrays.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(LINKLIB) allotarrays.c


ionz_main.o: 	ionz_main.c 
	$(CC) -c $(CFLAGS) $(INCLUDE) $(LINKLIB) ionz_main.c 


ionz_funcs.o:	ionz_funcs.c
	$(CC) -c $(CFLAGS) $(INCLUDE)  $(LINKLIB) ionz_funcs.c


ionz_main: ionz_main.o allotarrays.o  ionz_funcs.o  ionz_mpi.o read_param.o
	$(CC) $(CFLAGS) -o ionz_main  ionz_main.o allotarrays.o ionz_funcs.o ionz_mpi.o read_param.o $(LINKLIB) $(INCLUDE)


clean:
	rm -rf *.o
	rm -rf *~









