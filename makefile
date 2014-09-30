
LINKLIB=  -L/usr/local/lib/  -lsrfftw -lsfftw -lm
INCLUDE=-I/usr/local/include/
CFLAGS=-g
CC=gcc






allotarrays.o:	allotarrays.c
	$(CC) -c $(CFLAGS) $(INCLUDE) allotarrays.c


ionz_main.o:	ionz_main.c
	$(CC) -c $(CFLAGS) $(INCLUDE) ionz_main.c


ionz_funcs.o:	ionz_funcs.c
	$(CC) -c $(CFLAGS) $(INCLUDE) ionz_funcs.c


ionz_main: ionz_main.o allotarrays.o  ionz_funcs.o  
	$(CC) $(CFLAGS) -o ionz_main  ionz_main.o allotarrays.o ionz_funcs.o  $(LINKLIB)


clean:
	rm -rf *.o
	rm -rf *~










