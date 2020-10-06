override CPPINCLUDES += -DSMIOL_PNETCDF -DSMIOL_AGGREGATION

all: smiol

smiol:
	$(CC) $(CPPINCLUDES) $(CFLAGS) -c smiol_utils.c
	$(CC) $(CPPINCLUDES) $(CFLAGS) -c smiol.c
	$(FC) $(CPPINCLUDES) $(FFLAGS) -c smiolf.F90
	ar cr libsmiol.a smiol.o smiol_utils.o
	ar cr libsmiolf.a smiolf.o

clean:
	$(RM) -f smiol.o smiol_utils.o libsmiol.a
	$(RM) -f smiolf.o smiolf.mod libsmiolf.a
