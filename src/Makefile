OBJECTS=CRscalapack.o ParallelAgent.o \
		PAdistData.o PAcollectData.o callpdgesv.o callpdgeqrf.o \
		callpdgesvd.o callpdpotrf.o callpdpotri.o callpdsyevd.o \
		callpdgemm.o CRcollectData.o CRdistData.o 
PA_OBJECTS=ParallelAgent.o PAdistData.o PAcollectData.o
CR_OBJECTS=CRscalapack.o callpdgesv.o callpdgeqrf.o callpdgesvd.o callpdgemm.o\
		callpdpotrf.o callpdpotri.o callpdsyevd.o CRcollectData.o \
		CRdistData.o

all: CRDriver RScaLAPACK.so

RScaLAPACK.so: $(PA_OBJECTS)
	$(SHLIB_LINK) -o $@ $(PA_OBJECTS) $(ALL_LIBS) -lg2c $(ALL_CFLAGS) $(ALL_CPPFLAGS)

CRDriver: CRDriver.o $(CR_OBJECTS)
	$(CC) CRDriver.o $(CR_OBJECTS) $(ALL_LIBS) $(ALL_CFLAGS) $(ALL_CPPFLAGS) -lg2c $(PKG_CFLAGS) -o $@
	cp CRDriver ../exec

clean:
	rm -f CRDriver.o $(OBJECTS) CRDriver CRDriver2 RScaLAPACK.so
