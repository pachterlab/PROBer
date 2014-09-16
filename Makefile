CC = g++
CFLAGS = -Wall -c -I.
COFLAGS = -Wall -O3 -ffast-math -c -I.
PROGRAMS = 

.PHONY : all clean

all : $(PROGRAMS)

sam/libbam.a :
	cd sam ; $(MAKE) all

sampling.hpp : boost/random.hpp

DMSTransModel.hpp : sampling.hpp

DMSTransModel.cpp : sampling.hpp DMSTransModel.hpp

DMSTransModel.o : DMSTransModel.cpp boost/random.hpp sampling.hpp DMSTransModel.hpp
	$(CC) $(COFLAGS) $<

MyHeap.hpp : utils.h

DMSWholeModel.hpp : sam/bam.h sampling.hpp DMSTransModel.hpp

DMSWholeModel.cpp : sam/bam.h utils.h my_assert.h MyHeap.hpp DMSWholeModel.hpp

DMSWholeModel.o : DMSWholeModel.cpp sam/bam.h boost/random.hpp utils.h my_assert.h sampling.hpp MyHeap.hpp DMSTransModel.hpp DMSWholeModel.hpp
	$(CC) $(COFLAGS) $<

dms_single_transcript.o : dms_single_transcript.cpp boost/random.hpp sampling.hpp DMSTransModel.hpp
	$(CC) $(COFLAGS) $<

dms_single_transcript : DMSTransModel.o dms_single_transcript.o 
	$(CC) -O3 -o $@ $^

clean :
	rm -rf $(PROGRAMS) *.o *~
	cd sam ; $(MAKE) clean


