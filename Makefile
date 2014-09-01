CC = g++
CFLAGS = -Wall -c -I.
COFLAGS = -Wall -O3 -ffast-math -c -I.
PROGRAMS = extract_count_vector dms_single_transcript

.PHONY : all clean

all : $(PROGRAMS)

sam/libbam.a :
	cd sam ; $(MAKE) all

extract_count_vector : extract_count_vector.cpp sam/bam.h sam/sam.h sam/libbam.a
	$(CC) -O3 -ffast-math -I. $< sam/libbam.a -lz -lpthread -o $@

sampling.hpp : boost/random.hpp

DMSTransModel.hpp : sampling.hpp

DMSTransModel.cpp : sampling.hpp DMSTransModel.hpp

DMSTransModel.o : DMSTransModel.cpp boost/random.hpp sampling.hpp DMSTransModel.hpp
	$(CC) $(COFLAGS) $<

dms_single_transcript.o : dms_single_transcript.cpp boost/random.hpp sampling.hpp DMSTransModel.hpp
	$(CC) $(COFLAGS) $<

dms_single_transcript : DMSTransModel.o dms_single_transcript.o 
	$(CC) -O3 -o $@ $^

clean :
	rm -rf $(PROGRAMS) *.o *~
	cd sam ; $(MAKE) clean


