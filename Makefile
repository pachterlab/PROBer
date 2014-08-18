CC = g++
CFLAGS = -Wall -c -I.
COFLAGS = -Wall -O3 -ffast-math -c -I.
PROGRAMS = rsem-extract-reference-transcripts rsem-synthesis-reference-transcripts rsem-preref rsem-parse-alignments putTogether rsem-tbam2gbam

.PHONY : all clean

all : $(PROGRAMS)

sam/libbam.a :
	cd sam ; $(MAKE) all

Transcript.cpp : utils.h my_assert.h Transcript.hpp

Transcript.o : Transcript.cpp utils.h my_assert.h Transcript.hpp
	$(CC) $(COFLAGS) $<

Transcripts.hpp : Transcript.hpp

Transcripts.cpp : utils.h my_assert.h Transcript.hpp Transcripts.hpp

Transcripts.o : Transcripts.cpp utils.h my_assert.h Transcript.hpp Transcripts.hpp 
	$(CC) $(COFLAGS) $<

extractRef.o : extractRef.cpp utils.h my_assert.h GTFItem.h Transcript.hpp Transcripts.hpp
	$(CC) $(COFLAGS) $<

rsem-extract-reference-transcripts : Transcript.o Transcripts.o extractRef.o
	$(CC) -O3 -o $@ $^

synthesisRef.o : synthesisRef.cpp utils.h my_assert.h Transcript.hpp Transcripts.hpp
	$(CC) $(COFLAGS) $<

rsem-synthesis-reference-transcripts : Transcript.o Transcripts.o synthesisRef.o
	$(CC) -O3 -o $@ $^

RefSeq.cpp : RefSeq.hpp

RefSeq.o : RefSeq.cpp RefSeq.hpp
	$(CC) $(COFLAGS) $<

Refs.hpp : RefSeqPolicy.h PolyARules.h RefSeq.hpp

Refs.cpp : utils.h my_assert.h RefSeqPolicy.h PolyARules.h RefSeq.hpp Refs.hpp

Refs.o : Refs.cpp utils.h my_assert.h RefSeqPolicy.h PolyARules.h RefSeq.hpp Refs.hpp
	$(CC) $(COFLAGS) $<

preRef.o : preRef.cpp utils.h my_assert.h PolyARules.h RefSeqPolicy.h AlignerRefSeqPolicy.h RefSeq.hpp Refs.hpp
	$(CC) $(COFLAGS) $<

rsem-preref : RefSeq.o Refs.o preRef.o
	$(CC) -O3 -o $@ $^

CIGARstring.hpp : sam/bam.h

SEQstring.hpp : sam/bam.h

SEQstring.cpp : SEQstring.hpp

SEQstring.o : SEQstring.cpp sam/bam.h SEQstring.hpp
	$(CC) $(COFLAGS) $<

BamAlignment.hpp : sam/bam.h sam/sam.h my_assert.h CIGARstring.hpp SEQstring.hpp QUALstring.hpp

BamAlignment.cpp : sam/bam.h sam/sam.h my_assert.h SEQstring.hpp BamAlignment.hpp

BamAlignment.o : BamAlignment.cpp sam/bam.h sam/sam.h my_assert.h CIGARstring.hpp SEQstring.hpp QUALstring.hpp BamAlignment.hpp
	$(CC) $(COFLAGS) $<

AlignmentGroup.hpp : sam/sam.h SEQstring.hpp QUALstring.hpp BamAlignment.hpp

SamParser.hpp : sam/sam.h BamAlignment.hpp AlignmentGroup.hpp

SamParser.cpp : sam/sam.h my_assert.h SamParser.hpp

SamParser.o : SamParser.cpp sam/bam.h sam/sam.h my_assert.h CIGARstring.hpp SEQstring.hpp QUALstring.hpp BamAlignment.hpp AlignmentGroup.hpp SamParser.hpp
	$(CC) $(COFLAGS) $<

BamWriter.hpp : sam/bam.h sam/sam.h BamAlignment.hpp AlignmentGroup.hpp 

BamWriter.cpp : sam/bam.h sam/sam.h my_assert.h BamWriter.hpp

BamWriter.o : BamWriter.cpp sam/bam.h sam/sam.h my_assert.h CIGARstring.hpp SEQstring.hpp QUALstring.hpp BamAlignment.hpp AlignmentGroup.hpp BamWriter.hpp
	$(CC) $(COFLAGS) $<

MyHeap.hpp : utils.h

parseAlignments.o : parseAlignments.cpp sam/bam.h sam/sam.h utils.h my_assert.h Transcript.hpp Transcripts.hpp CIGARstring.hpp SEQstring.hpp QUALstring.hpp BamAlignment.hpp AlignmentGroup.hpp SamParser.hpp BamWriter.hpp MyHeap.hpp
	$(CC) $(COFLAGS) $<

rsem-parse-alignments : Transcript.o Transcripts.o SEQstring.o BamAlignment.o SamParser.o BamWriter.o parseAlignments.o sam/libbam.a 
	$(CC) -O3 -ffast-math -o $@ $^ -lz -lpthread  

putTogether.o : putTogether.cpp sam/bam.h sam/sam.h utils.h CIGARstring.hpp SEQstring.hpp QUALstring.hpp BamAlignment.hpp AlignmentGroup.hpp SamParser.hpp BamWriter.hpp
	$(CC) $(COFLAGS) $<

putTogether : SEQstring.o BamAlignment.o SamParser.o BamWriter.o putTogether.o sam/libbam.a
	$(CC) -O3 -ffast-math -o $@ $^ -lz -lpthread

TransBamAlignment.hpp : sam/bam.h Transcript.hpp BamAlignment.hpp

TransBamAlignment.cpp : sam/bam.h utils.h my_assert.h Transcript.hpp CIGARstring.hpp BamAlignment.hpp TransBamAlignment.hpp

TransBamAlignment.o : TransBamAlignment.cpp sam/bam.h sam/sam.h utils.h my_assert.h Transcript.hpp CIGARstring.hpp SEQstring.hpp QUALstring.hpp BamAlignment.hpp TransBamAlignment.hpp
	$(CC) $(COFLAGS) $<

TransAlignmentGroup.hpp : sam/sam.h Transcript.hpp Transcripts.hpp TransBamAlignment.hpp

BamConverter.hpp : sam/bam.h sam/sam.h Transcripts.hpp TransAlignmentGroup.hpp

BamConverter.cpp : sam/bam.h sam/sam.h utils.h my_assert.h Transcripts.hpp TransAlignmentGroup.hpp BamConverter.hpp 

BamConverter.o : BamConverter.cpp sam/bam.h sam/sam.h utils.h my_assert.h Transcript.hpp Transcripts.hpp CIGARstring.hpp SEQstring.hpp QUALstring.hpp BamAlignment.hpp TransBamAlignment.hpp TransAlignmentGroup.hpp BamConverter.hpp
	$(CC) $(COFLAGS) $<

tbam2gbam.o : tbam2gbam.cpp sam/bam.h sam/sam.h utils.h my_assert.h Transcript.hpp Transcripts.hpp CIGARstring.hpp SEQstring.hpp QUALstring.hpp BamAlignment.hpp TransBamAlignment.hpp TransAlignmentGroup.hpp BamConverter.hpp
	$(CC) $(COFLAGS) $<

rsem-tbam2gbam : Transcript.o Transcripts.o SEQstring.o BamAlignment.o TransBamAlignment.o BamConverter.o tbam2gbam.o sam/libbam.a
	$(CC) -O3 -ffast-math -o $@ $^ -lz -lpthread

wiggle.cpp : sam/bam.h sam/sam.h utils.h wiggle.h

wiggle.o : wiggle.cpp sam/bam.h sam/sam.h wiggle.h
	$(CC) $(COFLAGS) $<

sampling.hpp : boost/random.hpp

Orientation.hpp : sampling.hpp

Markov.hpp : utils.h sampling.hpp

Markov.cpp : utils.h Markov.hpp

Markov.o : Markov.cpp boost/random.hpp sampling.hpp utils.h Markov.hpp
	$(CC) $(COFLAGS) $<

QProfile.hpp : utils.h sampling.hpp

QProfile.cpp : utils.h QProfile.hpp

QProfile.o : QProfile.cpp boost/random.hpp sampling.hpp utils.h QProfile.hpp
	$(CC) $(COFLAGS) $<

Profile.hpp : utils.h sampling.hpp

Profile.cpp : utils.h Profile.hpp

Profile.o : Profile.cpp boost/random.hpp sampling.hpp utils.h Profile.hpp
	$(CC) $(COFLAGS) $<

SequencingModel.hpp : RefSeq.hpp CIGARstring.hpp SEQstring.hpp QUALstring.hpp Markov.hpp Profile.hpp QProfile.hpp

SequencingModel.cpp : Markov.hpp Profile.hpp Qprofile.hpp SequencingModel.hpp

SequencingModel.o : SequencingModel.cpp sam/bam.h boost/random.hpp sampling.hpp utils.h RefSeq.hpp CIGARstring.hpp SEQstring.hpp QUALstring.hpp Markov.hpp Profile.hpp QProfile.hpp SequencingModel.hpp
	$(CC) $(COFLAGS) $<

NoiseProfile.hpp : sampling.hpp SEQstring.hpp

NoiseProfile.cpp : utils.h NoiseProfile.hpp

NoiseProfile.o : NoiseProfile.cpp sam/bam.h boost/random.hpp sampling.hpp utils.h SEQstring.hpp NoiseProfile.hpp
	$(CC) $(COFLAGS) $<

QualDist.hpp : sampling.hpp QUALstring.hpp

QualDist.cpp : utils.h QualDist.hpp

QualDist.o : QualDist.cpp boost/random.hpp sampling.hpp utils.h QualDist.hpp
	$(CC) $(COFLAGS) $<

LenDist.hpp : utils.h sampling.hpp

LenDist.cpp : boost/math/distributions/normal.hpp utils.h my_assert.h LenDist.hpp

LenDist.o : LenDist.cpp boost/random.hpp boost/math/distributions/normal.hpp utils.h my_assert.h sampling.hpp LenDist.hpp
	$(CC) $(COFLAGS) $<

RSPD.hpp : sampling.hpp

RSPD.cpp : RSPD.hpp

RSPD.o : RSPD.cpp boost/random.hpp sampling.hpp RSPD.hpp
	$(CC) $(COFLAGS) $<

GroupInfo.cpp : my_assert.h GroupInfo.hpp

GroupInfo.o : GroupInfo.cpp my_assert.h GroupInfo.hpp
	$(CC) $(COFLAGS) $<

clean :
	rm -rf $(PROGRAMS) *.o *~
	cd sam ; $(MAKE) clean


