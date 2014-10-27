CC = g++
CFLAGS = -Wall -c -I.
COFLAGS = -Wall -O3 -ffast-math -c -I.
PROGRAMS = dms-seq-extract-reference-transcripts dms-seq-synthesis-reference-transcripts dms-seq-preref dms-seq-parse-alignments dms-seq-run-em

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

dms-seq-extract-reference-transcripts : Transcript.o Transcripts.o extractRef.o
	$(CC) -O3 -o $@ $^

synthesisRef.o : synthesisRef.cpp utils.h my_assert.h Transcript.hpp Transcripts.hpp
	$(CC) $(COFLAGS) $<

dms-seq-synthesis-reference-transcripts : Transcript.o Transcripts.o synthesisRef.o
	$(CC) -O3 -o $@ $^

RefSeq.hpp : utils.h 

RefSeq.cpp : RefSeq.hpp

RefSeq.o : RefSeq.cpp utils.h RefSeq.hpp
	$(CC) $(COFLAGS) $<

Refs.hpp : RefSeqPolicy.h PolyARules.h RefSeq.hpp

Refs.cpp : my_assert.h RefSeqPolicy.h PolyARules.h RefSeq.hpp Refs.hpp

Refs.o : Refs.cpp utils.h my_assert.h RefSeqPolicy.h PolyARules.h RefSeq.hpp Refs.hpp
	$(CC) $(COFLAGS) $<

preRef.o : preRef.cpp utils.h my_assert.h PolyARules.h RefSeqPolicy.h AlignerRefSeqPolicy.h RefSeq.hpp Refs.hpp
	$(CC) $(COFLAGS) $<

dms-seq-preref : RefSeq.o Refs.o preRef.o
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

dms-seq-parse-alignments : Transcript.o Transcripts.o SEQstring.o BamAlignment.o SamParser.o BamWriter.o parseAlignments.o sam/libbam.a
	$(CC) -O3 -ffast-math -o $@ $^ -lz -lpthread

sampling.hpp : boost/random.hpp

MateLenDist.hpp : sampling.hpp

MateLenDist.cpp : MateLenDist.hpp

MateLenDist.o : MateLenDist.cpp boost/random.hpp sampling.hpp MateLenDist.hpp
	$(CC) $(COFLAGS) $<

NoiseProfile.hpp : utils.h sampling.hpp SEQstring.hpp

NoiseProfile.cpp : utils.h NoiseProfile.hpp

NoiseProfile.o : NoiseProfile.cpp sam/bam.h boost/random.hpp sampling.hpp utils.h SEQstring.hpp NoiseProfile.hpp
	$(CC) $(COFLAGS) $<

QualDist.hpp : sampling.hpp QUALstring.hpp

QualDist.cpp : utils.h QualDist.hpp

QualDist.o : QualDist.cpp boost/random.hpp sampling.hpp utils.h QUALstring.hpp QualDist.hpp
	$(CC) $(COFLAGS) $<

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

InMemoryStructs.hpp : utils.h

InMemoryStructs.cpp : InMemoryStructs.hpp

InMemoryStructs.o : InMemoryStructs.cpp utils.h InMemoryStructs.hpp
	$(CC) $(COFLAGS) $<

DMSReadModel.hpp : utils.h sampling.hpp RefSeq.hpp Refs.hpp SEQstring.hpp QUALstring.hpp CIGARstring.hpp BamAlignment.hpp AlignmentGroup.hpp MateLenDist.hpp SequencingModel.hpp NoiseProfile.hpp QualDist.hpp InMemoryStructs.hpp 

DMSReadModel.cpp : utils.h sampling.hpp Refs.hpp MateLenDist.hpp SequencingModel.hpp NoiseProfile.hpp QualDist.hpp DMSReadModel.hpp 

DMSReadModel.o : DMSReadModel.cpp sam/bam.h sam/sam.h boost/random.hpp utils.h my_assert.h sampling.hpp RefSeq.hpp Refs.hpp SEQstring.hpp QUALstring.hpp CIGARstring.hpp BamAlignment.hpp AlignmentGroup.hpp MateLenDist.hpp Markov.hpp Profile.hpp QProfile.hpp SequencingModel.hpp NoiseProfile.hpp QualDist.hpp InMemoryStructs.hpp DMSReadModel.hpp
	$(CC) $(COFLAGS) $<

DMSTransModel.hpp : utils.h sampling.hpp InMemoryStructs.hpp

DMSTransModel.cpp : utils.h sampling.hpp DMSTransModel.hpp

DMSTransModel.o : DMSTransModel.cpp boost/random.hpp utils.h sampling.hpp InMemoryStructs.hpp DMSTransModel.hpp
	$(CC) $(COFLAGS) $<

DMSWholeModel.hpp : sam/bam.h sampling.hpp Transcripts.hpp InMemoryStructs.hpp DMSTransModel.hpp

DMSWholeModel.cpp : sam/bam.h utils.h my_assert.h Transcript.hpp Transcripts.hpp MyHeap.hpp DMSWholeModel.hpp

DMSWholeModel.o : DMSWholeModel.cpp sam/bam.h boost/random.hpp utils.h my_assert.h sampling.hpp Transcript.hpp Transcripts.hpp MyHeap.hpp InMemoryStructs.hpp DMSTransModel.hpp DMSWholeModel.hpp
	$(CC) $(COFLAGS) $<

EM.o : EM.cpp sam/bam.h sam/sam.h boost/random.hpp utils.h my_assert.h sampling.hpp RefSeq.hpp Refs.hpp SEQstring.hpp QUALstring.hpp CIGARstring.hpp BamAlignment.hpp AlignmentGroup.hpp MateLenDist.hpp Markov.hpp Profile.hpp QProfile.hpp SequencingModel.hpp NoiseProfile.hpp QualDist.hpp InMemoryStructs.hpp Transcript.hpp Transcripts.hpp MyHeap.hpp SamParser.hpp BamWriter.hpp DMSTransModel.hpp DMSWholeModel.hpp DMSReadModel.hpp
	$(CC) $(COFLAGS) $<

dms-seq-run-em : RefSeq.o Refs.o Transcript.o Transcripts.o SEQstring.o BamAlignment.o SamParser.o BamWriter.o MateLenDist.o Markov.o Profile.o QProfile.o SequencingModel.o NoiseProfile.o QualDist.o InMemoryStructs.o DMSTransModel.o DMSWholeModel.o DMSReadModel.o EM.o sam/libbam.a
	$(CC) -O3 -o $@ $^ -lz -lpthread

clean :
	rm -rf $(PROGRAMS) *.o *~
	cd sam ; $(MAKE) clean
