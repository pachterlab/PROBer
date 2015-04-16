CC = g++
CFLAGS = -Wall -c -I.
COFLAGS = -Wall -O3 -ffast-math -c -I.
PROGRAMS = PROBer-extract-reference-transcripts PROBer-synthesis-reference-transcripts PROBer-preref PROBer-parse-alignments PROBer-run-em PROBer-simulate-reads PROBer-run-em-separate PROBer_single_transcript

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

PROBer-extract-reference-transcripts : Transcript.o Transcripts.o extractRef.o
	$(CC) -O3 -o $@ $^

synthesisRef.o : synthesisRef.cpp utils.h my_assert.h Transcript.hpp Transcripts.hpp
	$(CC) $(COFLAGS) $<

PROBer-synthesis-reference-transcripts : Transcript.o Transcripts.o synthesisRef.o
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

PROBer-preref : RefSeq.o Refs.o preRef.o
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

SamParser.hpp : sam/sam.h sam/bam.h BamAlignment.hpp AlignmentGroup.hpp

SamParser.cpp : sam/sam.h sam/sam_header.h my_assert.h SamParser.hpp

SamParser.o : SamParser.cpp sam/bam.h sam/sam.h my_assert.h CIGARstring.hpp SEQstring.hpp QUALstring.hpp BamAlignment.hpp AlignmentGroup.hpp SamParser.hpp
	$(CC) $(COFLAGS) $<

BamWriter.hpp : sam/bam.h sam/sam.h BamAlignment.hpp AlignmentGroup.hpp

BamWriter.cpp : sam/bam.h sam/sam.h my_assert.h BamWriter.hpp

BamWriter.o : BamWriter.cpp sam/bam.h sam/sam.h my_assert.h CIGARstring.hpp SEQstring.hpp QUALstring.hpp BamAlignment.hpp AlignmentGroup.hpp BamWriter.hpp
	$(CC) $(COFLAGS) $<

MyHeap.hpp : utils.h

parseAlignments.o : parseAlignments.cpp sam/bam.h sam/sam.h utils.h my_assert.h Transcript.hpp Transcripts.hpp CIGARstring.hpp SEQstring.hpp QUALstring.hpp BamAlignment.hpp AlignmentGroup.hpp SamParser.hpp BamWriter.hpp MyHeap.hpp
	$(CC) $(COFLAGS) $<

PROBer-parse-alignments : Transcript.o Transcripts.o SEQstring.o BamAlignment.o SamParser.o BamWriter.o parseAlignments.o sam/libbam.a
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

SequencingModel.cpp : Markov.hpp Profile.hpp QProfile.hpp SequencingModel.hpp

SequencingModel.o : SequencingModel.cpp sam/bam.h boost/random.hpp sampling.hpp utils.h RefSeq.hpp CIGARstring.hpp SEQstring.hpp QUALstring.hpp Markov.hpp Profile.hpp QProfile.hpp SequencingModel.hpp
	$(CC) $(COFLAGS) $<

InMemoryStructs.hpp : utils.h

PROBerReadModel.hpp : utils.h sampling.hpp RefSeq.hpp Refs.hpp SEQstring.hpp QUALstring.hpp CIGARstring.hpp BamAlignment.hpp AlignmentGroup.hpp MateLenDist.hpp SequencingModel.hpp NoiseProfile.hpp QualDist.hpp InMemoryStructs.hpp 

PROBerReadModel.cpp : utils.h sampling.hpp Refs.hpp MateLenDist.hpp SequencingModel.hpp NoiseProfile.hpp QualDist.hpp PROBerReadModel.hpp 

PROBerReadModel.o : PROBerReadModel.cpp sam/bam.h sam/sam.h boost/random.hpp utils.h my_assert.h sampling.hpp RefSeq.hpp Refs.hpp SEQstring.hpp QUALstring.hpp CIGARstring.hpp BamAlignment.hpp AlignmentGroup.hpp MateLenDist.hpp Markov.hpp Profile.hpp QProfile.hpp SequencingModel.hpp NoiseProfile.hpp QualDist.hpp InMemoryStructs.hpp PROBerReadModel.hpp
	$(CC) $(COFLAGS) $<

PROBerTransModel.hpp : utils.h sampling.hpp InMemoryStructs.hpp

PROBerTransModel.cpp : utils.h sampling.hpp PROBerTransModel.hpp

PROBerTransModel.o : PROBerTransModel.cpp boost/random.hpp utils.h sampling.hpp InMemoryStructs.hpp PROBerTransModel.hpp
	$(CC) $(COFLAGS) $<

PROBerWholeModel.hpp : sam/bam.h sampling.hpp Transcripts.hpp InMemoryStructs.hpp PROBerTransModel.hpp

PROBerWholeModel.cpp : sam/bam.h utils.h my_assert.h Transcript.hpp Transcripts.hpp MyHeap.hpp PROBerWholeModel.hpp

PROBerWholeModel.o : PROBerWholeModel.cpp sam/bam.h boost/random.hpp utils.h my_assert.h sampling.hpp Transcript.hpp Transcripts.hpp MyHeap.hpp InMemoryStructs.hpp PROBerTransModel.hpp PROBerWholeModel.hpp
	$(CC) $(COFLAGS) $<

EM.o : EM.cpp sam/bam.h sam/sam.h boost/random.hpp utils.h my_assert.h sampling.hpp RefSeq.hpp Refs.hpp SEQstring.hpp QUALstring.hpp CIGARstring.hpp BamAlignment.hpp AlignmentGroup.hpp MateLenDist.hpp Markov.hpp Profile.hpp QProfile.hpp SequencingModel.hpp NoiseProfile.hpp QualDist.hpp InMemoryStructs.hpp Transcript.hpp Transcripts.hpp MyHeap.hpp SamParser.hpp BamWriter.hpp PROBerTransModel.hpp PROBerWholeModel.hpp PROBerReadModel.hpp
	$(CC) $(COFLAGS) $<

PROBer-run-em : RefSeq.o Refs.o Transcript.o Transcripts.o SEQstring.o BamAlignment.o SamParser.o BamWriter.o MateLenDist.o Markov.o Profile.o QProfile.o SequencingModel.o NoiseProfile.o QualDist.o PROBerTransModel.o PROBerWholeModel.o PROBerReadModel.o EM.o sam/libbam.a
	$(CC) -O3 -o $@ $^ -lz -lpthread

simulation.o : simulation.cpp sam/bam.h sam/sam.h boost/random.hpp utils.h my_assert.h sampling.hpp RefSeq.hpp Refs.hpp SEQstring.hpp QUALstring.hpp CIGARstring.hpp BamAlignment.hpp AlignmentGroup.hpp MateLenDist.hpp Markov.hpp Profile.hpp QProfile.hpp SequencingModel.hpp NoiseProfile.hpp QualDist.hpp InMemoryStructs.hpp Transcript.hpp Transcripts.hpp MyHeap.hpp SamParser.hpp BamWriter.hpp PROBerTransModel.hpp PROBerWholeModel.hpp PROBerReadModel.hpp
	$(CC) $(COFLAGS) $<

PROBer-simulate-reads : RefSeq.o Refs.o Transcript.o Transcripts.o SEQstring.o BamAlignment.o SamParser.o BamWriter.o MateLenDist.o Markov.o Profile.o QProfile.o SequencingModel.o NoiseProfile.o QualDist.o PROBerTransModel.o PROBerWholeModel.o PROBerReadModel.o simulation.o sam/libbam.a
	$(CC) -O3 -o $@ $^ -lz -lpthread




EM_separate.o : EM_separate.cpp sam/bam.h sam/sam.h boost/random.hpp utils.h my_assert.h sampling.hpp RefSeq.hpp Refs.hpp SEQstring.hpp QUALstring.hpp CIGARstring.hpp BamAlignment.hpp AlignmentGroup.hpp MateLenDist.hpp Markov.hpp Profile.hpp QProfile.hpp SequencingModel.hpp NoiseProfile.hpp QualDist.hpp InMemoryStructs.hpp Transcript.hpp Transcripts.hpp MyHeap.hpp SamParser.hpp BamWriter.hpp PROBerTransModel.hpp PROBerWholeModel.hpp PROBerReadModel.hpp
	$(CC) $(COFLAGS) $<

PROBer-run-em-separate : RefSeq.o Refs.o Transcript.o Transcripts.o SEQstring.o BamAlignment.o SamParser.o BamWriter.o MateLenDist.o Markov.o Profile.o QProfile.o SequencingModel.o NoiseProfile.o QualDist.o PROBerTransModel.o PROBerWholeModel.o PROBerReadModel.o EM_separate.o sam/libbam.a
	$(CC) -O3 -o $@ $^ -lz -lpthread




PROBerTransModelS.hpp : utils.h sampling.hpp InMemoryStructs.hpp

PROBerTransModelS.cpp : utils.h sampling.hpp PROBerTransModelS.hpp

PROBerTransModelS.o : PROBerTransModelS.cpp boost/random.hpp utils.h sampling.hpp InMemoryStructs.hpp PROBerTransModelS.hpp
	$(CC) $(COFLAGS) $<

PROBer_single_transcript.o : PROBer_single_transcript.cpp sam/bam.h sam/sam.h utils.h InMemoryStructs.hpp PROBerTransModelS.hpp
	$(CC) $(COFLAGS) $<

PROBer_single_transcript : PROBerTransModelS.o PROBer_single_transcript.o sam/libbam.a
	$(CC) -O3 -o $@ $^ -lz -lpthread





clean :
	rm -rf $(PROGRAMS) *.o *~ *.pyc
	cd sam ; $(MAKE) clean
