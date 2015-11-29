include_directories(${ZLIB_INCLUDE_DIRS})
include_directories(../ext)

file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)
list(REMOVE_ITEM sources extractRef.cpp synthesisRef.cpp preRef.cpp parseAlignments.cpp EM.cpp EM_separate.cpp simulation.cpp)

add_compile_options(-Wall -O3 -ffast-math)
add_library(PROBer_core ${sources} ${headers})

add_executable(PROBer-extract-reference-transcripts extractRef.cpp)
add_executable(PROBer-synthesis-reference-transcripts synthesisRef.cpp)
add_executable(PROBer-preref preRef.cpp)
add_executable(PROBer-parse-alignments parseAlignments.cpp)
add_dependencies(PROBer-parse-alignments samtools)
add_executable(PROBer-run-em EM.cpp)
add_dependencies(PROBer-run-em samtools)
add_executable(PROBer-simulate-reads simulation.cpp)
add_dependencies(PROBer-simulate-reads samtools)

target_link_libraries(PROBer-extract-reference-transcripts PROBer_core)
target_link_libraries(PROBer-synthesis-reference-transcripts PROBer_core)
target_link_libraries(PROBer-preref PROBer_core)
target_link_libraries(PROBer-parse-alignments PROBer_core ${SAMTOOLS_DIR}/libbam.a ${ZLIB_LIBRARIES} pthread)
target_link_libraries(PROBer-run-em PROBer_core ${SAMTOOLS_DIR}/libbam.a ${ZLIB_LIBRARIES} pthread)
target_link_libraries(PROBer-simulate-reads PROBer_core ${SAMTOOLS_DIR}/libbam.a ${ZLIB_LIBRARIES} pthread)


configure_file(PROBer PROBer COPYONLY)


install(TARGETS PROBer-extract-reference-transcripts PROBer-synthesis-reference-transcripts PROBer-preref PROBer-parse-alignments PROBer-run-em PROBer-simulate-reads DESTINATION bin)
install(PROGRAMS PROBer DESTINATION bin)



add_executable(PROBer-run-em-separate EM_separate.cpp)
add_dependencies(PROBer-run-em-separate samtools)
target_link_libraries(PROBer-run-em-separate PROBer_core ${SAMTOOLS_DIR}/libbam.a ${ZLIB_LIBRARIES} pthread)
configure_file(PROBer-estimate-parameters-separately PROBer-estimate-parameters-separately COPYONLY)
install(TARGETS PROBer-run-em-separate DESTINATION bin)
install(PROGRAMS PROBer-estimate-parameters-separately DESTINATION bin)
