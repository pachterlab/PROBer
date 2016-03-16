include_directories(${ZLIB_INCLUDE_DIRS})
include_directories(../ext)

file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)
list(REMOVE_ITEM sources extractRef.cpp synthesisRef.cpp preRef.cpp parseAlignments.cpp EM.cpp simulation.cpp)

add_compile_options(-Wall -O3 -ffast-math)
add_library(PROBer_core ${sources} ${headers})

add_executable(PROBer-extract-reference-transcripts extractRef.cpp)
add_executable(PROBer-synthesis-reference-transcripts synthesisRef.cpp)
add_executable(PROBer-preref preRef.cpp)
add_executable(PROBer-parse-alignments parseAlignments.cpp)
add_dependencies(PROBer-parse-alignments HTSlib)
add_executable(PROBer-run-em EM.cpp)
add_dependencies(PROBer-run-em HTSlib)
add_executable(PROBer-simulate-reads simulation.cpp)
add_dependencies(PROBer-simulate-reads HTSlib)

target_link_libraries(PROBer-extract-reference-transcripts PROBer_core)
target_link_libraries(PROBer-synthesis-reference-transcripts PROBer_core)
target_link_libraries(PROBer-preref PROBer_core)
target_link_libraries(PROBer-parse-alignments PROBer_core ${HTSLIB_DIR}/libhts.a ${ZLIB_LIBRARIES} pthread)
target_link_libraries(PROBer-run-em PROBer_core ${HTSLIB_DIR}/libhts.a ${ZLIB_LIBRARIES} pthread)
target_link_libraries(PROBer-simulate-reads PROBer_core ${HTSLIB_DIR}/libhts.a ${ZLIB_LIBRARIES} pthread)


configure_file(PROBer PROBer COPYONLY)


install(TARGETS PROBer-extract-reference-transcripts PROBer-synthesis-reference-transcripts PROBer-preref PROBer-parse-alignments PROBer-run-em PROBer-simulate-reads DESTINATION bin)
install(PROGRAMS PROBer DESTINATION bin)
