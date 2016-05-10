include_directories(${ZLIB_INCLUDE_DIRS})
include_directories(../ext)
include_directories(BEFORE ${HTSLIB_DIR})

file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)
list(REMOVE_ITEM sources buildRef.cpp parseAlignments.cpp EM.cpp simulation.cpp analyze_iCLIP.cpp PROBer_single_transcript_batch.cpp)

add_compile_options(-Wall -O3 -ffast-math)
add_library(PROBer_core ${sources} ${headers})

add_executable(PROBer-build-reference buildRef.cpp)
add_executable(PROBer-parse-alignments parseAlignments.cpp)
add_dependencies(PROBer-parse-alignments HTSlib)
add_executable(PROBer-run-em EM.cpp)
add_dependencies(PROBer-run-em HTSlib)
add_executable(PROBer-simulate-reads simulation.cpp)
add_dependencies(PROBer-simulate-reads HTSlib)
add_executable(PROBer-analyze-iCLIP analyze_iCLIP.cpp)
add_dependencies(PROBer-analyze-iCLIP HTSlib)

add_executable(PROBer-single-transcript-batch PROBer_single_transcript_batch.cpp)
add_dependencies(PROBer-single-transcript-batch HTSlib)

target_link_libraries(PROBer-build-reference PROBer_core)
target_link_libraries(PROBer-parse-alignments PROBer_core ${HTSLIB_DIR}/libhts.a ${ZLIB_LIBRARIES} pthread)
target_link_libraries(PROBer-run-em PROBer_core ${HTSLIB_DIR}/libhts.a ${ZLIB_LIBRARIES} pthread)
target_link_libraries(PROBer-simulate-reads PROBer_core ${HTSLIB_DIR}/libhts.a ${ZLIB_LIBRARIES} pthread)
target_link_libraries(PROBer-analyze-iCLIP PROBer_core ${HTSLIB_DIR}/libhts.a ${ZLIB_LIBRARIES} pthread)

target_link_libraries(PROBer-single-transcript-batch PROBer_core ${HTSLIB_DIR}/libhts.a ${ZLIB_LIBRARIES} pthread)
configure_file(utils.py utils.py COPYONLY)
configure_file(PROBer-single-batch-estimate PROBer-single-batch-estimate COPYONLY)

configure_file(PROBer PROBer COPYONLY)

install(TARGETS PROBer-build-reference PROBer-parse-alignments PROBer-run-em PROBer-simulate-reads PROBer-analyze-iCLIP DESTINATION bin)
install(PROGRAMS PROBer DESTINATION bin)
