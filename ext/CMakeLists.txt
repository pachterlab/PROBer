include(ExternalProject)

ExternalProject_Add(
	samtools
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/sam
	CONFIGURE_COMMAND ""
	BUILD_COMMAND make
	BUILD_IN_SOURCE 1
	INSTALL_COMMAND ""
)

ExternalProject_Get_Property(samtools SOURCE_DIR)
set(SAMTOOLS_DIR ${SOURCE_DIR} CACHE INTERNAL "Path to Samtools")
