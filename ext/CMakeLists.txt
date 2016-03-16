include(ExternalProject)

ExternalProject_Add(
	HTSlib
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/htslib-1.3
	CONFIGURE_COMMAND ""
	BUILD_COMMAND make
	BUILD_IN_SOURCE 1
	INSTALL_COMMAND ""
)

ExternalProject_Get_Property(HTSlib SOURCE_DIR)
set(HTSLIB_DIR ${SOURCE_DIR} CACHE INTERNAL "Path to HTSlib")
