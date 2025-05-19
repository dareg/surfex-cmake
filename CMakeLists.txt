cmake_minimum_required(VERSION 3.28)
find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )
project(surfex VERSION 1 LANGUAGES Fortran)
ecbuild_enable_fortran(REQUIRED NO_MODULE_DIRECTORY)
set(MODI_PATH ${CMAKE_BINARY_DIR}/modi)
file(MAKE_DIRECTORY ${MODI_PATH})

include(set_sources.cmake)

foreach(FILEPATH ${SRC})
	cmake_path(GET FILEPATH FILENAME FILENAME)
	set(OUTPUT_FILE ${MODI_PATH}/modi_${FILENAME})
	add_custom_command(
		OUTPUT ${OUTPUT_FILE}
		COMMAND ${CMAKE_SOURCE_DIR}/gen_modi.py -o ${MODI_PATH} ${CMAKE_SOURCE_DIR}/${FILEPATH}
		DEPENDS ${FILEPATH}
	)
	list(APPEND MODI_LIST ${OUTPUT_FILE})
endforeach()

ecbuild_add_option( FEATURE DOUBLE_PRECISION
	DEFAULT OFF
	DESCRIPTION "Support for Double Precision" )

ecbuild_find_package(fiat COMPONENTS REQUIRED)
ecbuild_find_package(NetCDF COMPONENTS Fortran REQUIRED)
ecbuild_find_package(eccodes REQUIRED)
ecbuild_find_package(MPI REQUIRED)
ecbuild_find_package(falfilfa REQUIRED)

ecbuild_add_library(
	TARGET surfex
	SOURCES ${SRC} ${MODI_LIST}
	TYPE STATIC
	PUBLIC_LIBS fiat
		MPI::MPI_Fortran
		NetCDF::NetCDF_Fortran
		eccodes
	)

if(${HAVE_DOUBLE_PRECISION})
	ecbuild_add_fortran_flags( "-r8" )
	target_link_libraries(surfex PUBLIC parkind_dp fa_dp)
else()
	target_link_libraries(surfex PUBLIC parkind_sp fa_sp)
	target_compile_definitions(surfex PUBLIC PARKIND1_SINGLE)
endif()

target_compile_definitions(surfex PUBLIC in_surfex)
target_compile_definitions(surfex PUBLIC SFX_ARO)
target_compile_definitions(surfex PUBLIC SFX_ASC)
target_compile_definitions(surfex PUBLIC SFX_OL)
target_compile_definitions(surfex PUBLIC SFX_TXT)
target_compile_definitions(surfex PUBLIC SFX_FA)
#target_compile_definitions(surfex PUBLIC SFX_LFI)
target_compile_definitions(surfex PUBLIC ARO)
target_compile_definitions(surfex PUBLIC OL)
target_compile_definitions(surfex PUBLIC ASC)
target_compile_definitions(surfex PUBLIC TXT)
target_compile_definitions(surfex PUBLIC FA)
target_compile_definitions(surfex PUBLIC LFI)

foreach(FILEPATH ${SRC_PROGRAMS})
	cmake_path(GET FILEPATH STEM EXENAME)
	ecbuild_add_executable(TARGET ${EXENAME} SOURCES ${FILEPATH})
	target_link_libraries( ${EXENAME} PUBLIC surfex fiat MPI::MPI_Fortran eccodes NetCDF::NetCDF_Fortran)
	if(${HAVE_DOUBLE_PRECISION})
		ecbuild_add_fortran_flags( "-r8" )
		target_link_libraries(${EXENAME} PUBLIC parkind_dp fa_dp)
	else()
		target_link_libraries(${EXENAME} PUBLIC parkind_sp fa_sp)
		target_compile_definitions(${EXENAME} PUBLIC PARKIND1_SINGLE)
	endif()
	message(STATUS "Added exe ${EXENAME}")
endforeach()

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION .)
install(TARGETS surfex LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})


ecbuild_install_project(NAME surfex)
ecbuild_print_summary()
