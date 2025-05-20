cmake_minimum_required(VERSION 3.28)
find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )
project(surfex VERSION 1 LANGUAGES Fortran)
#ecbuild_enable_fortran(REQUIRED NO_MODULE_DIRECTORY)
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
	DESCRIPTION "Support for double precision" )
ecbuild_add_option( FEATURE SINGLE_PRECISION
	DEFAULT OFF
	DESCRIPTION "Support for single precision" )

if(${HAVE_DOUBLE_PRECISION})
	list(APPEND PRECISIONS "_dp")
endif()
if(${HAVE_SINGLE_PRECISION})
	list(APPEND PRECISIONS "_sp")
endif()

ecbuild_find_package(fiat COMPONENTS REQUIRED)
ecbuild_find_package(NetCDF COMPONENTS Fortran REQUIRED)
ecbuild_find_package(eccodes REQUIRED)
ecbuild_find_package(MPI REQUIRED)
ecbuild_find_package(falfilfa REQUIRED)

foreach(PREC ${PRECISIONS})
	message(STATUS "building library surfex${PREC}")
	ecbuild_add_library(
		TARGET surfex${PREC}
		SOURCES ${SRC} ${MODI_LIST}
		TYPE STATIC
		PUBLIC_LIBS fiat
		MPI::MPI_Fortran
		NetCDF::NetCDF_Fortran
		eccodes
		PUBLIC_INCLUDES
		$<BUILD_INTERFACE:${ECCODES_INCLUDE_DIR}>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/includes>
		$<INSTALL_INTERFACE:include/${PREC}>
		$<INSTALL_INTERFACE:modules/${PREC}>
	)
	if(HAVE_SINGLE_PRECISION AND HAVE_DOUBLE_PRECISION)
		ecbuild_info("both precisions are activated code will be erased between the two compilations")
		set_target_properties(surfex${PREC} PROPERTIES CLEAN_DIRECT_OUTPUT 1)
		set_target_properties(surfex${PREC} PROPERTIES EXCLUDE_FROM_ALL FALSE)
	endif()
	set_target_properties(surfex${PREC} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module${PREC})
	target_include_directories(surfex${PREC} PUBLIC
		$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/module${PREC}>
	)

	target_compile_definitions(surfex${PREC} PUBLIC in_surfex)
	target_compile_definitions(surfex${PREC} PUBLIC SFX_ARO)
	target_compile_definitions(surfex${PREC} PUBLIC SFX_ASC)
	target_compile_definitions(surfex${PREC} PUBLIC SFX_OL)
	target_compile_definitions(surfex${PREC} PUBLIC SFX_TXT)
	target_compile_definitions(surfex${PREC} PUBLIC SFX_FA)
	#target_compile_definitions(surfex${PREC} PUBLIC SFX_LFI)
	target_compile_definitions(surfex${PREC} PUBLIC ARO)
	target_compile_definitions(surfex${PREC} PUBLIC OL)
	target_compile_definitions(surfex${PREC} PUBLIC ASC)
	target_compile_definitions(surfex${PREC} PUBLIC TXT)
	target_compile_definitions(surfex${PREC} PUBLIC FA)
	target_compile_definitions(surfex${PREC} PUBLIC LFI)
endforeach()

if(${HAVE_DOUBLE_PRECISION})
	target_compile_options(surfex_dp PUBLIC "-r8" )
	target_link_libraries(surfex_dp PUBLIC parkind_dp fa_dp)
endif()

if(${HAVE_SINGLE_PRECISION})
	target_link_libraries(surfex_sp PUBLIC parkind_sp fa_sp)
	target_compile_definitions(surfex_sp PUBLIC PARKIND1_SINGLE)
endif()

foreach(FILEPATH ${SRC_PROGRAMS})
foreach(PREC ${PRECISIONS})
	cmake_path(GET FILEPATH STEM STEM)
	set(EXENAME ${STEM}${PREC})
	ecbuild_add_executable(TARGET ${EXENAME} SOURCES ${FILEPATH})
	target_link_libraries(${EXENAME} PUBLIC fiat MPI::MPI_Fortran eccodes NetCDF::NetCDF_Fortran)
	message(STATUS "Added exe ${EXENAME}")
	if(${PREC} STREQUAL "_dp")
		target_compile_options(${EXENAME} PUBLIC "-r8" )
		target_link_libraries(${EXENAME} PUBLIC surfex_dp parkind_dp fa_dp)
	endif()
	if(${PREC} STREQUAL "_sp")
		target_link_libraries(${EXENAME} PUBLIC surfex_sp parkind_sp fa_sp)
		target_compile_definitions(${EXENAME} PUBLIC PARKIND1_SINGLE)
	endif()
endforeach()
endforeach()

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION .)
if(${HAVE_DOUBLE_PRECISION})
install(TARGETS surfex_dp LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
if(${HAVE_SINGLE_PRECISION})
install(TARGETS surfex_sp LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()


ecbuild_install_project(NAME surfex)
ecbuild_print_summary()
