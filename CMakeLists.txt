cmake_minimum_required(VERSION 3.28)
find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )
project(surfex VERSION 1 LANGUAGES Fortran)

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
	list(APPEND PRECISIONS "dp")
endif()
if(${HAVE_SINGLE_PRECISION})
	list(APPEND PRECISIONS "sp")
endif()

ecbuild_find_package(fiat COMPONENTS REQUIRED)
ecbuild_find_package(NetCDF COMPONENTS Fortran REQUIRED)
ecbuild_find_package(eccodes REQUIRED)
ecbuild_find_package(MPI REQUIRED)
ecbuild_find_package(falfilfa REQUIRED)

foreach(PRECISION ${PRECISIONS})
	set(LIBNAME surfex_${PRECISION})
	message(STATUS "Will build library ${LIBNAME}")
	ecbuild_add_library(
		TARGET ${LIBNAME}
		SOURCES ${SRC} ${MODI_LIST}
		TYPE STATIC
		PUBLIC_LIBS fiat
		MPI::MPI_Fortran
		NetCDF::NetCDF_Fortran
		eccodes
		parkind_${PRECISION}
		fa_${PRECISION}
		PUBLIC_INCLUDES
		$<BUILD_INTERFACE:${ECCODES_INCLUDE_DIR}>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
		$<INSTALL_INTERFACE:include/${PRECISION}>
		$<INSTALL_INTERFACE:module/${PRECISION}>
	)

	if(HAVE_SINGLE_PRECISION AND HAVE_DOUBLE_PRECISION)
		ecbuild_info("both precisions are activated code will be erased between the two compilations")
		set_target_properties(${LIBNAME} PROPERTIES CLEAN_DIRECT_OUTPUT 1)
		set_target_properties(${LIBNAME} PROPERTIES EXCLUDE_FROM_ALL FALSE)
	endif()
	set_target_properties(${LIBNAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/${PRECISION})
	target_include_directories(${LIBNAME} PUBLIC
		$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/module/${PRECISION}>
	)

	target_compile_definitions(${LIBNAME} PUBLIC in_surfex)
	target_compile_definitions(${LIBNAME} PUBLIC SFX_ARO)
	target_compile_definitions(${LIBNAME} PUBLIC SFX_ASC)
	target_compile_definitions(${LIBNAME} PUBLIC SFX_OL)
	target_compile_definitions(${LIBNAME} PUBLIC SFX_TXT)
	target_compile_definitions(${LIBNAME} PUBLIC SFX_FA)
	#target_compile_definitions(${LIBNAME} PUBLIC SFX_LFI)
	target_compile_definitions(${LIBNAME} PUBLIC ARO)
	target_compile_definitions(${LIBNAME} PUBLIC OL)
	target_compile_definitions(${LIBNAME} PUBLIC ASC)
	target_compile_definitions(${LIBNAME} PUBLIC TXT)
	target_compile_definitions(${LIBNAME} PUBLIC FA)
	target_compile_definitions(${LIBNAME} PUBLIC LFI)
endforeach()

if(${HAVE_DOUBLE_PRECISION})
	target_compile_options(surfex_dp PUBLIC "-r8" )
endif()

if(${HAVE_SINGLE_PRECISION})
	target_compile_definitions(surfex_sp PUBLIC PARKIND1_SINGLE)
endif()

foreach(FILEPATH ${SRC_PROGRAMS})
	foreach(PRECISION ${PRECISIONS})
		message(STATUS "Will build program ${EXENAME}")
		cmake_path(GET FILEPATH STEM STEM)
		set(EXENAME ${STEM}_${PRECISION})

		ecbuild_add_executable(TARGET ${EXENAME} SOURCES ${FILEPATH})
		target_link_libraries(${EXENAME}
			PUBLIC
			fiat
			MPI::MPI_Fortran
			eccodes
			NetCDF::NetCDF_Fortran
			surfex_${PRECISION}
			parkind_${PRECISION}
			fa_${PRECISION}
		)

		if(${PRECISION} STREQUAL "dp")
			target_compile_options(${EXENAME} PUBLIC "-r8" )
		endif()
		if(${PRECISION} STREQUAL "sp")
			target_compile_definitions(${EXENAME} PUBLIC PARKIND1_SINGLE)
		endif()
	endforeach()
endforeach()

if(${HAVE_DOUBLE_PRECISION})
	install(TARGETS surfex_dp LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
	install(DIRECTORY ${CMAKE_BINARY_DIR}/module/dp DESTINATION module/)
endif()

if(${HAVE_SINGLE_PRECISION})
	install(TARGETS surfex_sp LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
	install(DIRECTORY ${CMAKE_BINARY_DIR}/module/sp DESTINATION module/)
endif()

ecbuild_install_project(NAME surfex)
ecbuild_print_summary()
