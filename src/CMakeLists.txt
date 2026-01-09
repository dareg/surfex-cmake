cmake_minimum_required(VERSION 3.24)
find_package( ecbuild 3.11 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )
project(surfex VERSION 8.1.502 LANGUAGES Fortran)

set(MODI_PATH ${CMAKE_BINARY_DIR}/modi)
file(MAKE_DIRECTORY ${MODI_PATH})

include(set_sources.cmake)

foreach(FILEPATH ${SRC})
  cmake_path(GET FILEPATH FILENAME FILENAME)
  set(OUTPUT_FILE ${MODI_PATH}/modi_${FILENAME})
  add_custom_command(
    OUTPUT ${OUTPUT_FILE}
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/gen_modi.py -o ${MODI_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/${FILEPATH}
    DEPENDS ${FILEPATH}
  )
  list(APPEND MODI_LIST ${OUTPUT_FILE})
endforeach()

ecbuild_add_option( FEATURE DOUBLE_PRECISION
  DEFAULT ON
  DESCRIPTION "Support for double precision" )
ecbuild_add_option( FEATURE SINGLE_PRECISION
  DEFAULT OFF
  DESCRIPTION "Support for single precision" )
ecbuild_add_option( FEATURE SFX_MPI
  DEFAULT OFF
  DESCRIPTION "Activate the SFX_MPI preprocessor macro" )

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
ecbuild_find_package(pinuts REQUIRED)
find_package(LAPACK REQUIRED)

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
        LAPACK::LAPACK
        eccodes
        parkind_${PRECISION}
        fa_${PRECISION}
        pinuts_${PRECISION}
    PUBLIC_INCLUDES
        $<BUILD_INTERFACE:${ECCODES_INCLUDE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/${LIBNAME}>
        $<INSTALL_INTERFACE:module/${LIBNAME}>
  )

  if(HAVE_SINGLE_PRECISION AND HAVE_DOUBLE_PRECISION)
    ecbuild_info("both precisions are activated code will be erased between the two compilations")
    set_target_properties(${LIBNAME} PROPERTIES CLEAN_DIRECT_OUTPUT 1)
    set_target_properties(${LIBNAME} PROPERTIES EXCLUDE_FROM_ALL FALSE)
  endif()
  set_target_properties(${LIBNAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/${LIBNAME})
  target_include_directories(${LIBNAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/module/${LIBNAME}>
  )

  if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )   # TODO: CLEAN 
    if( CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL 10 )
          ecbuild_add_fortran_flags( -fallow-argument-mismatch NO_FAIL )
    endif()
  endif()

  target_compile_definitions(${LIBNAME} PRIVATE in_surfex)
  target_compile_definitions(${LIBNAME} PRIVATE SFX_ARO)
  target_compile_definitions(${LIBNAME} PRIVATE SFX_ASC)
  target_compile_definitions(${LIBNAME} PRIVATE SFX_OL)
  target_compile_definitions(${LIBNAME} PRIVATE SFX_TXT)
  target_compile_definitions(${LIBNAME} PRIVATE SFX_FA)
  target_compile_definitions(${LIBNAME} PRIVATE SFX_LFI)
  target_compile_definitions(${LIBNAME} PRIVATE ARO)
  target_compile_definitions(${LIBNAME} PRIVATE OL)
  target_compile_definitions(${LIBNAME} PRIVATE ASC)
  target_compile_definitions(${LIBNAME} PRIVATE TXT)
  target_compile_definitions(${LIBNAME} PRIVATE FA)
  target_compile_definitions(${LIBNAME} PRIVATE LFI)
  if(${HAVE_SFX_MPI})
    target_compile_definitions(${LIBNAME} PRIVATE SFX_MPI)
  endif()
endforeach()

if(${HAVE_DOUBLE_PRECISION})
  target_compile_options(surfex_dp PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${ECBUILD_Fortran_COMPILE_OPTIONS_REAL8}>)
endif()

if(${HAVE_SINGLE_PRECISION})
  target_compile_options(surfex_sp PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${ECBUILD_Fortran_COMPILE_OPTIONS_REAL4}>)
  target_compile_definitions(surfex_sp PUBLIC PARKIND1_SINGLE)
endif()

foreach(FILEPATH ${SRC_PROGRAMS})
  foreach(PRECISION ${PRECISIONS})
    message(STATUS "Will build program ${EXENAME}")
    cmake_path(GET FILEPATH STEM STEM)
    set(EXENAME ${STEM}_${PRECISION})

    ecbuild_add_executable(TARGET ${EXENAME} SOURCES ${FILEPATH}  LINKER_LANGUAGE  Fortran)
    target_link_libraries(${EXENAME}
      PUBLIC
      fiat
      MPI::MPI_Fortran
      eccodes
      NetCDF::NetCDF_Fortran
      LAPACK::LAPACK
      surfex_${PRECISION}
      parkind_${PRECISION}
      fa_${PRECISION}
      pinuts_${PRECISION}
    )

    if(${PRECISION} STREQUAL "dp")
      target_compile_options(${EXENAME} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${ECBUILD_Fortran_COMPILE_OPTIONS_REAL8}>)
    endif()
    if(${PRECISION} STREQUAL "sp")
      target_compile_options(${EXENAME} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${ECBUILD_Fortran_COMPILE_OPTIONS_REAL4}>)
      target_compile_definitions(${EXENAME} PUBLIC PARKIND1_SINGLE)
    endif()

  target_compile_definitions(${EXENAME} PRIVATE SFX_FA)
  if(${HAVE_SFX_MPI})
    target_compile_definitions(${EXENAME} PRIVATE SFX_MPI)
  endif()

  endforeach()
endforeach()

if(${HAVE_DOUBLE_PRECISION})
  install(TARGETS surfex_dp LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY ${CMAKE_BINARY_DIR}/module/surfex_dp DESTINATION module/)
endif()

if(${HAVE_SINGLE_PRECISION})
  install(TARGETS surfex_sp LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY ${CMAKE_BINARY_DIR}/module/surfex_sp DESTINATION module/)
endif()

ecbuild_install_project(NAME surfex)
ecbuild_print_summary()
