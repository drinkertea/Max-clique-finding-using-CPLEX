
if(WIN32)
    execute_process(COMMAND cmd /C set CPLEX_STUDIO_DIR OUTPUT_VARIABLE CPLEX_STUDIO_DIR_VAR ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(NOT CPLEX_STUDIO_DIR_VAR)
        MESSAGE(FATAL_ERROR "Unable to find CPLEX: environment variable CPLEX_STUDIO_DIR<VERSION> not set.")
    endif()

    STRING(REGEX REPLACE "^CPLEX_STUDIO_DIR" "" CPLEX_STUDIO_DIR_VAR ${CPLEX_STUDIO_DIR_VAR})
    STRING(REGEX MATCH "^[0-9]+" CPLEX_WIN_VERSION ${CPLEX_STUDIO_DIR_VAR})
    STRING(REGEX REPLACE "^[0-9]+=" "" CPLEX_STUDIO_DIR_VAR ${CPLEX_STUDIO_DIR_VAR})
    file(TO_CMAKE_PATH "${CPLEX_STUDIO_DIR_VAR}" CPLEX_ROOT_DIR_GUESS)

    set(CPLEX_WIN_VERSION "${CPLEX_WIN_VERSION}" CACHE STRING "CPLEX version to be used.")
    set(CPLEX_ROOT_DIR "${CPLEX_ROOT_DIR_GUESS}" CACHE PATH "CPLEX root directory.")

    MESSAGE(STATUS "Found CLPEX version ${CPLEX_WIN_VERSION} at '${CPLEX_ROOT_DIR}'")


    # For release build set(CPLEX_WIN_PLATFORM "x64_windows_msvc14/stat_mda")
    set(CPLEX_WIN_PLATFORM "x64_windows_msvc14/stat_mdd")

else()

    set(CPLEX_ROOT_DIR "" CACHE PATH "CPLEX root directory.")
    set(CPLEX_WIN_PLATFORM "")

endif()


FIND_PATH(CPLEX_INCLUDE_DIR
        ilcplex/cplex.h
        HINTS ${CPLEX_ROOT_DIR}/cplex/include
        ${CPLEX_ROOT_DIR}/include
        PATHS ENV C_INCLUDE_PATH
        ENV C_PLUS_INCLUDE_PATH
        ENV INCLUDE_PATH
        )

message(STATUS "CPLEX INCLUDE: ${CPLEX_INCLUDE_DIR}")


FIND_PATH(CPLEX_CONCERT_INCLUDE_DIR
        ilconcert/iloenv.h
        HINTS ${CPLEX_ROOT_DIR}/concert/include
        ${CPLEX_ROOT_DIR}/include
        PATHS ENV C_INCLUDE_PATH
        ENV C_PLUS_INCLUDE_PATH
        ENV INCLUDE_PATH
        )

message(STATUS "concert INCLUDE: ${CPLEX_CONCERT_INCLUDE_DIR}")


FIND_LIBRARY(CPLEX_LIBRARY
        NAMES cplex${CPLEX_WIN_VERSION}0 cplex
        HINTS ${CPLEX_ROOT_DIR}/cplex/lib/${CPLEX_WIN_PLATFORM} #windows
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_debian4.0_4.1/static_pic #unix
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_sles10_4.1/static_pic #unix 
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_linux/static_pic #unix 
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_osx/static_pic #osx 
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_darwin/static_pic #osx 
        PATHS ENV LIBRARY_PATH #unix
        ENV LD_LIBRARY_PATH #unix
        )
message(STATUS "CPLEX Library: ${CPLEX_LIBRARY}")

FIND_LIBRARY(CPLEX_ILOCPLEX_LIBRARY
        ilocplex
        HINTS ${CPLEX_ROOT_DIR}/cplex/lib/${CPLEX_WIN_PLATFORM} #windows
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_debian4.0_4.1/static_pic #unix 
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_sles10_4.1/static_pic #unix 
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_linux/static_pic #unix 
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_osx/static_pic #osx 
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_darwin/static_pic #osx 
        PATHS ENV LIBRARY_PATH
        ENV LD_LIBRARY_PATH
        )
message(STATUS "ILOCPLEX Library: ${CPLEX_ILOCPLEX_LIBRARY}")

FIND_LIBRARY(CPLEX_CONCERT_LIBRARY
        concert
        HINTS ${CPLEX_ROOT_DIR}/concert/lib/${CPLEX_WIN_PLATFORM} #windows
        ${CPLEX_ROOT_DIR}/concert/lib/x86-64_debian4.0_4.1/static_pic #unix 
        ${CPLEX_ROOT_DIR}/concert/lib/x86-64_sles10_4.1/static_pic #unix 
        ${CPLEX_ROOT_DIR}/concert/lib/x86-64_linux/static_pic #unix 
        ${CPLEX_ROOT_DIR}/concert/lib/x86-64_osx/static_pic #osx 
        ${CPLEX_ROOT_DIR}/concert/lib/x86-64_darwin/static_pic #osx 
        PATHS ENV LIBRARY_PATH
        ENV LD_LIBRARY_PATH
        )
message(STATUS "CONCERT Library: ${CPLEX_CONCERT_LIBRARY}")

if(WIN32)
    FIND_PATH(CPLEX_BIN_DIR
            cplex${CPLEX_WIN_VERSION}0.dll
            HINTS ${CPLEX_ROOT_DIR}/cplex/bin/x64_win64 #windows
            )
else()
    FIND_PATH(CPLEX_BIN_DIR
            cplex
            HINTS ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_sles10_4.1 #unix
            ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_debian4.0_4.1 #unix
            ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_linux #unix
            ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_osx #osx
            ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_darwin #osx
            ENV LIBRARY_PATH
            ENV LD_LIBRARY_PATH
            )
endif()
message(STATUS "CPLEX Bin Dir: ${CPLEX_BIN_DIR}")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CPLEX DEFAULT_MSG
        CPLEX_LIBRARY CPLEX_INCLUDE_DIR CPLEX_ILOCPLEX_LIBRARY CPLEX_CONCERT_LIBRARY CPLEX_CONCERT_INCLUDE_DIR)

IF(CPLEX_FOUND)
    SET(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR} ${CPLEX_CONCERT_INCLUDE_DIR})
    SET(CPLEX_LIBRARIES ${CPLEX_CONCERT_LIBRARY} ${CPLEX_ILOCPLEX_LIBRARY} ${CPLEX_LIBRARY} )
    IF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        SET(CPLEX_LIBRARIES "${CPLEX_LIBRARIES};m;pthread")
    ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
ENDIF(CPLEX_FOUND)

MARK_AS_ADVANCED(CPLEX_LIBRARY CPLEX_INCLUDE_DIR CPLEX_ILOCPLEX_LIBRARY CPLEX_CONCERT_INCLUDE_DIR CPLEX_CONCERT_LIBRARY)