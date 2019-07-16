cmake_minimum_required(VERSION 3.5)

project(seedlib
    VERSION 0.0.1
    LANGUAGES CXX)
set(PROJECT_HOMEPAGE_URL "https://github.com/Piezoid/seedlib")
set(PROJECT_DESCRIPTION "An genomic sequences index allowing query with 1 error")
set(PROJECT_VENDOR "Maël Kerbiriou")
set(PROJECT_CONTACT "mael.kerbiriou@free.fr")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


set(thirdparty_dir "thirdparty/")

find_library(SDSL_LIBRARY NAMES sdsl sdsl-lite)
if(NOT SDSL_LIBRARY)
    FILE(TO_CMAKE_PATH "${thirdparty_dir}/sdsl-lite" sdsl_project_dir)
    if(NOT EXISTS "${sdsl_project_dir}/CMakeLists.txt")
        execute_process(COMMAND git submodule update --init -- ${sdsl_project_dir}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	endif()
    set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)
    FILE(TO_CMAKE_PATH "${sdsl_project_dir}/include" SDSL_INCLUDE_DIR)
	FILE(TO_CMAKE_PATH "${sdsl_project_dir}/lib" SDSL_LIB_DIR)
endif()

find_path(GATBL_INCLUDE_DIR "common.hpp" PATHS gatbl gatb-lite ENV GATBL_INCLUDE_DIR)
if(NOT GATBL_INCLUDE_DIR)
    FILE(TO_CMAKE_PATH "${thirdparty_dir}/gatb-lite" gatbl_project_dir)
    if(NOT EXISTS "${gatbl_project_dir}/CMakeLists.txt")
        execute_process(COMMAND git submodule update --init -- ${gatbl_project_dir}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
    FILE(TO_CMAKE_PATH "${gatbl_project_dir}/include" GATBL_INCLUDE_DIR)
endif()

find_path(CLIPP_INCLUDE_DIR "clipp.h" ENV CLIPP_INCLUDE_DIR)
if(NOT CLIPP_INCLUDE_DIR)
    FILE(TO_CMAKE_PATH "${thirdparty_dir}/clipp" clipp_project_dir)
    if(NOT EXISTS "${clipp_project_dir}/README.md")
        execute_process(COMMAND git submodule update --init -- ${clipp_project_dir}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
    FILE(TO_CMAKE_PATH "${clipp_project_dir}/include" CLIPP_INCLUDE_DIR)
endif()

include_directories(SYSTEM ${CLIPP_INCLUDE_DIR})
include_directories(SYSTEM ${SDSL_INCLUDE_DIR})
include_directories(${GATBL_INCLUDE_DIR})

set( sdsl_SRCS
	"${SDSL_LIB_DIR}/bits.cpp"
	"${SDSL_LIB_DIR}/util.cpp"
	"${SDSL_LIB_DIR}/ram_fs.cpp"
	"${SDSL_LIB_DIR}/memory_management.cpp"
	"${SDSL_LIB_DIR}/rrr_vector_15.cpp"
)

add_executable(${PROJECT_NAME} "main.cpp" "indexer.cpp" ${sdsl_SRCS})
