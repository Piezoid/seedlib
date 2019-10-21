cmake_minimum_required(VERSION 3.5)

project(seedlib
    VERSION 0.0.1
    LANGUAGES CXX)
set(PROJECT_HOMEPAGE_URL "https://github.com/Piezoid/seedlib")
set(PROJECT_DESCRIPTION "An index for approximate 1-lev mapping of genomic sequences")
set(PROJECT_VENDOR "Maël Kerbiriou")
set(PROJECT_CONTACT "m431.kerbiriou@gmail.com")

set(CMAKE_CXX_STANDARD 11)

set(thirdparty_dir "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/")
set(include_dir "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(src_dir "${CMAKE_CURRENT_SOURCE_DIR}/src")

include(GNUInstallDirs)
add_subdirectory(thirdparty)

add_library(${PROJECT_NAME} STATIC "src/seedlib.cpp")
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 14)
target_link_libraries(${PROJECT_NAME} PRIVATE gatbl sdsl)
target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<BUILD_INTERFACE:${include_dir}>
    PRIVATE
        ${src_dir}
)
install(DIRECTORY ${include_dir}/seedlib DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

set(BINARY_TARGET "${PROJECT_NAME}-bin")
add_executable(${BINARY_TARGET} "${src_dir}/main.cpp")
set_target_properties(${BINARY_TARGET}
    PROPERTIES
        OUTPUT_NAME ${PROJECT_NAME}
)
target_include_directories(${BINARY_TARGET}
    PRIVATE
        ${src_dir}
        ${include_dir}
)
target_link_libraries(${BINARY_TARGET} PRIVATE ${PROJECT_NAME} gatbl clipp)


install(TARGETS seedlib seedlib-bin
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	)
