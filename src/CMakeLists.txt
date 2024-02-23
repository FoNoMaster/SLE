if (USE_GLOB)
    file(GLOB_RECURSE SOURCES *.hpp)
    add_library(${PROJECT_NAME} INTERFACE ${SOURCES}) 
else () 
	set(HEADERS
		code/CSR.hpp
		code/dense.hpp
		code/vector.hpp
	)

	set(SOURCES
        )

    source_group("Header Files" FILES ${HEADERS})
    source_group("Source Files" FILES ${SOURCES})

    add_library(${PROJECT_NAME} INTERFACE ${HEADERS} ${SOURCES})
endif ()

target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)
