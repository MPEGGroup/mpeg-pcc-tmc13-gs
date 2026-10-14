cmake_minimum_required(VERSION 3.0)
project (TMC3 CXX)

# Note: This is a temporary macro. It must be removed after the conformance bitstream is regenerated.
option(TEST_CONFORMANCE "Conformance testing" ON)

include( CheckCXXCompilerFlag )
CHECK_CXX_COMPILER_FLAG( "-std=c++17" COMPILER_SUPPORTS_CXX17 )
CHECK_CXX_COMPILER_FLAG( "-std=c++0x" COMPILER_SUPPORTS_CXX0X )

if (COMPILER_SUPPORTS_CXX17)
 set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17" )
elseif (COMPILER_SUPPORTS_CXX0X)
 set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x" )
else ()
 MESSAGE(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++17 support. Please use a different C++ compiler.")
endif ()

if(EMSCRIPTEN)
  set(WASM_MODE 1)
endif ()

add_subdirectory(tmc3)
