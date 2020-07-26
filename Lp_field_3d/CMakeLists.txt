CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

PROJECT( lp_field )

set( lp_field_SRCS
     src/GeometryUtils.cpp
     src/IO.cpp
     src/Distance.cpp
     src/predicates.cpp
     src/LinearAlgebra.cpp
     src/TriangleIntegral.cpp
     src/TrianglePotential.cpp
)


FIND_PACKAGE (OpenMP)

# For OpenMP
if (OPENMP_FOUND)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)


set (CMAKE_CXX_STANDARD 11)


# Use sse2 for math
IF (UNIX)
   set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2 -mfpmath=sse" )
   set( PLATFORM_LIBS m )
ENDIF (UNIX)


IF(MSVC)
	ADD_DEFINITIONS(/arch:SSE2)
ENDIF(MSVC)


INCLUDE_DIRECTORIES( ${PROJECT_BINARY_DIR} )

add_library ( All ${lp_field_SRCS} )

add_executable( SmoothDist src/SmoothDist.cpp )
target_link_libraries( SmoothDist All ${PLATFORM_LIBS} )

add_executable( Psi1 src/Psi1.cpp )
target_link_libraries( Psi1 All ${PLATFORM_LIBS} )

add_executable( SmoothDistNumerical src/SmoothDistNumerical.cpp )
target_link_libraries( SmoothDistNumerical All ${PLATFORM_LIBS} )

add_executable( SmoothDistHybrid src/SmoothDistHybrid.cpp )
target_link_libraries( SmoothDistHybrid All ${PLATFORM_LIBS} )

add_executable( SmoothDistNumericalFast src/SmoothDistNumericalFast.cpp )
target_link_libraries( SmoothDistNumericalFast All ${PLATFORM_LIBS} )

add_executable( Test_SmoothDist src/Test_SmoothDist.cpp )
target_link_libraries( Test_SmoothDist All ${PLATFORM_LIBS} )
