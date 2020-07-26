CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

PROJECT( lp_field )

# Platform specific libraries
IF ( UNIX )
   SET ( PLATFORM_LIBS m )
ENDIF ()


INCLUDE_DIRECTORIES( ${PROJECT_BINARY_DIR} )


# ADD_LIBRARY( AllLib m )


ADD_EXECUTABLE( lp_field src/lp_field.cpp )
TARGET_LINK_LIBRARIES ( lp_field ${PLATFORM_LIBS} )

