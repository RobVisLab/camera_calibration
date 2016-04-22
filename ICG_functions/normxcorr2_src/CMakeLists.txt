project (cudaApplications CXX C)

cmake_minimum_required(VERSION 2.6)

FILE(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/" OT_CMAKE_MODULE_PATH)
SET(CMAKE_MODULE_PATH ${OT_CMAKE_MODULE_PATH})
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../private")

# COMPILER FLAGS
IF(WIN32)
  ADD_DEFINITIONS(-DWIN32)
ELSE(WIN32)
  ADD_DEFINITIONS( -Wall)
  ADD_DEFINITIONS(-O3 -march=core2 )
ENDIF(WIN32)
#SET(CMAKE_OSX_ARCHITECTURES i386)

##-----------------------------------------------------------------------------
## OpenCV
FIND_PACKAGE(OpenCV REQUIRED)
INCLUDE_DIRECTORIES( ${OpenCV_INCLUDE_DIRS} )

##-----------------------------------------------------------------------------
## MATLAB
FIND_PACKAGE(Matlab_ICG)
INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIR})


SET( SOURCES
  normxcorr2_mex.cpp
  )
  
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  )
  
ADD_LIBRARY(normxcorr2_mex SHARED ${HEADERS} ${SOURCES})

SET_TARGET_PROPERTIES(normxcorr2_mex PROPERTIES PREFIX "")
SET_TARGET_PROPERTIES(normxcorr2_mex PROPERTIES SUFFIX ${MATLAB_SUFFIX})
if(WIN32)
	SET_TARGET_PROPERTIES(normxcorr2_mex PROPERTIES LINK_FLAGS "/export:mexFunction")
endif(WIN32)
if(WIN32)
target_link_libraries(normxcorr2_mex ${OpenCV_LIBRARIES} ${MATLAB_LIBRARIES})
else(WIN32)
target_link_libraries(normxcorr2_mex ${OpenCV_LIBRARIES} ${MATLAB_LIBRARIES} m)
endif(WIN32)
