# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

# This project is designed to be built outside the Insight source tree.
project(USshoulder)

# Tell Cmake to find and process ITK
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

# Tell Cmake to find and process SimpleITK
find_package(SimpleITK REQUIRED)
include(${SimpleITK_USE_FILE})

ENABLE_TESTING()
INCLUDE(Dart)

INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

SET(Required_Packages 
ITK)
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

ADD_EXECUTABLE(USshoulder 
USshoulder_codeForTesting.cpp 
)
							 
TARGET_LINK_LIBRARIES(USshoulder ${ITK_LIBRARIES} ${SimpleITK_LIBRARIES} ${Libraries})
