# CMakeList.txt : CMake project for MaterialSimulator, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("MaterialSimulator")

# Add source to this project's executable.
find_package(FFTW3 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
include_directories(C:/dev/spectra/include)



add_executable (MaterialSimulator "MaterialSimulator.cpp" "MaterialSimulator.h"  "Fields.h" "Fields.cpp" "LinearOperator.h"  "DFT.h" "DFT.cpp" "DebuggingWindow.cpp" "DebuggingWindow.h")
if(FFTW3_FOUND)
    target_link_libraries(MaterialSimulator PRIVATE FFTW3::fftw3)
    target_link_libraries(MaterialSimulator PRIVATE Eigen3::Eigen)
endif() 

if(OpenMP_CXX_FOUND)
    target_link_libraries(MaterialSimulator PRIVATE OpenMP::OpenMP_CXX)
endif()

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET MaterialSimulator PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.



