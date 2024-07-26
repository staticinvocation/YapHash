cmake_minimum_required(VERSION 3.1...3.26)
project(YapHash VERSION 1.0 LANGUAGES CXX)
find_package(ipp REQUIRED)
find_package(MKL CONFIG REQUIRED)

add_library(YapHash src/Array2d_32f.cpp src/Audio.cpp src/melFB.cpp src/parseConfig.cpp src/Stft.cpp
#src/Stwt.cpp
src/VAD.cpp src/VIATUtilities.cpp src/Wave.cpp 
#src/WaveLash.cpp
src/YapHash.cpp)
target_link_libraries(YapHash PUBLIC ${IPP_LIBRARIES} MKL::MKL)
target_compile_definitions(YapHash PRIVATE _CRT_SECURE_NO_WARNINGS)