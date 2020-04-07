cmake_minimum_required(VERSION 3.0)
project(gLike VERSION 0.8.0)

# load ROOT libraries
set(ROOT_REQUIRED_COMPONENTS Minuit Fitsio)
find_package(ROOT REQUIRED ${ROOT_REQUIRED_COMPONENTS})

# include dirs
include_directories(${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include)

# build one shared lib with all sources
add_library(glike
    ${PROJECT_SOURCE_DIR}/src/Lkl.cc
    ${PROJECT_SOURCE_DIR}/src/PoissonLkl.cc
)
target_link_libraries(glike ${ROOT_LIBRARIES})

add_executable(testPoissonLkl ${PROJECT_SOURCE_DIR}/scripts/testPoissonLkl.cc)
target_link_libraries(testPoissonLkl glike ${ROOT_LIBRARIES})
