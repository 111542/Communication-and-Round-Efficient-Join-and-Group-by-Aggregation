find_package(ABY QUIET)
if(ABY_FOUND)
   message(STATUS "Found ABY")
elseif (NOT ABY_FOUND AND NOT TARGET ABY::aby)
   message("ABY was not found: add ABY subdirectory")
   add_subdirectory(ABY)
endif()

# Adding MPI Library
find_package(MPI REQUIRED)
if (MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
else (MPI_FOUND)
    message(SEND_ERROR "This application cannot compile without MPI")
endif (MPI_FOUND)

#find_package(LibNL REQUIRED)
find_package(Threads REQUIRED)

# 包含自定义的查找模块
#include(cmake/FindSodium.cmake)


find_package(PkgConfig REQUIRED)
pkg_check_modules(_LIBSODIUM REQUIRED libsodium)

find_path(SODIUM_INCLUDE_DIR sodium.h HINTS ${LIBSODIUM_INCLUDE_DIRS} /usr/local/include /opt/local/include /opt/include)
find_library(SODIUM_LIBRARY NAMES sodium HINTS ${LIBSODIUM_LIBRARY_DIRS} /usr/local/lib /opt/local/lib /opt/lib)

add_executable(test_bool test_bool.cpp)


target_link_libraries(test_bool PUBLIC ABY::aby ${MPI_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${SODIUM_LIBRARY})



target_include_directories(test_bool PUBLIC ${SODIUM_INCLUDE_DIR})


