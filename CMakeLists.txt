cmake_minimum_required(VERSION 3.5.0)
project(cbne VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

 # either set it here or from the command line  
set(BOOST_ROOT "/home/abhayat/boost_libraries/")
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost REQUIRED COMPONENTS graph)

# change FALSE to TRUE in order to build and run tests
if(FALSE)

  find_package(Catch2 3 REQUIRED)

  add_executable(tests tests/test.cpp Complex.h Complex.cpp Config.h Statistics.h)
  target_link_libraries(tests PRIVATE ${Boost_LIBRARIES} Catch2::Catch2WithMain)

  add_executable(benchmarks tests/benchmark.cpp Complex.h Complex.cpp Config.h Statistics.h)
  target_link_libraries(benchmarks PRIVATE ${Boost_LIBRARIES} Catch2::Catch2WithMain)

  include(Catch)
  catch_discover_tests(tests)
  catch_discover_tests(benchmarks)

endif()


add_executable(cbne Cbne.h Complex.h Complex.cpp main.cpp Config.h Statistics.h)
target_include_directories(cbne PUBLIC ${Boost_INCLUDE_DIRS}) 
target_link_libraries(cbne ${Boost_LIBRARIES})
target_include_directories(cbne PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
