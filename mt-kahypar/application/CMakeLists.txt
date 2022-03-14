message(STATUS "app CS:" ${CMAKE_CURRENT_SOURCE_DIR})
message(STATUS "app CL:" ${CMAKE_CURRENT_LIST_DIR})

add_library(MtKaHyParStrongAsync INTERFACE )
target_include_directories(MtKaHyParStrongAsync INTERFACE ../../external_tools/kahypar)
target_include_directories(MtKaHyParStrongAsync INTERFACE ../../)

target_sources(MtKaHyParStrongAsync INTERFACE kahypar.cc)
target_link_libraries(MtKaHyParStrongAsync INTERFACE  ${Boost_LIBRARIES})
target_link_libraries(MtKaHyParStrongAsync INTERFACE  pthread tbb)
#set_property(TARGET MtKaHyParStrongAsync PROPERTY CXX_STANDARD 17)
#set_property(TARGET MtKaHyParStrongAsync PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_definitions(MtKaHyParStrongAsync INTERFACE  USE_STRONG_PARTITIONER)
target_compile_definitions(MtKaHyParStrongAsync INTERFACE  USE_ASYNC_UNCOARSENING)

if(ENABLE_PROFILE MATCHES ON)
  target_link_libraries(MtKaHyParStrongAsync ${PROFILE_FLAGS})
endif()

add_executable(MtKaHyParMQRandomC4)
target_link_libraries(MtKaHyParMQRandomC4 PRIVATE MtKaHyParStrongAsync mq_random benchmark_base)
target_compile_definitions(MtKaHyParMQRandomC4 PRIVATE MQ_C=4)
target_compile_features(knapsack_heap INTERFACE cxx_std_17)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} MtKaHyParMQRandomC4 PARENT_SCOPE)
