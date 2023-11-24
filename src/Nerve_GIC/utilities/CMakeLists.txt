project(Nerve_GIC_examples)

add_executable ( Nerve Nerve.cpp )
add_executable ( VoronoiGIC VoronoiGIC.cpp )

if(TARGET TBB::tbb)
  target_link_libraries(Nerve TBB::tbb)
  target_link_libraries(VoronoiGIC TBB::tbb)
endif()

file(COPY KeplerMapperVisuFromTxtFile.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
# Copy files for not to pollute sources when testing
file(COPY "${CMAKE_SOURCE_DIR}/data/points/human.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_test(NAME Nerve_GIC_utilities_nerve COMMAND $<TARGET_FILE:Nerve>
    "human.off" "2" "10" "0.3")

add_test(NAME Nerve_GIC_utilities_VoronoiGIC COMMAND $<TARGET_FILE:VoronoiGIC>
    "human.off" "100")

install(TARGETS Nerve DESTINATION bin)
install(TARGETS VoronoiGIC DESTINATION bin)
install(FILES KeplerMapperVisuFromTxtFile.py DESTINATION bin)
