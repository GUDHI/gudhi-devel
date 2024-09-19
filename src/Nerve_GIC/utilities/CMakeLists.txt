add_executable_with_targets(Nerve Nerve.cpp TBB::tbb)
add_executable_with_targets(VoronoiGIC VoronoiGIC.cpp TBB::tbb)

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
