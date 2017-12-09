# Copy files, otherwise result files are created in sources
file(COPY "${CMAKE_SOURCE_DIR}/data/persistence_diagram/first.pers" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
file(COPY "${CMAKE_SOURCE_DIR}/data/persistence_diagram/second.pers" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

function(add_persistence_representation_creation_utility creation_utility)
    add_executable ( ${creation_utility} ${creation_utility}.cpp )

    # as the function is called in a subdirectory level, need to '../' to find persistence files
    # ARGN will add all the other arguments (except creation_utility) sent to the CMake functions
    add_test(NAME Persistence_representation_utilities_${creation_utility} COMMAND $<TARGET_FILE:${creation_utility}>
        ${ARGN} "${CMAKE_CURRENT_BINARY_DIR}/../first.pers"
        "${CMAKE_CURRENT_BINARY_DIR}/../second.pers")
endfunction(add_persistence_representation_creation_utility)

function(add_persistence_representation_plot_utility plot_utility tool_extension)
    add_executable ( ${plot_utility} ${plot_utility}.cpp )

    # as the function is called in a subdirectory level, need to '../' to find persistence heat maps files
    add_test(NAME Persistence_representation_utilities_${plot_utility}_first COMMAND $<TARGET_FILE:${plot_utility}>
        "${CMAKE_CURRENT_BINARY_DIR}/../first.pers${tool_extension}")
    add_test(NAME Persistence_representation_utilities_${plot_utility}_second COMMAND $<TARGET_FILE:${plot_utility}>
        "${CMAKE_CURRENT_BINARY_DIR}/../second.pers${tool_extension}")
    if(GNUPLOT_PATH)
        add_test(NAME Persistence_representation_utilities_${plot_utility}_first_gnuplot COMMAND ${GNUPLOT_PATH}
            "-e" "load '${CMAKE_CURRENT_BINARY_DIR}/../first.pers${tool_extension}_GnuplotScript'")
        add_test(NAME Persistence_representation_utilities_${plot_utility}_second_gnuplot COMMAND ${GNUPLOT_PATH}
            "-e" "load '${CMAKE_CURRENT_BINARY_DIR}/../second.pers${tool_extension}_GnuplotScript'")
    endif()
endfunction(add_persistence_representation_plot_utility)

function(add_persistence_representation_function_utility function_utility tool_extension)
    add_executable ( ${function_utility} ${function_utility}.cpp )

    # ARGV2 is an optional argument
    if (${ARGV2})
        # as the function is called in a subdirectory level, need to '../' to find persistence heat maps files
        add_test(NAME Persistence_representation_utilities_${function_utility} COMMAND $<TARGET_FILE:${function_utility}>
            "${ARGV2}"
            "${CMAKE_CURRENT_BINARY_DIR}/../first.pers${tool_extension}"
            "${CMAKE_CURRENT_BINARY_DIR}/../second.pers${tool_extension}")
    else()
        # as the function is called in a subdirectory level, need to '../' to find persistence heat maps files
        add_test(NAME Persistence_representation_utilities_${function_utility} COMMAND $<TARGET_FILE:${function_utility}>
            "${CMAKE_CURRENT_BINARY_DIR}/../first.pers${tool_extension}"
            "${CMAKE_CURRENT_BINARY_DIR}/../second.pers${tool_extension}")
    endif()
endfunction(add_persistence_representation_function_utility)

add_subdirectory(persistence_heat_maps)
add_subdirectory(persistence_intervals)
add_subdirectory(persistence_landscapes)
add_subdirectory(persistence_landscapes_on_grid)
add_subdirectory(persistence_vectors)

