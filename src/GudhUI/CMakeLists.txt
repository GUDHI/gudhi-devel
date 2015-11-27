cmake_minimum_required(VERSION 2.8)
project(GudhUI)

find_package(CGAL  COMPONENTS Qt4)
find_package(Qt4)
find_package(QGLViewer)
find_package(OpenGL)

if ( CGAL_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND )
  set( QT_USE_QTXML    TRUE )
  set( QT_USE_QTMAIN   TRUE )
  set( QT_USE_QTSCRIPT  TRUE )
  set( QT_USE_QTOPENGL  TRUE )
  SET(Boost_USE_STATIC_LIBS ON)
  SET(Boost_USE_MULTITHREAD OFF)

  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
  
  include(${QT_USE_FILE})

  include(${CGAL_USE_FILE})
  # In CMakeLists.txt, when include(${CGAL_USE_FILE}), CXX_FLAGS are overwritten.
  # cf. http://doc.cgal.org/latest/Manual/installation.html#title40
  # A workaround is to add "-std=c++11" again.
  # A fix would be to use https://cmake.org/cmake/help/v3.1/prop_gbl/CMAKE_CXX_KNOWN_FEATURES.html 
  # or even better https://cmake.org/cmake/help/v3.1/variable/CMAKE_CXX_STANDARD.html
  # but it implies to use cmake version 3.1 at least.
  if(NOT MSVC)
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG(-std=c++11 COMPILER_SUPPORTS_CXX11)
    if(COMPILER_SUPPORTS_CXX11)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()
  endif()
  # - End of workaround

  include_directories (${QGLVIEWER_INCLUDE_DIR})
  include_directories(.)	
 
  # qt : ui file, created wih Qt Designer ###############
  set(CMAKE_CURRENT_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gui")
  qt4_wrap_ui( uis 
    gui/main_window.ui
    gui/MenuEdgeContraction.ui
    gui/KNearestNeighborsMenu.ui
    gui/UniformNeighborsMenu.ui
    gui/PersistenceMenu.ui
  )	

  set(CMAKE_CURRENT_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gui")
  qt4_automoc(
    gui/MainWindow.cpp
    gui/Menu_k_nearest_neighbors.cpp
    gui/Menu_uniform_neighbors.cpp
    gui/Menu_edge_contraction.cpp
    gui/Menu_persistence.cpp    
  )
	
  set(CMAKE_CURRENT_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/view") 
  qt4_automoc(view/Viewer_instructor.cpp		
    view/Viewer.cpp
  )	
  #####################################################################

  add_executable  ( GudhUI 
    gui/gudhui.cpp
    gui/MainWindow.cpp
    gui/Menu_k_nearest_neighbors.cpp	
    gui/Menu_uniform_neighbors.cpp
    gui/Menu_edge_contraction.cpp  
    gui/Menu_persistence.cpp
    view/Viewer_instructor.cpp				
    view/Viewer.cpp
    ${uis}  
  )	
	
  target_link_libraries( GudhUI ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} )
  target_link_libraries( GudhUI ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )	

###############################################################################
  if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    include( ${CGAL_USE_FILE} )
     # In CMakeLists.txt, when include(${CGAL_USE_FILE}), CXX_FLAGS are overwritten.
    # cf. http://doc.cgal.org/latest/Manual/installation.html#title40
    # A workaround is to add "-std=c++11" again.
    # A fix would be to use https://cmake.org/cmake/help/v3.1/prop_gbl/CMAKE_CXX_KNOWN_FEATURES.html 
    # or even better https://cmake.org/cmake/help/v3.1/variable/CMAKE_CXX_STANDARD.html
    # but it implies to use cmake version 3.1 at least.
    if(NOT MSVC)
      include(CheckCXXCompilerFlag)
      CHECK_CXX_COMPILER_FLAG(-std=c++11 COMPILER_SUPPORTS_CXX11)
      if(COMPILER_SUPPORTS_CXX11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
      endif()
    endif()
    # - End of workaround

    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      
      add_executable (acp alpha_complex_persistence.cpp)
      target_link_libraries(acp ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY} ${QT_LIBRARIES})

    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
###############################################################################

else()
  message(STATUS "NOTICE: GudhUI requires CGAL, the QGLViewer, OpenGL and Qt4, and will not be compiled.")
endif()
