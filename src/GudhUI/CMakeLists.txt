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

if(MSVC)
  # Turn off some VC++ warnings
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11 -Wall -Wpedantic -Wsign-compare")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb -O0")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endif()

message("CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}")
message("CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE}")
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

else()
  message(STATUS "NOTICE: GudhUI requires CGAL, the QGLViewer, OpenGL and Qt4, and will not be compiled.")
endif()
