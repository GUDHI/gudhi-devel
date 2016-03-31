cmake_minimum_required(VERSION 2.8)
project(GudhUI)

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

  include(${QT_USE_FILE})

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
if (TBB_FOUND)
  target_link_libraries( GudhUI ${TBB_RELEASE_LIBRARY})
endif()

else()
  message(STATUS "NOTICE: GudhUI requires CGAL, the QGLViewer, OpenGL and Qt4, and will not be compiled.")
endif()
