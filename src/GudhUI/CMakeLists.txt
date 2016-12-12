cmake_minimum_required(VERSION 2.8)
project(GudhUI)

find_package(CGAL COMPONENTS Qt5)
find_package(Qt5 COMPONENTS Widgets Xml OpenGL)
find_package(QGLViewer)
find_package(OpenGL)


if (NOT CGAL_FOUND)
  message(ERROR " GudhUI requires CGAL and will not be compiled.")
endif()

if (NOT Qt5_FOUND)
  message(ERROR " GudhUI requires Qt5 and will not be compiled.")
endif()

if (NOT OPENGL_FOUND)
  message(ERROR " GudhUI requires OpenGL and will not be compiled.")
endif()

if (NOT QGLVIEWER_FOUND)
  message(ERROR " GudhUI requires QGLViewer and will not be compiled.")
endif()



if ( CGAL_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND )

  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)

  SET(Boost_USE_STATIC_LIBS ON)
  SET(Boost_USE_MULTITHREAD OFF)
  include_directories (${QGLVIEWER_INCLUDE_DIR})
 
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
  )	
  target_link_libraries( GudhUI Qt5::Widgets Qt5::Xml Qt5::OpenGL )
  target_link_libraries( GudhUI ${QGLVIEWER_LIBRARIES} )
  target_link_libraries( GudhUI ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )	
if (TBB_FOUND)
  target_link_libraries( GudhUI ${TBB_LIBRARIES})
endif()

###############################################################################

else()
  message(STATUS "NOTICE: GudhUI requires CGAL, the QGLViewer, OpenGL and Qt5, and will not be compiled.")
endif()
