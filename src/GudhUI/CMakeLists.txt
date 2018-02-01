cmake_minimum_required(VERSION 2.8)
project(GudhUI)

# Need to find OpenGL first as find_package(Qt5) tries to #include"GL/gl.h" on some platforms
find_package(OpenGL)

if (OPENGL_FOUND)
  find_package(Qt5 COMPONENTS Widgets Xml OpenGL)
  find_package(QGLViewer)
  
  if ( CGAL_FOUND AND Qt5_FOUND AND QGLVIEWER_FOUND AND NOT CGAL_VERSION VERSION_EQUAL 4.8.0)
  
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
  
    SET(Boost_USE_STATIC_LIBS ON)
    SET(Boost_USE_MULTITHREAD OFF)
    include_directories (${QGLVIEWER_INCLUDE_DIR})
   
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
  
    install(TARGETS GudhUI DESTINATION bin)
  
  endif()
endif(OPENGL_FOUND)