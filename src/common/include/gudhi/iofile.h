//todo remove and use off_reader instead

#ifndef GUDHI_IOFILE_H_
#define GUDHI_IOFILE_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <list>
#include <cassert>



//todo use my new visitor based off instead
/**
 * @brief OFF reader, save the content of the file (vertices and maximal faces) in 'complex'.
 * The class Complex has to handle the following operations:
 * - void add_vertex(double[dim],int dim)
 * - void add_face(int dimension ,int[dimension] vertices)
 * Source from : http://www.holmes3d.net/graphics/offfiles/OFFLoading.txt
 *
 * @todo todo ignore comments in the file -> # comment
 */
template<typename Complex> inline
bool general_read_off_file(const std::string & file_name, Complex& complex){
  // Declare temporary variables to read data into.
  // If the read goes well, we'll copy these into
  // our class variables, overwriting what used to
  // be there. If it doesn't, we won't have messed up
  // our previous data structures.
  int tempNumPoints   = 0;  // Number of x,y,z coordinate triples
  int tempNumFaces    = 0;  // Number of polygon sets
  int tempNumEdges    = 0;  // Unused, except for reading.
  int tempDimPoints   = 0;
  double** tempPoints = NULL;  // An array of x,y,z coordinates.
  int** tempFaces = NULL;    // An array of arrays of point
  // pointers. Each entry in this
  // is an array of integers. Each
  // integer in that array is the
  // index of the x, y, and z
  // coordinates in the corresponding
  // arrays.
  int*  tempFaceSizes = NULL;  // An array of polygon point counts.
  // Each of the arrays in the tempFaces
  // array may be of different lengths.
  // This array corresponds to that
  // array, and gives their lengths.
  int i;        // Generic loop variable.
  bool goodLoad = true;    // Set to false if the file appears
  // not to be a valid OFF file.
  char tempBuf[128];    // A buffer for reading strings
  // from the file.

  // Create an input file stream for the file the CArchive
  // is connected to. This allows use of the overloaded
  // extraction operator for conversion from string
  // data to doubles and ints.
  std::ifstream ifs (file_name.c_str(), std::ios::in);
  if(!ifs.is_open()) {
    std::cerr << "Unable to open file " << file_name << std::endl;
    return false;
  }

  // Grab the first string. If it's "OFF", we think this
  // is an OFF file and continue. Otherwise we give up.
  ifs >> tempBuf;
  if (strcmp(tempBuf, "OFF") != 0) {
    goodLoad = false;
    std::cerr << "No OFF preambule\n";
  }

  // Read the sizes for our two arrays, and the third
  // int on the line. If the important two are zero
  // sized, this is a messed up OFF file. Otherwise,
  // we setup our temporary arrays.
  if (goodLoad) {
    ifs >> tempNumPoints >> tempNumFaces >> tempNumEdges;
    if (tempNumPoints < 1 || tempNumFaces < 0) {
      // If either of these were negative, we make
      // sure that both are set to zero. This is
      // important for later deleting our temporary
      // storage.
      goodLoad      = false;
      std::cerr << "tempNumPoints < 1 || tempNumFaces < 0\n";
      tempNumPoints = 0;
      tempNumFaces  = 0;
    } else {
      tempPoints = new double*[tempNumPoints];
      tempFaces = new int*[tempNumFaces];
      tempFaceSizes = new int[tempNumFaces];
    }
  }

  if (goodLoad) {
    // Load all of the points.

    // we start by loading the first one
    // the case is difference because then we dont know the dimension by advance
    // we discover the point dimension by reading the first line
    std::string lineFirstPoint;
	std::getline(ifs, lineFirstPoint);
    while(lineFirstPoint.size()<3){
    	std::getline(ifs, lineFirstPoint);
    }

    // we store the first point in a temporary list
    std::istringstream lineFirstPointStream(lineFirstPoint);
    std::list<double> firstTempPoint;
    double coord;
    while(lineFirstPointStream>>coord){
      firstTempPoint.push_back(coord);
      ++tempDimPoints;
    }
    // we store the point in our points array
    tempPoints[0]=new double[tempDimPoints];
    for( int j = 0 ; j<tempDimPoints; ++j){
      tempPoints[0][j] = firstTempPoint.front();
      firstTempPoint.pop_front();
    }

    // now, we know the dimension and can read safely other points
    for (i = 1; i < tempNumPoints; i++) {
      tempPoints[i] = new double[tempDimPoints];
      for (int j = 0 ; j<tempDimPoints ; ++j){
        ifs>>tempPoints[i][j];
      }
    }

    // Load all of the faces.
    for (i = 0; i < tempNumFaces; i++) {
      // This tells us how many points make up
      // this face.
      ifs >> tempFaceSizes[i];
      // So we declare a new array of that size
      tempFaces[i] = new int[tempFaceSizes[i]];
      // And load its elements with the vertex indices.
      for (int j = 0; j < tempFaceSizes[i]; j++) {
        ifs >> tempFaces[i][j];
      }
      // Clear out any face color data by reading up to
      // the newline. 128 is probably considerably more
      // space than necessary, but better safe than
      // sorry.
      ifs.getline(tempBuf, 128);
    }
  }

  // Here is where we copy the data from the temp
  // structures into our permanent structures. We
  // probably will do some more processing on the
  // data at the same time. This code you must fill
  // in on your own.
  if (goodLoad) {
    // we save vertices first in the complex
    for (i = 0; i < tempNumPoints; i++)
      complex.add_vertex(tempPoints[i],tempDimPoints);

    // we save faces
    for (i = 0; i < tempNumFaces; i++) {
      for (int j = 0; j < tempFaceSizes[i]; j++)
        complex.add_face(tempFaceSizes[i],tempFaces[i]);
    }
  }

  // Now that we're done, we have to make sure we
  // free our dynamic memory.
  for (i = 0; i < tempNumPoints; i++) {
    delete []tempPoints[i];
  }
  delete []tempPoints;

  for (i = 0; i < tempNumFaces; i++) {
    delete tempFaces[i];
  }
  delete []tempFaces;
  delete []tempFaceSizes;

  // Clean up our ifstream. The MFC framework will
  // take care of the CArchive.
  ifs.close();

  return goodLoad;
}


template<typename Complex>
class Geometric_flag_complex_wrapper{
  Complex& complex_;
  typedef typename Complex::Vertex_handle Vertex_handle;
  typedef typename Complex::Point Point;

  const bool load_only_points_;

public:
  Geometric_flag_complex_wrapper(Complex& complex,bool load_only_points = false):
    complex_(complex),
    load_only_points_(load_only_points)
{}


  void add_vertex(double* xyz,int dim){
    Point p(dim);
    for(int i=0;i<dim;++i)
      p[i] = xyz[i];
    complex_.add_vertex(p);
  }

  void add_face(int dimension ,int* vertices){
    if (!load_only_points_){
      for (int i = 0; i<dimension ; ++i)
        for (int j = i+1; j<dimension ; ++j)
          complex_.add_edge(Vertex_handle(vertices[i]),Vertex_handle(vertices[j]));
    }
  }
};




/**
 * @brief Read a mesh into a OFF file
 * load_only_points should be true if only the points have to be loaded.
 */
template<typename Complex>
bool read_off_file(std::string file_name,Complex &complex,bool load_only_points = false){
  complex.clear();
  Geometric_flag_complex_wrapper<Complex> complex_wrapper(complex,load_only_points);
  return general_read_off_file(file_name,complex_wrapper);

}


#endif
