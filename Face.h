#ifndef FACE_H
#define FACE_H

#include "vertex.h"
#include <vector>
using namespace std;
// This is the Face class representing a face.

class vertex;

/**

\brief The Face class representing a face.
This class represents a face in a planar graph with a boundary defined by a
sequence of vertices in clockwise order.
*/
class Face {
public:
  /// Vector of pointers to the vertices that form the boundary of the face in
  /// clockwise order.
  vector<Vertex *> vertices = {};
};


#endif