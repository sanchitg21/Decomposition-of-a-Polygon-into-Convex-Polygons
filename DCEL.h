#ifndef DCEL_H
#define DCEL_H

#include "Face.h"
#include "Edge.h"
#include "vertex.h"
#include<map>

class Face;
class vertex;
class Edge;

// This is the DCEL class representing a doubly connected edge list.
/**

\brief The DCEL class representing a doubly connected edge list.

This class represents a doubly connected edge list (DCEL) data structure, which
is used to represent planar subdivisions of a polygon.

The DCEL consists of a set of vertices, edges and faces. The vertices and edges
have pointers to each other to form a connected structure.

The faces are bounded by a set of edges.
*/
class DCEL {
public:
  ///\brief Vector of pointers to the vertices in the DCEL.
  vector<Vertex *> Vert;

  ///\brief Vector of pointers to the edges in the DCEL.
  vector<Edge *> Edges;

  ///\brief Vector of pointers to the faces in the DCEL.
  vector<Face *> faces;

  ///\brief Map of pointers to the vertices in the DCEL to the edge that starts
  /// at that vertex and the next vertex in the polygon.
  map<Vertex *, pair<Edge *, Vertex *>> mp;

  ///\brief Constructor of DCEL class.
  ///
  ///\param vertices Vector of pointers to the vertices in the polygon.

  // Constructor of DCEL class
    DCEL(vector<Vertex *> &vertices);
    DCEL* build_polygon(vector<Vertex *> &vertices);
};


#endif