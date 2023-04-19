#ifndef EDGE_H
#define EDGE_H
#include <stddef.h>
#include "vertex.h"
using namespace std;
// This is the Edge class representing an edge.

class vertex;

/**

\brief The Edge class representing an edge.

This class represents an edge in a graph with two vertices, origin and
destination, and two adjacent edges, prev and next.

It also has a twin edge, which represents the same edge in the opposite
direction.
*/
class Edge {
public:
  /// Pointer to the vertex at the start of the edge.
  Vertex *origin = nullptr;

  /// Pointer to the twin edge, representing the same edge in the opposite
  /// direction.
  Edge *twin = nullptr;

  /// Pointer to the previous edge.
  Edge *prev = nullptr;

  /// Pointer to the next edge.
  Edge *next = nullptr;

  /// Pointer to the vertex at the end of the edge.
  Vertex *dest = nullptr;
};


#endif