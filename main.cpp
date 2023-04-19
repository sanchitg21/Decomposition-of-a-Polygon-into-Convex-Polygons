#include <bits/stdc++.h>
#include <chrono>
using namespace std;

#include "DCEL.h"
#include "Face.h"
#include "Edge.h"
#include "vertex.h"
class Vertex;
class Edge;
class Face;

DCEL::DCEL(vector<Vertex *> &vertices)
{
  int numVertices = vertices.size();
  Vert = vertices;
  for (int i = 0; i < numVertices; i++)
  {
    Edge *edge = new Edge();
    edge->origin = vertices[i];
    edge->dest = vertices[(i + 1) % numVertices];
    Edge *twin1 = new Edge();
    twin1->origin = vertices[(i + 1) % numVertices];
    twin1->dest = vertices[i];
    twin1->twin = edge;
    edge->twin = twin1;
    Edges.push_back(edge);
  }

  for (int i = 0; i < numVertices; i++)
  {
    Edges[i]->next = Edges[(i + 1) % numVertices];
    Edges[i]->prev = Edges[(i - 1 + numVertices) % numVertices];
    Edges[i]->twin->next = Edges[(i - 1 + numVertices) % numVertices]->twin;
    Edges[i]->twin->prev = Edges[(i + 1) % numVertices]->twin;
  }
  for (int i = 0; i < numVertices; i++)
  {
    mp[vertices[i]] = {Edges[i], vertices[(i + 1) % numVertices]};
  }
}

/**

@brief Vector of all notches in the polygon.
This vector contains pointers to all the notches in the polygon.
Each notch is represented by a Vertex object pointer.
*/
vector<Vertex *> total_notches;
/**

@brief Vector of all diagonals in the polygon.
This vector contains pairs of pointers to the vertices that make up each
diagonal in the polygon. Each pair is represented by a std::pair object with two
Vertex object pointers.
*/
vector<pair<Vertex *, Vertex *>> diagonals;

/**

    @brief Calculates the vector of all notches in the given polygon.
    This function calculates and returns a vector of pointers to all notches in the
    given polygon. A notch is defined as a vertex that makes a concave angle in the
    polygon.
    @param polygon A pointer to a DCEL object representing the polygon.
    @return vector<Vertex*> A vector of pointers to all notches in the polygon.
*/

vector<Vertex *> calc_notches(DCEL *polygon)
{
  /// Calculate the number of vertices in the polygon
  int n = polygon->Vert.size();
  /// Initialize the vector of notches
  vector<Vertex *> ans;

  /// Iterate over all vertices in the polygon
  for (int i = 0; i < n; i++)
  {
    /// Get pointers to the three adjacent vertices
    Vertex *v1 = polygon->Vert[(i - 1 + n) % n];
    Vertex *v2 = polygon->Vert[i];
    Vertex *v3 = polygon->Vert[(i + 1) % n];

    /// Check if v2 is a notch by checking if the cross product of the vectors
    /// (v2 - v1) and (v3 - v2) is positive
    if ((v2->x - v1->x) * (v3->y - v2->y) - (v3->x - v2->x) * (v2->y - v1->y) >
        0)
    {
      ans.push_back(v2);
    }
  }

  // Return the vector of notches
  return ans;
}

/**
    @brief Check if the given point lies inside the given polygon.
    *Implementing a point-in-polygon test using the "crossing number"
    algorithm. Given a list of polygon vertices, the function checks if a point with
    coordinates (x,y) lies inside the polygon. The algorithm determines if a point
    is inside a polygon by checking the number of times a horizontal line passing
    through the point intersects with the polygon's edges. If the number of
    intersections is odd, the point is inside the polygon; otherwise, it is outside.
    @param vertices A vector of pointers to the vertices of the polygon.
    @param testx The x-coordinate of the point to test.
    @param testy The y-coordinate of the point to test.
    @return An integer value of 1 if the point is inside the polygon, and 0 if it is
    outside.
*/
int pnpoly(vector<Vertex *> vertices, double testx, double testy)
{
  int vert_num = vertices.size();
  vertices.push_back(vertices[0]);
  vector<double> vertx;
  vector<double> verty;
  for (int i = 0; i < vertices.size(); i++)
  {
    vertx.push_back(vertices[i]->x);
    verty.push_back(vertices[i]->y);
  }
  int i, j, c = 0;
  for (i = 0, j = vert_num - 1; i < vert_num; j = i++)
  {
    if (((verty[i] > testy) != (verty[j] > testy)) &&
        (testx <
         (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) +
             vertx[i]))
      c = !c;
  }
  return c;
}

/**

    @brief Checks whether the given three vertices form a notch in clockwise order.
    @param v1 Pointer to the first vertex.
    @param v2 Pointer to the second vertex.
    @param v3 Pointer to the third vertex.
    @return Returns true if the vertices form a notch in clockwise order, otherwise
    returns false.
*/
bool is_notch(Vertex *v1, Vertex *v2, Vertex *v3)
{
  if ((((v2->x - v1->x) * (v3->y - v2->y)) -
       ((v3->x - v2->x) * (v2->y - v1->y))) <= 0)
  {
    return false;
  }
  else
  {
    return true;
  }
}

/**

    @brief Builds a doubly connected edge list (DCEL) representing the polygon
    formed by the given vertices.
    @param vertices The vector of pointers to the vertices of the polygon.
    @return A pointer to the DCEL representing the polygon.

*/
DCEL *DCEL::build_polygon(vector<Vertex *> &vertices)
{
  DCEL *polygon = new DCEL(vertices);
  return polygon;
}

/**

    @brief Creates a minimal rectangle from the given vertices and returns the
    boundary of the rectangle.
    @param L A vector of pointers to the vertices of a polygon.
    @return A vector containing the coordinates of the top right and bottom left
    corners of the minimal rectangle. The first two elements represent the x
    coordinates of the top right and bottom left corners, respectively, while the
    last two elements represent the y coordinates of the top right and bottom left
    corners, respectively.
*/
vector<double> calc_rectangle(vector<Vertex *> L)
{
  vector<double> rec(4);
  rec[0] = INT_MIN;
  rec[1] = INT_MAX;
  rec[2] = INT_MIN;
  rec[3] = INT_MAX;
  for (int i = 0; i < L.size(); i++)
  {
    if (L[i]->x > rec[0])
    {
      rec[0] = L[i]->x;
    }
    if (L[i]->x < rec[1])
    {
      rec[1] = L[i]->x;
    }
    if (L[i]->y > rec[2])
    {
      rec[2] = L[i]->y;
    }
    if (L[i]->y < rec[3])
    {
      rec[3] = L[i]->y;
    }
  }
  return rec;
}

/**

    @brief This function checks whether the given current vertex lies to the right
    of the line formed by the given notch and vertex v1.
    @param v1 Pointer to the first vertex.
    @param notch Pointer to the second vertex which forms the notch with v1.
    @param curr Pointer to the current vertex to check if it lies to the right of
    the line.
    @return true if the curr vertex lies to the right of the line formed by the
    notch and v1, false otherwise.
*/
bool on_right(Vertex *v1, Vertex *notch, Vertex *curr)
{
  double a = notch->x - v1->x;
  double b = notch->y - v1->y;
  double c = curr->x - v1->x;
  double d = curr->y - v1->y;

  if (a * d - b * c <= 0)
  {
    return true;
  }
  else
  {
    return false;
  }
}

/**
 * @brief This function checks whether the given vertex lies inside the
 * rectangle
 *
 * @param v given vertex
 * @param rec vector containing the coordinates of the top right and bottom left
 * @return true
 * @return false
 */
// This function checks whether the given vertex lies inside the rectangle
//  formed by the given vertices.
bool rectangle(Vertex *v, vector<double> rec)
{
  if (v->x <= rec[0] && v->x >= rec[1] && v->y <= rec[2] && v->y >= rec[3])
  {
    return true;
  }
  else
  {
    return false;
  }
}

/**
 * @brief This function implements the given alogrithm to convert non-convex
 * polygon to convex polygons.
 *
 * @param polygon as a DCEL
 */
// This function implements the given alogrithm to convert non-convex polygon to
// convex polygons.
void solve(DCEL *&polygon)
{
  vector<Vertex *> L;                    ///> It contains the vertices of the polygon in clockwise
                                         ///> order in current iteration.
  vector<Vertex *> list = polygon->Vert; ///> It contains the vertices left in
                                         /// the > polygon to be processed.
  Vertex *init_vertex = list[0];         ///> This is the first vertex of the polygon.

  /// 2
  Vertex *v1 = init_vertex; ///> This is the first vertex of the current edge.
  Vertex *v2 = polygon->mp[v1]
                   .second; ///> This is the second vertex of the current edge.

  /// 3
  while (list.size() > 3)
  {
    // 3.1 and 3.2
    L.push_back(init_vertex); ///> The first vertex is added to the list.
    v1 = init_vertex;
    v2 = polygon->mp[v1].second;
    L.push_back(v2);

    // 3.3
    ///  We check whether the next vertex creates a notch or not. If it does
    ///  not, we add it to the list.
    while (
        !is_notch(L[L.size() - 2], L[L.size() - 1],
                  polygon->mp[L[L.size() - 1]].second) &&
        !is_notch(L[L.size() - 1], polygon->mp[L[L.size() - 1]].second, v1) &&
        !is_notch(polygon->mp[L[L.size() - 1]].second, v1, v2) &&
        L.size() < list.size())
    {
      L.push_back(polygon->mp[L[L.size() - 1]].second);
    }

    // 3.4
    ///  We check whether L and list are of same size. If yes, then the polygon
    /// remaining polygon is convex and we can stop the algorithm.
    if (list.size() != L.size())
    {
      // 3.4.1
      //  We find the notches which are part of list and not part of L.
      vector<Vertex *> notches;

      for (int i = 0; i < total_notches.size(); i++)
      {
        int flag = 0;
        for (int j = 0; j < L.size(); j++)
        {
          if (total_notches[i]->x == L[j]->x &&
              total_notches[i]->y == L[j]->y)
          {
            flag--;
          }
        }
        for (int j = 0; j < list.size(); j++)
        {
          if (total_notches[i]->x == list[j]->x &&
              total_notches[i]->y == list[j]->y)
          {
            flag++;
            break;
          }
        }

        if (flag == 1)
        {
          notches.push_back(total_notches[i]);
        }
      }

      // 3.4.2
      ///  We remove the notches which lie outside the rectangle.
      vector<double> rec = calc_rectangle(L);
      for (int i = 0; i < notches.size(); i++)
      {
        if (!rectangle(notches[i], rec))
        {
          notches.erase(notches.begin() + i, notches.begin() + i + 1);
          i--;
        }
      }

      /// We remove the vertices which lie to the right of the line formed by
      /// the notches and v1.
      for (int i = 0; i < notches.size(); i++)
      {
        // Calling the point-in-polygon test
        if (!pnpoly(L, notches[i]->x, notches[i]->y))
        {
          continue;
        }

        for (int j = 1; j < L.size() && L.size() > 2; j++)
        {
          bool right = on_right(v1, notches[i], L[j]);
          if (right == true)
          {
            L.erase(L.begin() + j, L.begin() + j + 1);
            j--;
          }
        }
      }
    }
    else
    {
      /// Stop the algorithm after its completion.
      break;
    }

    /// Decompose the polygon by creating a new diagonal.
    // 3.5
    if (L.size() >= 3)
    {
      /// Changing the edges of the original polygon after decomposing it.
      /// 3.5.1
      Edge *edge = new Edge();
      edge->origin = L[0];
      edge->dest = L[L.size() - 1];
      Edge *twin1 = new Edge();
      twin1->origin = L[L.size() - 1];
      twin1->dest = L[0];
      twin1->twin = edge;
      edge->twin = twin1;

      edge->next = polygon->mp[L[L.size() - 1]].first;
      edge->prev = polygon->mp[L[0]].first->prev;
      polygon->mp[L[0]].first->prev->next = edge;
      polygon->mp[L[L.size() - 1]].first->prev = edge;

      edge->twin->next = edge->prev->twin;
      edge->twin->prev = edge->next->twin;
      edge->prev->twin->prev = edge->twin;
      edge->next->twin->next = edge->twin;

      polygon->mp[L[0]] = {edge, L[L.size() - 1]};

      // Adding the new diagonal to the list of diagonals.
      pair<Vertex *, Vertex *> diagonal = {L[0], L[L.size() - 1]};
      diagonals.push_back(diagonal);

      // 3.5.2
      ///  Removing the vertices which are part of L from the list and adding
      ///  the first and last edge of L. list <- (list-L) U {L[0],
      ///  L[L.size()-1]}
      for (int i = 1; i + 1 < L.size(); i++)
      {
        for (int j = 0; j < list.size(); j++)
        {
          if (L[i]->x == list[j]->x && L[i]->y == list[j]->y)
          {
            list.erase(list.begin() + j, list.begin() + j + 1);

            j--;
          }
        }
      }

      /// Creating new face and adding the new face to the list of faces.
      Face *new_face = new Face();
      for (int i = 0; i < L.size(); i++)
      {
        new_face->vertices.push_back(L[i]);
      }
      polygon->faces.push_back(new_face);

      /// Initializing init_vertex for the next iteration.
      init_vertex = L[L.size() - 1];

      /// Clearing L for the next iteration.
      L.clear();
    }
    else
    {
      /// Initializing init_vertex for the next iteration.
      init_vertex = L[L.size() - 1];

      /// Clearing L for the next iteration.
      L.clear();
    }
  }

  /// Adding the last face to the list of faces.
  if (list.size() >= 3)
  {
    Face *new_face = new Face();
    for (int i = 0; i < list.size(); i++)
    {
      new_face->vertices.push_back(list[i]);
    }
    polygon->faces.push_back(new_face);
  }
}
/**
 * @brief This is the merging algorithm when a lies above b and f1 lies to the
 *
 * @param f1 is the face to the left of diagonal ab.
 * @param f2 is the face to the right of diagonal ab.
 * @param a the first vertex of the diagonal ab.
 * @param b the second vertex of the diagonal ab.
 * @return Face* the new face formed after merging f1 and f2.
 */
/// This is the merging algorithm when a lies above b and f1 lies to the left of
/// diagonal ab and f2 lies to the right of diagonal ab.
Face *do_merge_ab(vector<Vertex *> f1, vector<Vertex *> f2, Vertex *a,
                  Vertex *b)
{
  Face *f = new Face();
  vector<Vertex *> ans;
  int flag = 0;
  /// First we go clockwise in f2 and add all the vertices in clockwise order
  /// from a to b in ans.
  for (int i = 0; i < f2.size(); i++)
  {
    if (f2[i]->x == a->x && f2[i]->y == a->y)
    {
      ans.push_back(a);
      int j = (i + 1) % f2.size();
      while (!(f2[j]->x == b->x && f2[j]->y == b->y))
      {
        ans.push_back(f2[j]);
        j = (j + 1) % f2.size();
      }
      ans.push_back(b);
      flag = 1;
    }
    if (flag == 1)
    {
      break;
    }
  }
  flag = 0;
  /// Now we go clockwise in f1 and add all the vertices in clockwise order from
  /// b to a in ans.
  for (int i = 0; i < f1.size(); i++)
  {
    if (f1[i]->x == b->x && f1[i]->y == b->y)
    {

      int j = (i + 1) % f1.size();
      while (!(f1[j]->x == a->x && f1[j]->y == a->y))
      {
        ans.push_back(f1[j]);
        j = (j + 1) % f1.size();
      }
      flag = 1;
    }
    if (flag == 1)
    {
      break;
    }
  }
  f->vertices = ans;
  return f;
}
/**
 * @brief This is the merging algorithm when b lies above a and f1 lies to the
 * left of diagonal ba and f2 lies to the right of diagonal ba.
 *
 * @param f1 the face to the left of diagonal ba.
 * @param f2 the face to the right of diagonal ba.
 * @param a the second vertex of the diagonal ba.
 * @param b the first vertex of the diagonal ba.
 * @return Face* face formed after merging f1 and f2.
 */
/// This is the merging algorithm when b lies above a and f1 lies to the left of
/// diagonal ba and f2 lies to the right of diagonal ba.
Face *do_merge_ba(vector<Vertex *> f1, vector<Vertex *> f2, Vertex *a,
                  Vertex *b)
{
  Face *f = new Face();
  vector<Vertex *> ans;
  int flag = 0;
  // First we go clockwise in f2 and put the vertices in clockwise order from b
  // to a in ans.
  for (int i = 0; i < f2.size(); i++)
  {
    if (f2[i]->x == b->x && f2[i]->y == b->y)
    {
      ans.push_back(b);
      int j = (i + 1) % f2.size();
      while (!(f2[j]->x == a->x && f2[j]->y == a->y))
      {
        ans.push_back(f2[j]);
        j = (j + 1) % f2.size();
      }
      ans.push_back(a);
      flag = 1;
    }
    if (flag == 1)
    {
      break;
    }
  }
  flag = 0;
  /// Then we go clockwise in f1 and put the vertices in clockwise order from a
  /// to b in ans.
  for (int i = 0; i < f1.size(); i++)
  {
    if (f1[i]->x == a->x && f1[i]->y == a->y)
    {

      int j = (i + 1) % f1.size();
      while (!(f1[j]->x == b->x && f1[j]->y == b->y))
      {
        ans.push_back(f1[j]);
        j = (j + 1) % f1.size();
      }
      flag = 1;
    }
    if (flag)
    {
      break;
    }
  }
  f->vertices = ans;
  return f;
}
/**
 * @brief This function merges any two polygons having a common diagonal having
 * no notches at both vertices connected by the diagonal.
 *
 * @param polygon the polygon in which we have to merge the diagonals.
 */
// This function merges any two polygons having a common diagonal having no
// notches at both vertices connected by the diagonal.
void merge(DCEL *polygon)
{
  int m = diagonals.size(); // m is the number of diagonals in the polygon.

  for (int i = 0; i < m; i++)
  {
    Vertex *a = diagonals[i].first; ///> a is the first vertex of the diagonal.
    Vertex *b =
        diagonals[i].second; ///> b is the second vertex of the diagonal.
    Vertex *c = NULL;
    Vertex *d = NULL;
    Vertex *e = NULL;
    Vertex *f = NULL;
    Face *f1 = NULL;
    Face *f2 = NULL;

    /// We find the two faces f1 and f2 which have the diagonal ab/ba as their
    /// common edge.
    for (int j = 0; j < polygon->faces.size(); j++)
    {
      int count = 0;

      for (int k = 0; k < polygon->faces[j]->vertices.size(); k++)
      {
        if (polygon->faces[j]->vertices[k]->x == a->x &&
            polygon->faces[j]->vertices[k]->y == a->y)
        {
          count++;
        }
        else if (polygon->faces[j]->vertices[k]->x == b->x &&
                 polygon->faces[j]->vertices[k]->y == b->y)
        {
          count++;
        }
      }

      if (count == 2 && !f1)
      {
        f1 = polygon->faces[j];
      }
      else if (count == 2)
      {
        f2 = polygon->faces[j];
      }
    }

    int k = f1->vertices.size();

    bool ab = false;

    /// We find the vertices c and d which are adjacent to a and b in f1.
    for (int j = 0; j < f1->vertices.size(); j++)
    {
      if ((f1->vertices[j]->x == a->x && f1->vertices[j]->y == a->y) &&
          (f1->vertices[(j + 1) % k]->x == b->x &&
           f1->vertices[(j + 1) % k]->y == b->y))
      {
        c = f1->vertices[(j - 1 + k) % k];
        d = f1->vertices[(j + 2) % k];
        ab = true;
      }
      else if ((f1->vertices[(j + 1) % k]->x == a->x &&
                f1->vertices[(j + 1) % k]->y == a->y) &&
               (f1->vertices[j]->x == b->x && f1->vertices[j]->y == b->y))
      {
        e = f1->vertices[(j + 2) % k];
        f = f1->vertices[(j - 1 + k) % k];
      }
    }
    k = f2->vertices.size();

    /// We find the vertices e and f which are adjacent to a and b in f2.
    for (int j = 0; j < f2->vertices.size(); j++)
    {
      if ((f2->vertices[j]->x == a->x && f2->vertices[j]->y == a->y) &&
          (f2->vertices[(j + 1) % k]->x == b->x &&
           f2->vertices[(j + 1) % k]->y == b->y))
      {
        c = f2->vertices[(j - 1 + k) % k];
        d = f2->vertices[(j + 2) % k];
      }
      else if ((f2->vertices[(j + 1) % k]->x == a->x &&
                f2->vertices[(j + 1) % k]->y == a->y) &&
               (f2->vertices[j]->x == b->x && f2->vertices[j]->y == b->y))
      {
        e = f2->vertices[(j + 2) % k];
        f = f2->vertices[(j - 1 + k) % k];
      }
    }

    /// If there are no notches at both the vertices connected by the diagonal,
    /// we merge the two faces.
    if (!(is_notch(c, a, e) || is_notch(f, b, d)))
    {
      Face *f;
      // We merge the two faces according to the spatial arrangement of
      // a,b,f1,f2.
      if (ab == true)
      {
        f = do_merge_ab(f1->vertices, f2->vertices, a, b);
      }
      else
      {
        f = do_merge_ba(f1->vertices, f2->vertices, a, b);
      }

      // We remove the two faces f1 and f2 from the list of faces of the
      // polygon.
      for (int j = 0; j < polygon->faces.size(); j++)
      {
        if (polygon->faces[j] == f1)
        {
          polygon->faces.erase(polygon->faces.begin() + j);
          j--;
        }
      }
      for (int j = 0; j < polygon->faces.size(); j++)
      {
        if (polygon->faces[j] == f2)
        {
          polygon->faces.erase(polygon->faces.begin() + j);
          j--;
        }
      }

      // We add the new face f to the list of faces of the polygon.
      polygon->faces.push_back(f);
    }
  }
}

int main()
{
  freopen("input.txt", "r", stdin); // Reading from the file.
  int n;
  cin >> n;

  // Creating Vertices for the given input polygon from the input taken from the
  // file.
  vector<Vertex *> vertices;
  for (int i = 0; i < n; i++)
  {
    Vertex *v = new Vertex();
    double w, u;
    cin >> w >> u;
    v->x = w;
    v->y = u;
    vertices.push_back(v);
  }

  reverse(vertices.begin(), vertices.end());   // Comment if vertices are to be taken in clockwise order.

  // Creating a DCEL pointer named polygon
  DCEL *poly = new DCEL(vertices);
  DCEL *polygon = poly->build_polygon(vertices);

  // Writing the coordinates of the initial polygon face to a file.
  ofstream file1("Polygon.txt");
  for (int i = 0; i < vertices.size(); i++)
  {
    file1 << vertices[i]->x << " " << vertices[i]->y << " ";
  }
  file1 << "\n";
  file1.close();

  // Calculating the total notches of the polygon.
  vector<Vertex *> t = calc_notches(polygon);
  for (int i = 0; i < t.size(); i++)
  {
    total_notches.push_back(t[i]);
  }

  // Running the Decomposition algorithm to create convex polygons.
  solve(polygon);

  // Writing the coordinates of the decomposed faces to a file.
  ofstream file2("Decomposed.txt");
  for (int i = 0; i < polygon->faces.size(); i++)
  {
    for (int j = 0; j < polygon->faces[i]->vertices.size(); j++)
    {
      file2 << polygon->faces[i]->vertices[j]->x << " "
            << polygon->faces[i]->vertices[j]->y << " ";
    }
    file2 << "\n";
  }
  file2.close();

  // Running the Merging algorithm to merge the convex polygons.
  merge(polygon);

  // Writing the coordinates of the merged faces to a file.
  ofstream file3("Merged.txt");

  for (int i = 0; i < polygon->faces.size(); i++)
  {
    for (int j = 0; j < polygon->faces[i]->vertices.size(); j++)
    {
      file3 << polygon->faces[i]->vertices[j]->x << " "
            << polygon->faces[i]->vertices[j]->y << " ";
    }
    file3 << "\n";
  }
  file3.close();

  return 0;
}
