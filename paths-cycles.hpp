/**
   Copyright (C) 2016 Diego Rubert
   
   This file is part of the O(k)-approximation algorithm implementation
   for the DCJ distance on linear unichromosomal genomes in:

   Diego P. Rubert, Pedro Feijão, Marília D. V. Braga, Jens Stoye and Fábio V. Martinez
   A Linear Time Approximation Algorithm for the DCJ Distance for Genomes with Bounded Number of Duplicates

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
  This file implements a representation of (possibily closed) paths
  and procedures that collect them. Specific to the Adjacency Graph.

  Here is also defined the CycleGraph, a graph where every vertex
  represents a cycle in an adjcency graph.
*/

#ifndef _PATHS_CYCLES_HPP

#define _PATHS_CYCLES_HPP 1

#include <vector>
#include <utility>
#include <forward_list>

#include "graph.hpp"



/****************
 ** PATH CLASS **
 ****************/
// This class stores a consistent path (initial container size/capacity = 4)
// If we are storing a cycle, we don't store the first vertex again at end of the closed path
class Path {
private:
  int l;              // Length (number of **vertices**, not edges)
  int le;             // Length in edges (number of edges)
  std::vector<Vertex *> v; // Vertices in path  
  std::vector<Edge *> e;   // Edges in path (we need to store edges used too, because we may be working on a multigraph)

  const static int INITIAL_CAPACITY = 4;

public:
  enum Type {
    EVEN = 0,
    ODD,
    ANY
  };

  // Default empty constructor
  Path();

  // Copy constructor
  Path(const Path &obj);

  // Constructor that receives one vertex as an initial path of length 0
  Path(Vertex *v);

  // Destructor
  ~Path();
  
  // Checks if vertex is already in path
  inline bool inPath(Vertex *v);

  // Checks if vertex is already in path
  inline bool inPath(int id);

  // Checks if edge is already in path
  inline bool inPath(Edge *e);

  // Checks if exists some edge with both extremities in path
  inline bool inPath(Extremity ex1, Extremity ex2);

  // Adds vertex to path, returning new length
  inline int addVertex(Vertex *v);

  // Removes vertex from path, returning new length
  inline int removeVertex(void);
  
  // Adds vertex + edge to path, returnin new length (in vertices)
  inline int add(Vertex *v, Edge *e);
    
  // Replaces some vertex (use with caution, pos must be within bounds)
  inline void replace(int pos, Vertex *v);
  
  // Adds edge to path, returning new edges count
  inline int addEdge(Edge *e);

  // Removes last edge from path, returning new edges count
  inline int removeEdge(void);
  
  // Returns n-th (starting in 0) vertex id in path
  inline Vertex *nth(int n);

  // Returns last vertex id in path
  inline Vertex *last(void);

  // Returns first vertex id in path
  inline Vertex *first(void);

  // Returns path length (in **vertices**, not edges)
  inline int len(void);

  // Returns n-th (starting in 0) edge in path
  inline Edge *nthE(int n);

  // Returns last edge in path
  inline Edge *lastE(void);

  // Returns first edge in path
  inline Edge *firstE(void);

  // Returns path length in edges (edges count, may be zero)
  inline int lenE(void);
  
  // Returns vertices vector
  std::vector<Vertex *> getVertices(void);

  // Returns edges vector
  std::vector<Edge *> getEdges(void);

  // Returns the number of null extremities used in this path
  int countNullExtremities(void);

  // Returns the number of null adjacencies in this path
  int countNullAdjacencies(void);

  // Returns true if path is a cycle or false otherwise
  inline bool isCycle(void);

  // Returns true if current path + edge is a cycle or false otherwise
  inline bool isCycle(Edge *e);
  
  // Prints path
  void print(void);

  // Prints path edges
  void printEdges(void);
  
  // Tests whether current path/cycle is consistent or not
  bool consistent(void);
  
  // Tests whether current path + edge is consistent or not
  bool consistent(Edge *edge);

  // Tests whether (this path U path p) is consistent (**assuming** that each one is consistent and they aren't the same)
  bool consistent(Path *p);

  // Generates and returns string with path signature (sorted edge labels)
  std::string signature(void);

  // += edge operator overload
  inline Path &operator+=(const Edge &rhs);

  // += edge* operator overload
  inline Path &operator+=(Edge *rhs);

  // + edge operator overload
  inline const Path operator+(const Edge &other) const;

  // + edge* operator overload
  inline const Path operator+(Edge *other) const;

  // += vertex operator overload
  inline Path &operator+=(const Vertex &rhs);

  // += vertex* operator overload
  inline Path &operator+=(Vertex *rhs);

  // + vertex operator overload
  inline const Path operator+(const Vertex &other) const;

  // + vertex* operator overload
  inline const Path operator+(Vertex *other) const;

  // Iterator on path vertices (if it's a cycle, we usually don't store the last: it is equal to the first)
  typedef std::vector<Vertex *>::iterator iterator;
  iterator begin() { return v.begin(); }
  iterator end() { return v.begin() + l; }
};


/**********************
 ** CYCLEGRAPH CLASS **
 **********************/
// Graph where every vertex represents a cycle in an adjcency graph
class CyclesGraph : public Graph
{
private:
  /* Receives an Adjacency Graph, building a graph whose vertices
   represent cycles of lenght len. In paths used internally, **we don't
   add the first vertex again at end of path**. General steps:

   1. Starting in each vertex, we find all simple paths of lenght len-1
   2. We allow to close cycles just whe their length = len
   2. For each consistent cycle resulting, we generate its signature
   4. We add signatures to a hash table, so we don't add the same
   cycle 2 times to final set

   We generate at most 2 times each cycle for each vertex in the cycle from first part

   Optimization:
   * We just close cycles when the label of closing edge is greater than the first edge,
   since we build every cycle in both directions

   Notes:
   * Its hard to detect duplicated cycles, here we generate a signature by
   ordering the labels of edges in cycle together with a hash table
  */
  void buildCyclesGraph(Graph *ag, int len);

  // Auxiliary function, receives a list of cycles and build a graph
  // representing the packing of cycles
  void buildCyclesGraph(Graph *ag, std::forward_list<Path *> *cycles);
  
public:
  // Default constructor, receives the corresponding adjacency graph,
  // the label and the length of cycles we want to pack
  CyclesGraph(Graph *ag, const char *label = 0x0, int len = 0);

  // Destructor
  ~CyclesGraph();
};


/*******************
 ** OTHER METHODS **
 *******************/
// An interactive walk in the graph, for debugging purposes
void walk(Graph *ag, Vertex *v);



/*************************
 ** PATH INLINE METHODS **
 *************************/
inline bool Path::inPath(Vertex *v)
{
  for (int i = 0; i < l; i++) 
    if (this->v[i] == v)
      return true;
  return false;
}

inline bool Path::inPath(int id)
{
  for (int i = 0; i < l; i++) 
    if (v[i]->getId() == id)
      return true;
  return false;
}

inline bool Path::inPath(Edge *e)
{
  for (int i = 0; i < le; i++) 
    if (this->e[i] == e || this->e[i] == e->getAdjRef())
      return true;
  return false;
}

inline bool Path::inPath(Extremity ex1, Extremity ex2)
{
  for (int i = 0; i < le; i++) 
    if ( (e[i]->getExtremityFrom() == ex1 && e[i]->getExtremityTo() == ex2)
         || (e[i]->getExtremityFrom() == ex2 && e[i]->getExtremityTo() == ex1) )
      return true;
  return false;
}

inline int Path::addVertex(Vertex *v)
{
  if ((unsigned)l == this->v.size())
    this->v.resize(this->v.size()*2, 0);
      
  this->v[l++] = v;
  return l;
}

inline int Path::removeVertex(void)
{
  if ((unsigned)l < this->v.size() / 2)
    this->v.resize(this->v.size() / 2, 0);
      
  return --l;
}
  
inline int Path::add(Vertex *v, Edge *e)
{
  addEdge(e);
  return addVertex(v);
}
    
inline void Path::replace(int pos, Vertex *v)
{
  this->v[pos] = v;
}
  
inline int Path::addEdge(Edge *e)
{
  if ((unsigned)le == this->e.size())
    this->e.resize(this->e.size()*2, 0);
    
  this->e[le++] = e;
  return le;
}

inline int Path::removeEdge(void)
{
  if ((unsigned)le < this->e.size() / 2)
    this->e.resize(this->e.size() / 2, 0);
    
  return --le;
}
  
inline Vertex *Path::nth(int n)
{
  return v[n];
}

inline Vertex *Path::last(void)
{
  return v[l-1];
}

inline Vertex *Path::first(void)
{
  return v[0];
}

inline int Path::len(void)
{
  return l;
}

inline Edge *Path::nthE(int n)
{
  return e[n];
}

inline Edge *Path::lastE(void)
{
  return e[le-1];
}

inline Edge *Path::firstE(void)
{
  return e[0];
}

inline int Path::lenE(void)
{
  return le;
}

inline bool Path::isCycle(void)
{
  // we do not check firstEx() == lastEx() because we may return to first by a different extremity
  if (len() > 1 && len() == lenE()+1 && first() == last()) // when we add first vertex again at end, |V| = |E|-1
    return true;
  else
    if (lenE() > 1 && len() == lenE() && lastE()->getAdj() == first()) // when we don't add first vertex again at end, |V| = |E|
      return true;
    else
      return false;
}

inline bool Path::isCycle(Edge *e)
{
  if (len() == lenE()+1 && e->getAdj() == first())
    return true;
  else
    return false;
}

inline Path &Path::operator+=(const Edge &rhs)
{
  addEdge((Edge *)&rhs);
  return *this;
}

inline Path &Path::operator+=(Edge *rhs)
{
  addEdge(rhs);
  return *this;
}

inline const Path Path::operator+(const Edge &other) const
{
  return Path(*this) += other;
}

inline const Path Path::operator+(Edge *other) const
{
  return Path(*this) += other;
}

inline Path &Path::operator+=(const Vertex &rhs)
{
  addVertex((Vertex *)&rhs);
  return *this;
}

inline Path &Path::operator+=(Vertex *rhs)
{
  addVertex(rhs);
  return *this;
}

inline const Path Path::operator+(const Vertex &other) const
{
  return Path(*this) += other;
}

inline const Path Path::operator+(Vertex *other) const
{
  return Path(*this) += other;
}




#endif /* paths-cycles.hpp  */
