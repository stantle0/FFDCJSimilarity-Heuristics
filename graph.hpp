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
  Unweighted undirected graph library. Vertices and edges can be
  labeled. For each edge, both endpoints store the edge. This is NOT a
  general purpose library (although it may be customized to do
  so). For that, you can see HyperGraph.h from OGDF (Open Graph
  Drawing Framework). This library aims to be used for graphs of DCJ
  problems (gene graph, adjacency graph), and to be efficient in terms
  of inserting/removing vertices and edges. The edge removal can be
  faster if instead of a linked links we use a tree (or C++ std::set)
*/


#ifndef _GRAPH_HPP

#define _GRAPH_HPP 1

#define _GRAPH_MAX_LABEL 100

#include <iterator>
#include <vector>


/* Some forward-declaration */
class Edge;
class Vertex;
class Graph;


/*********************
 ** EXTREMITY CLASS **
 *********************/
class Extremity { /* To store, in adjacency graph, which extremities a vertex
                     represent or to which extremity some edge is incident */
  friend class Edge;
  friend class Vertex;

public:
  enum Type : char {
    TAIL = 't', /* Represents some gene tail */
    HEAD = 'h', /* Represents some gene head */
    UNDEF = '_' /* Undefined extremity, id is irrelevant */
  };

private:
  int id; /* The unique ID (not the family) of the gene this extremity represents */
  Type t; /* The type: tail, head or undefined (see enum above) */

public:
  /* Constructor that receives (optionally) the ID of the gene this extremity represents */
  Extremity(int id = 0, Type t = UNDEF) : id(id), t(t) {}

  /* Returns the gene which this extremity represents */
  inline int getId(void) const { return id; }

  /* Returns the type of the extremity: tail, head or undefined (used in null extremities) */
  inline Type getType(void) const { return t; }

  /* == operator overload */
  inline bool operator==(const Extremity &other) const {
    return (id == other.id && t == other.t) || (t == UNDEF && other.t == UNDEF);
  }

  /* != operator overload */
  inline bool operator!=(const Extremity &other) const {
    return !(*this == other);
  }

  /* Returns an extremity with inverse type */
  inline Extremity operator!() const {
    return Extremity(id, t == UNDEF ? UNDEF : (t == TAIL ? HEAD : TAIL));
  }

  /* Prints the extremity */
  void print(void);
};


/****************
 ** EDGE CLASS **
 ****************/
class Edge {           /* To be used as a doubly linked list, */
  friend class Vertex; /* an edge in a graph is actually two objects, */
  friend class Graph;  /* each one stored in one of the two vertcies */

private:
  Edge *next;    /* Next on list */
  Edge *prev;    /* Previous on list */
  Vertex *adj;   /* Vertex adjacent to */
  Edge *adjRef;  /* Reference to this edge on adjacent vertex's edge list */
  char *label;   /* Stores label */
  Extremity ex1; /* Extremity of vertex where this edge is stored */
  Extremity ex2; /* Extremity of adjacent vertex to the one this edge is stored */
  Edge *sibling; /* Stores this edge's sibling (used on adjacency graph) */

public:
  /* Default constructor */
  Edge(Vertex *adj = 0x0, const char *label = 0x0);

  /* Default destructor, clean edge data and references */
  ~Edge();

  /* Prints an edge */
  void print(bool printAdj = true);

  /* Returns adjacent vertex referenced by this edge */
  Vertex *getAdj(void) const;

  /* Returns this edge, but the one stored in neighbor vertex  */
  Edge *getAdjRef(void) const;

  /* Returns label */
  inline const char *getLabel(void) { return this->label; }

  /* Sets label */
  void setLabel(const char *label);

  /* Sets its extremities */
  void setExtremities(int id1, Extremity::Type t1, int id2, Extremity::Type t2);

  /* Returns first extremity */
  Extremity getExtremityFrom(void);

  /* Returns second extremity */
  Extremity getExtremityTo(void);

  /* Returns true if extremities of this edge conflicts with the ones of passed edge
     Two extremities i and j conflicts with extremities k and l if (we ignore if tail/head):
      * i = k XOR j = l, or
      * i = l XOR j = k
     For example, 1t2t conflicts with 2h5h but not with 2h1h or 3t4t
  */
  bool incompatible(Edge *e);

  /* Sets this edge's sibling */
  void setSibling(Edge *e);

  /* Gets this edge's sibling */
  Edge *getSibling(void);

  /* Returns true if this edge is incident to v */
  bool incident(Vertex *v);

  /* Overload for comparing edge extremities alphabetically */
  inline bool operator<=(const Edge &other) const;

  /* Overload for comparing edge extremities alphabetically */
  inline bool operator<(const Edge &other) const;

  /* Overload for comparing edge extremities alphabetically */
  inline bool operator>(const Edge &other) const;

  /* Overload for comparing edge extremities alphabetically */
  inline bool operator>=(const Edge &other) const;

  /* Overload for comparing if the edges are the same (NOT EXTREMITIES) */
  inline bool operator==(const Edge &other) const;
};



/******************
 ** VERTEX CLASS **
 ******************/
class Vertex {   /* To be used as an array */
  friend class Edge;
  friend class Graph;

public:

private:
  int id;                     /* Vertex id (should be equal to array index) */
  unsigned short degree;      /* Vertex degree, optional */
  char direction;             /* Gene direction, optional (1: -->, -1: <--, 0: unoriented */
  unsigned char part;         /* Which part of graph this vertex belongs, optional */
  unsigned short family;      /* Family id, 0 = no family */
  char *label;                /* Stores string */
  Edge *edges;                /* Doubly-linked list with a head */
  void *data;                 /* Arbitrary satellite data, user must destroy it
                                 since we can't call delete to a void pointer */
  Extremity ex1;              /* Left extremity */
  Extremity ex2;              /* Right extremity */

public:
  /* Default constructor */
  Vertex(int id = -1, char direction = 0, int family = 0);

  /*
    Default destructor. This shouldn't be called by user, just by
    graph class for proper removal of all edge references (the caller
    must also remove the edges to it from other endpoints)
  */
  ~Vertex();

  /* Prints a vertex */
  void print(bool printEdges = true, const char *fname = 0x0);

  /* Gets vertex id */
  inline int getId(void) const { return id; }

  /* Returns the first extremity of the adjacency (left extremity) */
  inline Extremity getExtremityLeft(void) const { return ex1; }

  /* Returns the second extremity of the adjacency (right extremity) */
  inline Extremity getExtremityRight(void) const { return ex2; }

  /* Sets its extremities */
  Vertex *setExtremities(int id1, Extremity::Type t1, int id2, Extremity::Type t2);

  /* Returns true if the vertex has some extremity equal to ex */
  bool hasExtremity(Extremity ex);

  /* Returns the direction of the gene */
  inline char getDirection(void) const { return direction; }

  /* Sets the direction of the gene */
  inline void setDirection(char direction) { this->direction = direction; }

  /* Returns the vertex degree */
  inline int getDegree(void) const { return degree; }

  /* Returns the part of the bipartite graph this vertex belongs (not used if not bipartite) */
  inline char getPart(void) const { return part; }

  /* Sets the part of the bipartite graph this vertex belongs */
  inline void setPart(char part) { this->part = part; }

  /* Returns the pointer to the arbitrary data stored by the void pointer */
  inline void *getData(void) const { return data; }

  /* Sets the pointer to the arbitrary data */
  inline void setData(void *data) { this->data = data; }

  /* Returns the vertex label */
  inline const char *getLabel(void) const { return this->label; }

  /* Sets the vertex label */
  void setLabel(const char *label);

  /* Add an edge to this vertex (the caller must also add the edge to other endpoint) */
  Edge *addEdge(Vertex *adj, const char *label = 0x0);

  /* Remove an edge from this vertex (the caller must also remove the edge from other endpoint) */
  void removeEdge(Edge *e);

  /* Iterator (over vertices) class and associated methods */
  class iterator : public std::iterator<std::forward_iterator_tag, Edge>
  {
  private:
    Edge *cur;

  public:
    inline iterator(Edge *cur);
    inline iterator(const iterator& i);
    inline iterator& operator=(const iterator& i);
    inline iterator& operator++();
    inline iterator operator++(int);
    inline Edge *operator*() const;
    inline Edge *operator->() const;
    inline bool operator==(const iterator& i) const;
    inline bool operator!=(const iterator& i) const;
  };

  inline iterator begin();
  inline iterator end();
};



/****************************
 ** UNDIRECTED GRAPH CLASS **
 ****************************/
class Graph {

private:
  int n;                          /* Number of vertices */
  int maxn;                       /* Max number of vertices (also represents greater-id-possible + 1)*/
  int m;                          /* Number of edges */
  int lastVid;                    /* Last (greater) id used on adding a vertex (incremental) */
  std::vector<Vertex *> vertices; /* Vertex array */
  char *label;                    /* Optional graph label */
  int npart[128];                 /* Number of vertices on each part */
  std::vector<int> fsize;         /* Size of each family */
  std::vector<char *> fname;      /* Name of each family */

public:
  /*
    Initialize an empty graph. A initial max number of vertices must
    be provided, vertices array will be resized (automatically) beyond
    this limit.
  */
  Graph (const char *label = 0x0, int maxvertices = 0);

  /* Copy constructor */
  Graph (Graph &g);

  /* Destroy a graph, freeing all allocated memory */
  ~Graph ();

  /* Print a graph, use carefully with big graphs */
  void print();

  /* Returns the number of vertices */
  int getN(void);

  /* Returns the number of edges */
  int getM(void);

  /* Returns the greater vertex id */
  int getMaxVertexId(void);

  /* Returns a pointer to vertex with id */
  Vertex *getVertex(int id);

  /* Returns a pointer to vertex with label */
  Vertex *getVertex(char label[]);

  /* Returns label */
  const char *getLabel(void);

  /* Sets label */
  void setLabel(const char *label);

  /* Add an edge, returning it (v1 to v2) if added or NULL */
  Edge *addEdge(int id1, int id2, const char *label = 0x0);

  /* Add an edge, returning it (v1 to v2) if added or NULL */
  Edge *addEdge(Vertex *v1, Vertex *v2, const char *label = 0x0);

  /* Remove edge from graph */
  void removeEdge(Edge *e);

  /* Remove edge with this extremities from graph (not so cheap) */
  void removeEdge(Extremity ex1, Extremity ex2);

  /*
    Add to graph a vertex. The chosen id it the next not used
    (possibly lastVid + 1), part = 0 means NO SPECIFIC PART

    Returns:
      NOT NULL - vertex added (and the address returned is the new vertex)
      NULL - vertex not added (out of memory or max vertices reached)
  */
  Vertex *addVertex(const char *label = 0x0, char part = 0, unsigned int family = 0);

  /*
   Add to graph a vertex with specific id. Argument id must be between
   0 and maxn-1. lastVid is updated to id if it's greater.

   Returns:
     NOT NULL - vertex added (and the address returned is the new vertex)
     NULL - vertex not added (out of memory or max vertices reached)
  */
  Vertex *addVertex(int id, const char *label = 0x0, char part = 0, unsigned int family = 0);

  /*
    Remove from graph a vertex, cleaning all it's data, edge
    endpoints and cross references
  */
  void removeVertex(int id);

  /*
    Remove from graph a vertex, cleaning all it's data, edge
    endpoints and cross references
  */
  void removeVertex(Vertex *v);

  /*
    Returns part size
  */
  int partSize(char part);

  /*
    Returns family size (optionally, just in some part)
  */
  int familySize(unsigned int family, char part = -1);

  /*
    Returns family name or NULL if family name still not defined
  */
  const char *familyName(unsigned int family);

  /*
    (Re)Sets family name
  */
  void setFamilyName(unsigned int family, const char *name);


  /* Iterator (over vertices) class and associated methods */
  class iterator : public std::iterator<std::forward_iterator_tag, Vertex>
  {
  private:
    Graph *g;
    int cur;
    char part;
    int family;

  public:
    inline iterator(Graph *g, int cur, char part = -1, int family = -1);
    inline iterator(const iterator& i);
    inline iterator& operator=(const iterator& i);
    inline iterator& operator++();
    inline iterator operator++(int);
    inline Vertex *operator*() const;
    inline Vertex *operator->() const;
    inline bool operator==(const iterator& i) const;
    inline bool operator!=(const iterator& i) const;
  };

  inline iterator begin();
  inline iterator begin(char part);
  inline iterator begin(char part, unsigned int family);
  inline iterator begin(unsigned int family);
  inline iterator begin(Vertex *);
  inline iterator begin(int id);
  inline iterator end();

private:
  inline iterator _begin(char part, int family, int id = 0); // private, so users won't call with family < 0 (-1 = any)
};


/***************************************************
 ** GRAPH ITERATOR (OVER VERTICES) INLINE METHODS **
 ***************************************************/

inline Graph::iterator Graph::begin()
{
  return _begin(-1, -1);
}

inline Graph::iterator Graph::begin(char part)
{
  return _begin(part, -1);
}

inline Graph::iterator Graph::begin(char part, unsigned int family)
{
  return _begin(part, family);
}

inline Graph::iterator Graph::begin(unsigned int family)
{
  return _begin(-1, family);
}

inline Graph::iterator Graph::begin(Vertex *v)
{
  return begin(v->getId());
}

inline Graph::iterator Graph::begin(int id)
{
  return _begin(-1, -1, id >= 0 ? id : 0);
}

inline Graph::iterator Graph::_begin(char part, int family, int id)
{
  for (; id < maxn; id++)
    if (vertices[id] != NULL &&
        (part == -1 || part == vertices[id]->part) &&
        (family == -1 || family == vertices[id]->family) )
      break;
  return iterator(this, id, part, family);
}

inline Graph::iterator Graph::end()
{
  return iterator(this, maxn, -1, -1);
}

inline Graph::iterator::iterator(Graph *g, int cur, char part, int family) :
  g(g),
  cur(cur),
  part(part),
  family(family)
{}

inline Graph::iterator::iterator(const iterator& i) :
  g(i.g),
  cur(i.cur),
  part(i.part),
  family(i.family)
{}

inline Graph::iterator& Graph::iterator::operator=(const iterator& i)
{
  *this=i; // copy all contents
  return *this;
}

inline Graph::iterator& Graph::iterator::operator++()
{
  for (++cur; cur < g->maxn; ++cur)
    if (g->vertices[cur] != NULL &&
        (part == -1 || part == g->vertices[cur]->part) &&
        (family == -1 || family == g->vertices[cur]->family) )
      break;
  return *this;
}

inline Graph::iterator Graph::iterator::operator++(int)
{
  iterator tmp(*this);
  for (++cur; cur < g->maxn; ++cur)
    if (g->vertices[cur] != NULL &&
        (part == -1 || part == g->vertices[cur]->part) &&
        (family == -1 || family == g->vertices[cur]->family) )
      break;
  return tmp;
}

inline Vertex* Graph::iterator::operator*() const
{
  return g->vertices[cur];
}

inline Vertex* Graph::iterator::operator->() const
{
  return g->vertices[cur];
}

inline bool Graph::iterator::operator==(const iterator& i) const
{
  return g == i.g && cur == i.cur; // part or family doesn't matter
}

inline bool Graph::iterator::operator!=(const iterator& i) const
{
  return g != i.g || cur != i.cur; // part or family doesn't matter
}


/*************************************************
 ** VERTEX ITERATOR (OVER EDGES) INLINE METHODS **
 *************************************************/

inline Vertex::iterator Vertex::begin()
{
  return iterator(edges->next);
}

inline Vertex::iterator Vertex::end()
{
  return iterator(NULL);
}

inline Vertex::iterator::iterator(Edge *cur) :
  cur(cur)
{}

inline Vertex::iterator::iterator(const iterator& i) :
  cur(i.cur)
{}

inline Vertex::iterator& Vertex::iterator::operator=(const iterator& i)
{
  cur=i.cur;
  return *this;
}

inline Vertex::iterator& Vertex::iterator::operator++()
{
  cur = cur->next;
  return *this;
}

inline Vertex::iterator Vertex::iterator::operator++(int)
{
  iterator tmp(*this);
  cur = cur->next;
  return tmp;
}

inline Edge* Vertex::iterator::operator*() const
{
  return cur;
}

inline Edge* Vertex::iterator::operator->() const
{
  return cur;
}

inline bool Vertex::iterator::operator==(const iterator& i) const
{
  return cur == i.cur;
}

inline bool Vertex::iterator::operator!=(const iterator& i) const
{
  return cur != i.cur;
}


/*************************
 ** EDGE INLINE METHODS **
 *************************/

inline bool Edge::operator<=(const Edge &other) const
{
  if (this == &other || this == other.getAdjRef())
    return true;

  return *this < other;
}

inline bool Edge::operator<(const Edge &other) const
{
  if (this == &other || this == other.getAdjRef())
    return true;

  int e1[2] = {ex1.id, ex2.id}, e2[2] = {other.ex1.id, other.ex2.id};

  if (e1[0] > e1[1]) {
    e1[0] = ex2.id;
    e1[1] = ex1.id;
  }
  if (e2[0] > e2[1]) {
    e2[0] = other.ex2.id;
    e2[1] = other.ex1.id;
  }

  // so we won't have problems with null adjacencies
  if (ex1.t ==  Extremity::UNDEF)
    return true;
  if (other.ex1.t ==  Extremity::UNDEF)
    return false;

  if (e1[0] < e2[0])
    return true;
  else if (e1[0] > e2[0])
    return false;
  else if (e1[1] < e2[1]) // here and below, e1[0] == e2[0]
    return true;
  else if (e1[1] > e2[1])
    return false;
  else                    // here and below, e1[1] == e2[1] too
    return ex1.t == Extremity::TAIL; // then, first is less than second if its tail
}

inline bool Edge::operator>(const Edge &other) const
{
  return !(*this <= other);
}

inline bool Edge::operator>=(const Edge &other) const
{
  if (this == &other || this == other.getAdjRef())
    return true;

  return *this > other;
}

inline bool Edge::operator==(const Edge &other) const
{
  return (this == &other || this == other.getAdjRef());
}

#endif /* graph.hpp  */
