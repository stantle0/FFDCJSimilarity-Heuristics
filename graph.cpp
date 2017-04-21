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

#include "graph.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iterator>


/**********************
 ** ETREMITY METHODS **
 **********************/
void Extremity::print(void)
{
  if (t == UNDEF) // usually a telomere
    printf("T_");
  else
    printf("%d%c", id, t);
}


/******************
 ** EDGE METHODS **
 ******************/
Edge::Edge(Vertex *adj, const char *label) :
  next(NULL),
  prev(NULL),
  adj(adj),
  adjRef(NULL),
  label(NULL),
  ex1(0, Extremity::Type::UNDEF),
  ex2(0, Extremity::Type::UNDEF),
  sibling(NULL)
{
  setLabel(label);
}

Edge::~Edge()
{
  delete[] label; // No problem if NULL
}

void Edge::print(bool printAdj)
{
  if (label)
    printf("%s", label);
  else if (adj->label)
    printf("%s", adj->label);
  else
    printf("%d", adj->id);

  if (label && adj->label && printAdj)
    printf("(%s)", adj->label);
}

Vertex *Edge::getAdj(void) const
{
  return adj;
}

Edge *Edge::getAdjRef(void) const
{
  return adjRef;
}

void Edge::setExtremities(int id1, Extremity::Type t1, int id2, Extremity::Type t2)
{
  ex1.id = adjRef->ex2.id = id1;
  ex1.t = adjRef->ex2.t = t1;
  ex2.id = adjRef->ex1.id = id2;
  ex2.t = adjRef->ex1.t = t2;
}

Extremity Edge::getExtremityFrom(void)
{
  return ex1;
}

Extremity Edge::getExtremityTo(void)
{
  return ex2;
}

void Edge::setLabel(const char *label)
{
  if (this->label)
    delete[] this->label;

  this->label = NULL;
  if (label) {
    int len = strnlen(label, _GRAPH_MAX_LABEL) + 1;
    this->label = new char[len];
    memcpy(this->label, label, len * sizeof(char));
  }
}

bool Edge::incompatible(Edge *e)
{
  if ((ex1.id == e->ex1.id) ^ (ex2.id == e->ex2.id))
    return true;

  if ((ex1.id == e->ex2.id) ^ (ex2.id == e->ex1.id))
    return true;
  
  return false;
}

void Edge::setSibling(Edge *e)
{
  sibling = adjRef->sibling = e; // here we set the sibling of this edge for the 2 objects representing this edge (at it's two endpoints)
}

Edge *Edge::getSibling(void)
{
  return sibling;
}

bool Edge::incident(Vertex *v)
{
  if (adj == v)
    return true;
  if (adjRef->adj == v)
    return true;
  else
    return false;
}


/********************
 ** VERTEX METHODS **
 ********************/
Vertex::Vertex(int id, char direction, int family) :
  id(id),
  degree(0),
  direction(direction),
  part(0),
  family(family),
  label(NULL),
  data(NULL),
  ex1(0, Extremity::Type::UNDEF),
  ex2(0, Extremity::Type::UNDEF)
{
  edges = new Edge(0); // head
}

void Vertex::print(bool printEdges, const char *fname)
{
  Edge *e;

  if (direction)
    putchar(direction > 0 ? '+' : '-');
    
  if (label)
    printf("%s", label);
    //printf("%s[%d]", label, id);
  else
    printf("%d", id);

  if (fname)
    printf("[%s]", fname);
  else if (family)
    printf("[%d]", family);

  if (part)
    printf("(%c)", part);

  if (!printEdges)
    return;
  
  printf(": ");
  for (e = edges->next; e != NULL; e = e->next) {
    e->print();
    if (e->next != NULL) putchar(',');
  }
  putchar('\n');
}

Vertex::~Vertex()
{
  /* this destructor won't remove edges from neighbors lists BUT, when
     this function is called from graph (removeVertex or destructor),
     that is handled properly */
  while(edges->next != NULL)
    removeEdge(edges->next);
  
  delete[] label; // No problem if NULL
  delete edges;   // head

  // IMPORTANT: user must delete void *data contents, if used
}

Edge *Vertex::addEdge(Vertex *adj, const char *label)
{
  Edge *e = new Edge(adj, label);
  e->next = edges->next;
  e->prev = edges;
  edges->next = e;
  if (e->next != NULL)
    e->next->prev = e;
  degree++;
  return e;
}

void Vertex::removeEdge(Edge *e)
{
  e->prev->next = e->next;
  if (e->next != NULL) e->next->prev = e->prev;
  delete e;
  degree--;
}

void Vertex::setLabel(const char *label)
{
  if (this->label)
    delete[] this->label;

  this->label = NULL;
  if (label) {
    int len = strnlen(label, _GRAPH_MAX_LABEL) + 1;
    this->label = new char[len];
    memcpy(this->label, label, len * sizeof(char));
  }
}

Vertex *Vertex::setExtremities(int id1, Extremity::Type t1, int id2, Extremity::Type t2)
{
  ex1.id = id1;
  ex1.t =  t1;
  ex2.id = id2;
  ex2.t = t2;
  return this;
}

bool Vertex::hasExtremity(Extremity ex)
{
  return (ex == ex1 || ex == ex2);
}


/*******************
 ** GRAPH METHODS **
 *******************/
Graph::Graph(const char *label, int maxvertices) :
  n(0),
  maxn(maxvertices),
  m(0),
  lastVid(-1),
  vertices(maxvertices, NULL),
  label(NULL),
  npart{0}
{
  setLabel(label);

  if (maxn < 1) {
    maxn = 128;
    vertices.resize(maxn, NULL);
  }
  fsize.resize(128);
  fname.resize(128);
}

Graph::Graph(Graph &g) :
  Graph(g.label, g.maxn)
{
  for (const auto v : g) {
    Vertex *newv = addVertex(v->id, v->label, v->part, v->family); // there is no way we could copy void *data content
    newv->setDirection(v->getDirection());
    
    Extremity e1 = v->getExtremityLeft();
    Extremity e2 = v->getExtremityRight();
    
    newv->setExtremities(e1.getId(), e1.getType(), e2.getId(), e2.getType());
  }

  for (const auto v : g)
    for (const auto e : *v)
      if (e->adj->id > v->id) { // so we don't add an edge two times
        Edge *sibling = e->getSibling();

        if (sibling && e->getExtremityFrom().getType() == Extremity::HEAD) // so we don't add and edge and its sibling two times
          continue;
        
        Edge *newe = addEdge(v->id, e->adj->id, e->label);
        
        Extremity e1 = e->getExtremityFrom();
        Extremity e2 = e->getExtremityTo();
        newe->setExtremities(e1.getId(), e1.getType(), e2.getId(), e2.getType()); // must use this function to add extremities to crossref

        if (sibling) {
          Edge *newe_sibling = addEdge(sibling->adjRef->adj->id, sibling->adj->id, sibling->label);

          e1 = e->getExtremityFrom();
          e2 = e->getExtremityTo();
          newe_sibling->setExtremities(e1.getId(), e1.getType(), e2.getId(), e2.getType());

          newe->setSibling(newe_sibling);
          newe_sibling->setSibling(newe);
        }
      }
}

Graph::~Graph()
{
  for (int id = 0; id < maxn; id++)
    delete vertices[id]; // no problem if null
  delete[] label;

  for (auto it : fname)
    delete[] it;
}

void Graph::print()
{
  int i;

  if(label)
    printf("##%s##\n", label);
  for (i = 0; i <= lastVid; i++)
    if (vertices[i] != NULL)
      vertices[i]->print(true, familyName(vertices[i]->family));
}

int Graph::getN(void)
{
  return n;
}

int Graph::getM(void)
{
  return m;
}

int Graph::getMaxVertexId(void)
{
  int i;
  for (i = lastVid; i >= 0 && vertices[i] == NULL; i--)
    ;
  return i;
}

const char *Graph::getLabel(void)
{
  return label;
}

void Graph::setLabel(const char *label)
{
  if (this->label)
    delete[] this->label;

  this->label = NULL;
  if (label) {
    int len = strnlen(label, _GRAPH_MAX_LABEL) + 1;
    this->label = new char[len];
    memcpy(this->label, label, len * sizeof(char));
  }
}
  
Vertex *Graph::getVertex(int id)
{
  return vertices[id];
}

Vertex *Graph::getVertex(char label[])
{
  for (int i = 0; i <= lastVid; i++)
    if (vertices[i] != NULL && strncmp(vertices[i]->label, label, _GRAPH_MAX_LABEL) == 0)
      return vertices[i];
  return NULL;
}

Edge *Graph::addEdge(int id1, int id2, const char *label)
{
  if (id1 >= maxn || id2 > maxn)
    return NULL;
  return addEdge(getVertex(id1), getVertex(id2), label);
}

Edge *Graph::addEdge(Vertex *v1, Vertex *v2, const char *label)
{
  Edge *e1, *e2;

  if (v1 == v2 || v1 == NULL || v2 == NULL) // self edges are not allowed, but duplicated edges are
    return NULL;
  
  m++;
  e1 = v1->addEdge(v2, label);
  e2 = v2->addEdge(v1, label);
  e1->adjRef = e2;
  e2->adjRef = e1;
  return e1;
}

Vertex *Graph::addVertex(const char *label, char part, unsigned int family)
{
  int id;
  if (lastVid < maxn - 1) /* if there is empty space at end */
    id = lastVid + 1;
  else                    /* else, we try to find am empty space */
    for (id = 0; id < maxn && vertices[id] != NULL; id++)
      ;

  return addVertex(id, label, part, family);
}

Vertex *Graph::addVertex(int id, const char *label, char part, unsigned int family)
{
  Vertex *v;
  int len;

  if (part < 0) part = 0;
  
  if (id >= maxn) { // Need to resize vertices vector
    while (id >= maxn)
      maxn *= 2;
    vertices.resize(maxn, NULL);
  }
  
  if (vertices[id] != NULL) // vertex with this id already exists
    return NULL;

  while (family >= fsize.size()) // Need to resize family sizes vector
    fsize.resize(fsize.size() * 2);
  fsize[family]++;
  
  if (id > lastVid)
    lastVid = id;

  n++;
  npart[(unsigned int)part]++;
  v = vertices[id] = new Vertex(id);
  v->part = part;
  v->family = family;
  //v->edges and v->degree should be ok

  if (label) {
    len = strnlen(label, _GRAPH_MAX_LABEL) + 1;
    v->label = new char[len];
    memcpy(v->label, label, len * sizeof(char));
  }

  return v;
}

void Graph::removeVertex(int id)
{
  removeVertex(getVertex(id));
}

void Graph::removeVertex(Vertex *v)
{
  if (v == NULL)
    return;

  while (v->edges->next != NULL) // this remove edges from BOTH endpoints
    removeEdge(v->edges->next);
  
  vertices[v->id] = NULL;
  npart[(unsigned int)v->part]--;
  fsize[v->family]--;
  delete v;
  // user must be sure that v belongs to this graph
  n--;
}

void Graph::removeEdge(Edge *e)
{
  Vertex *v1, *v2;
  Edge *e1, *e2;

  if (!e)
    return;
  
  e1 = e;
  e2 = e->adjRef;

  v1 = e2->adj;
  v2 = e1->adj;

  if (e1->getSibling())
    e1->getSibling()->setSibling(NULL); // this sets the sibling for the two endpoints of the edge
  
  v1->removeEdge(e1);
  v2->removeEdge(e2);
  
  m--;
}
 
void Graph::removeEdge(Extremity ex1, Extremity ex2)
{
  for (auto v : *this)
    if (v->hasExtremity(ex1) || v->hasExtremity(ex2)) {
      for (auto e = v->begin(); e != v->end(); )
        if ((e->getExtremityFrom() == ex1 && e->getExtremityTo() == ex2)
            || (e->getExtremityFrom() == ex2 && e->getExtremityTo() == ex1)){
          Edge *rem = *e;
          e++;
          removeEdge(rem);
        }
        else
          e++;
    }
}

int Graph::partSize(char part)
{
  if (part < 0 || part > 127)
    return 0;
  return npart[(unsigned int)part];
}

int Graph::familySize(unsigned int family, char part)
{
  int size = 0;
  for (auto it = begin(family, part); it != end(); it++)
    size++;
  return size;
}

const char *Graph::familyName(unsigned int family)
{
  if (family >= fname.size())
    return NULL;
  return fname[family];
}

void Graph::setFamilyName(unsigned int family, const char *name)
{
  while (family > fname.size() - 1)
    fname.resize(fname.size()*2);

  if (fname[family] != NULL)
    delete[] fname[family];

  fname[family] = new char[strlen(name) + 1];
  strcpy(fname[family], name);
}
