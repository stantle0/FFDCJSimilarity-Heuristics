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
*/

#include <vector>
#include <utility>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <forward_list>

#include "graph.hpp"
#include "paths-cycles.hpp"


using std::vector;
using std::unordered_set;
using std::unordered_map;
using std::set;
using std::map;
using std::forward_list;
using std::string;



/******************
 ** PATH METHODS **
 ******************/
Path::Path() :
  l(0),
  le(0)
{
  v.resize(INITIAL_CAPACITY, 0);
  e.resize(INITIAL_CAPACITY, 0);
}

// On copy constructor, we don't copy label, since it is temporary
Path::Path(const Path &obj) :
  l(obj.l),
  le(obj.le),
  v(obj.v),
  e(obj.e)
{ }
  
Path::Path(Vertex *v) :
  Path()
{
  addVertex(v);
}

Path::~Path()
{
}

vector<Vertex *> Path::getVertices(void)
{
  vector<Vertex *> copy = v;
  copy.resize(l);
  return copy;
}

vector<Edge *> Path::getEdges(void)
{
  vector<Edge *> copy = e;
  copy.resize(le);
  return copy;
}

int Path::countNullExtremities(void)
{
  int nulls = 0;
  for(int i = 0; i < le; i++) {
    if (e[i]->getExtremityFrom().getType() == Extremity::UNDEF) // If edge links vertices by null a extremity
      nulls++;
    if (e[i]->getExtremityTo().getType() == Extremity::UNDEF) // If edge links vertices by null a extremity
      nulls++;
  }
  return nulls;
}

int Path::countNullAdjacencies(void)
{
  int nulls = 0;

  for(int i = 0; i < l; i++)
    if (v[i]->getExtremityLeft().getType() == Extremity::UNDEF
        && v[i]->getExtremityRight().getType() == Extremity::UNDEF)
      nulls++;
    
  return nulls;
} 
  
void Path::print(void)
{
  for (int i = 0; i < l; i++) {
    v[i]->print(false);
    if (i < le) { printf("--<"); e[i]->print(false); putchar('>'); }
    if (i < l-1) printf("--");
  }
  printf(",(l:%d,%s)\n", l, isCycle() ? "cycle" : "path");
}

void Path::printEdges(void)
{
  for (int i = 0; i < le; i++) {
    e[i]->print();
    if (i < le-1) printf(",");
  }
}

bool Path::consistent(void)
{
  for (int i = 0; i < le; i++) {  // for each edge i in the path
      
    for (int j = 0; j < le; j++)  // and for each edge j != i in the path
      if (j > i && (e[i]->incompatible(e[j]) || // we check if i and j conflicts (e.g. 1h2h and 1h5h)
                    *e[i] == *e[j]))           // or are the same (appears twice)
        return false;
  }

  return true;
}
  
bool Path::consistent(Edge *edge)
{
  bool result;

  if (inPath(edge)) // small optimization
    return false;
    
  addEdge(edge);             // Add edge
  addVertex(edge->getAdj()); // Add vertex
  result = consistent();     // Test consistency
  removeEdge();              // Remove edge
  removeVertex();            // Remove vertex
    
  return result;
}

bool Path::consistent(Path *p)
{
  for (int i = 0; i < le; i++) // for each edge i in the path
    for (int j = 0; j < p->le; j++)  // and for each edge j != i in the path
      if (e[i]->incompatible(p->e[j])) // we check if i and j conflicts (e.g. 1h2h and 1h5h)
        return false;
    
  return true;
}

string Path::signature(void)
{
  int i, j;
  string s;
  vector<Edge *> edges(le);
  Edge *x;

  for (i = 0; i < le; i++) // new vector with edge pointers
    edges[i] = e[i];

  for (i = 1; i < le; i++) { // sort edge pointers by edge label
    x = edges[i];
    for (j = i - 1; j >= 0 && *x < *edges[j]; j--)
      edges[j+1] = edges[j];
    edges[j+1] = x;
  }

  s.reserve(le * 10);
  for (i = 0; i < le; i++) // print on a string the sorted labels
    s.append(edges[i]->getLabel());
    
  return s;
}





/*************************
 ** CYCLESGRAPH METHODS **
 *************************/
CyclesGraph::CyclesGraph(Graph *ag, const char *label, int len) :
  Graph(label, ag->getN())
{
  buildCyclesGraph(ag, len);
}

CyclesGraph::~CyclesGraph()
{
  for (auto v : *this) {
    Path *d = (Path *) v->getData();
    delete d;
  }
}

void CyclesGraph::buildCyclesGraph(Graph *ag, int len)
//CyclesGraph *cyclesGraph(Graph *ag, const char glabel[], int len)//TODO: remove
{
  Vertex *v;
  char part;

  if (ag->getN() < 1 || len < 2) { // We can't find cycles when there are no vertices or the length of cycles is less than 2 (we have no self-edges)
    forward_list<Path *> empty;
    buildCyclesGraph(ag, &empty);
    return;
  }
  
  forward_list<Path *> *newlist, *list = new forward_list<Path *>;       // temporary lists containing paths
  forward_list<Path *> *cycles = new forward_list<Path *>;               // list containing cycles found
  unordered_set<string> cycle_signatures((ag->getN()/2)*(ag->getN()/2)); // hash table size: (n/2)^2

  // assuming we have at least 1 vertex
  part = ag->begin()->getPart();

  // We try to find cycles starting just in one part
  for (auto it = ag->begin(part); it != ag->end(); ++it) {
    v = *it;
    list->push_front(new Path(v));

    for (int i = 0; i < len; i++) {
      newlist = new forward_list<Path *>;
      for (auto p : *list) {          // for each path in list
        for (auto e : *p->last()) {   // we try to add each edge incident to last vertex in path
          
          if (!p->consistent(e))
            continue;
          
          bool cycle = p->isCycle(e);
          if (i < len-1 && !cycle)                               // if not in desired lenght
            newlist->push_front(new Path(*p + e->getAdj() + e));   // for now, |V| = |E+1| in the cycle
          else if (i == len-1 && cycle && *e > *p->firstE())     // if it may close the cycle of desired lenght (optimization)
            newlist->push_front(new Path(*p + e));                 // at end, |V| = |E| in the cycle
        }
      }
      // clear old list
      for (auto p : *list) // just deleting list results in memory leaks, as list elements are pointers and it's destructor won't destroy objects
        delete p;
      delete list;
      list = newlist; // in newlist, all paths have lenght equal to (list lenght) + 1
    }

    for (auto c : *list) {
      string sign = c->signature();
      
      if (cycle_signatures.find(sign) == cycle_signatures.end()) { // cycle generated for the first time
        cycle_signatures.insert(sign); // add signature to hash
        cycles->push_front(c);         // add cycle to set
      }
      else
        delete c; // cycle previously generated
    }
    list->clear();
  }

  delete list;
  buildCyclesGraph(ag, cycles);
  delete cycles; // we don't delete cycles as they are satellite data of cg vertices, we delete just the container
}

void CyclesGraph::buildCyclesGraph(Graph *ag, std::forward_list<Path *> *cycles)
//CyclesGraph *cyclesGraph(Graph *ag, forward_list<Path *> *cycles, const char glabel[])//TODO: remove
{
  
  // 1st level hash map (unordered, faster to access):
  // * key: gene a
  // * value: 2nd level hash map (ordered, because we must iterate on it):
  //   * key: gene b
  //   * value: a list with cycles containing edges associating a with b
  unordered_map<int, map<int, forward_list<Vertex *>>> associations(ag->getN()/2);
    
  for (auto c : *cycles) {
    Vertex *v = addVertex(c->signature().c_str()); // we add vertex representing cycle to CG
    v->setData(c);                                     // store the cycle it represents

    // we keep track of vertices we already added an edge, so we don't
    // add duplicate edges (unordered_set uses a hash table, and since the
    // number of elements in this set will probably be small, we don't
    // want to waste time creating that table, otherwise would be more
    // efficient to create a vector and lookup all elements every time)
    set<Vertex *> added;
        
    // add edges linking vertices that represents cycles C and C' when C U C' is inconsistent
    // (if two cycles share the same edge, there will be other inconsistent edges)
    for (int i = 0; i < c->lenE(); i++) {

      Edge *e = c->nthE(i);
      Extremity from = e->getExtremityFrom(), to = e->getExtremityTo();

      if (from.getType() == Extremity::UNDEF || to.getType() == Extremity::UNDEF)
        continue;
      
      for (auto &&table = associations[from.getId()].cbegin(); table != associations[from.getId()].cend(); ++table) {
         // avoid adding an edge between a pair of vertices representing cycles with siblings edges
        if (table->first == to.getId())
          continue;

        // using auto&& so that I can change the elements by acessing them by reference
        for (auto &&w : table->second)
          if (added.find(w) == added.end()) {
            addEdge(v, w);
            added.insert(w);
          }
      }
      associations[from.getId()][to.getId()].push_front(v);

      for (auto &&table = associations[to.getId()].cbegin(); table != associations[to.getId()].cend(); ++table) {
         // avoid adding an edge between a pair of vertices representing cycles with siblings edges
        if (table->first == from.getId())
          continue;

        for (auto &&w : table->second)
          if (added.find(w) == added.end()) {
            addEdge(v, w);
            added.insert(w);
          }
      }
      associations[to.getId()][from.getId()].push_front(v);
    }
  }

  // Old simpler but inefficient code
  /*
  for (auto c : *cycles) {
    Vertex *v = addVertex(c->signature().c_str()); // we add vertex representing cycle to CG
    v->setData(c);                                     // store the cycle it represents

    // add edges linking vertices that represents cycles C and C' when C U C' is inconsistent
    for (auto w : *cg) {
      Path *d = (Path *) w->getData();
      if (w->getId() != v->getId() && !c->consistent(d)) // we skip this when v = w
        addEdge(v, w);
    }
  }
  */
}



/*******************
 ** OTHER METHODS **
 *******************/
void walk(Graph *ag, Vertex *v)
{
  int op, i, to;
  Path p(v);
  Edge *e = NULL;
 
  do {

    v = p.last();
    printf("\n");
    p.print();
    printf("\n1: list adjacencies\n2: test consistency\n3: walk\n4: print path\n0: exit\n\n");
    scanf("%d", &op);
    
    switch (op) {

    case 1:
      i = 0;
      for (auto e : *v) {
        printf("\t%d: ", i++);
        e->print(); putchar('\n');
      }
      break;

    case 2:
      i = 0;
      for (auto e : *v) {
        printf("\t%d: ", i++);
        e->print(); putchar('\n');
      }
      printf("\t%d: return\n\nWhich? ", i);
      scanf("%d", &to);
      i = 0;
      for (auto e2 = v->begin(); e2 != v->end() && i <= to; e2++, i++)
        e = *e2;
      if (e != NULL)
        printf(p.consistent(e) ? "CONSISTENT\n" : "INCONSISTENT\n");
      break;
      
    case 3:
      i = 0;
      for (auto e : *v) {
        printf("\t%d: ", i++);
        e->print(); putchar('\n');
      }
      printf("\t%d: return\n\nWhere to? ", i);
      scanf("%d", &to);
      i = 0;
      for (auto e2 = v->begin(); e2 != v->end() && i <= to; e2++, i++)
        e = *e2;
      if (e != NULL)
        p.add(e->getAdj(), e);
      break;

    case 4:
      p.print();
      break;
      
    case 0:
      break;
      
    default:
      break;
    }
  } while(op != 0);
}
