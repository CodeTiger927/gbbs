#pragma once

#include <iostream>
#include <map>
#include <set>
#include <limits>
#include <vector>
#include <math.h>

#include "ligra/ligra.h"
#include "dynamic_symmetric_graph.h"
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "pbbslib/sequence.h"


size_t INIT_DEG_SIZE = 4;
size_t INIT_C_SIZE = 10;

class HSet {

  public:
    dynamic_symmetric_graph<dynamic_symmetric_vertex, pbbs::empty>* G;
    size_t hindex;

    HSet(dynamic_symmetric_graph<dynamic_symmetric_vertex, pbbs::empty>* _G) {
      G = _G;
      hindex = 0;
    }

    virtual pbbs::sequence<uintE> getH() = 0;
    virtual bool contains(uintE target) = 0;

    virtual uintE insertVertices(pbbs::sequence<uintE> vertices) = 0;
    virtual uintE eraseVertices(pbbs::sequence<uintE> vertices) = 0;
    virtual uintE insertEdges(pbbs::sequence<std::pair<uintE, uintE>> edges) = 0;
    virtual uintE eraseEdges(pbbs::sequence<std::pair<uintE, uintE>> edges) = 0;
};








