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

/**
 * Abstract class that dynamically maintains H in parallel
 *
 * H stores h vertices, where h is the largest number such
 * that there are at least h vertices with degree greater than or equal to h
 */
class HSet {

  public:
    dynamic_symmetric_graph<dynamic_symmetric_vertex, pbbs::empty>* G;
    size_t hindex;

    /**
     * Constructs HSet given a pointer to a dynamic graph
     *
     * @param _G, an unweighted dynamic_symmetric_graph
     *     graph does not have to be empty
     */
    HSet(dynamic_symmetric_graph<dynamic_symmetric_vertex, pbbs::empty>* _G) {
      G = _G;
      hindex = 0;
    }

    /**
     * Returns a sequence with all of the vertices in H
     *
     * @return sequence of all the vertices in H
     */
    virtual pbbs::sequence<uintE> getH() = 0;

    /**
     * Returns a boolean that determines if a vertex belongs to H or not
     *
     * @param target, the vertex in question
     * @return whether or not target is in H
     */
    virtual bool contains(uintE target) = 0;

    /**
     * Adds all of the new vertices in a given batch to HSet in parallel
     *
     * @param vertices, sequence of vertices to be added
     *     can contain existing vertices (will be ignored)
     *     CANNOT contain duplicate vertices in batch
     * @return the h-index after all the vertex insertions
     */
    virtual uintE insertVertices(pbbs::sequence<uintE> vertices) = 0;

    /**
     * Deletes all of the existing vertices in batch from HSet in parallel
     *
     * @param vertices, sequence of vertices to be deleted
     *     can contain vertices that don't exist yet (will be ignored)
     *     CANNOT contain duplicate vertices in batch
     * @return the h-index after all the vertex deletions
     */
    virtual uintE eraseVertices(pbbs::sequence<uintE> vertices) = 0;

    /**
     * Given a batch edges, inserts all of the new edges in parallel
     * Automatically adds any new vertices in the edge list
     *
     * @param edges, sequence of edges to be added
     *      Adding edge u, v also adds edge v, u since the graph is symmetric
     *      CANNOT contain duplicate  edges
     *          Filter with the getEdges() function from SubgraphCounting.cc
     *      Can contain edges that already exist (will just be ignored)
     *      Edges can contain new vertices (will be added automatically)
     * @return the h-index after adding all the edges
     */
    virtual uintE insertEdges(pbbs::sequence<std::pair<uintE, uintE>> edges) = 0;

    /**
     * Given a batch edges, deletes all of the new edges in parallel
     * Automatically removes any zero degree vertices after the deletion
     *
     * @param edges, sequence of edges to be erased
     *      Erasing edge u, v also erases edge v, u since the graph is symmetric
     *      CANNOT contain duplicate edges
     *          Filter with the getEdges() function from SubgraphCounting.cc
     *      Can contain edges that don't exist (will be ignored)
     * @return the h-index after deleting all the edges
     */
    virtual uintE eraseEdges(pbbs::sequence<std::pair<uintE, uintE>> edges) = 0;

    /**
     * Frees all data structures used in HSet
     */
    virtual void del() = 0;
};








