#pragma once

#include "hindex.h"

/**
 * HSet implementation that dynamically maintains H in parallel using a dyn_arr to bucket the vertices
 *
 * H stores h vertices, where h is the largest number such
 * that there are at least h vertices with degree greater than or equal to h
 */
class HSetDynArr : public HSet {

  private:

    uintE bSize; //Number of elements with degree hindex and are inside H
    pbbslib::dyn_arr<pbbslib::dyn_arr<uintE>*> C; //Sorts all vertices by degree

    /**
     * Given a batch of vertices, adds all of them to HSet
     * Automatically resizes everything by itself
     *
     * @param batch, sequence of vertices to be added to HSet
     *     All vertices in batch must be in the graph
     *     Can be empty
     *     CANNOT add vertex if it is already tracked by HSet (e.g. present in C)
     * @param sorted, optional oolean determining if batch is already sorted by degree in nonascending order
     *     Set to false by default
     *     Saves time by not resorting batch if it's already sorted
     */
    uintE insert(sequence<uintE> batch, bool sorted = false) {

      //Exits if there are no vertices in batch
      if (batch.size() == 0) return this->hindex;

      //Sortes vertices by degree in descending order  
      pbbs::sequence<uintE> sortedBatch;
      if (!sorted) {
        sortedBatch = integer_sort(batch, [&] (uintE v) { return this->G->get_vertex(v).degree; }).rslice();
      }

      else {
        sortedBatch = batch;
      }
      batch.clear();

      //Stores degrees in descending order
      sequence<uintE> deg = pbbs::sequence<uintE>(sortedBatch.size());
      par_for(0, sortedBatch.size(), [&] (size_t i) {
        deg[i] = deg[i] = this->G->get_vertex(sortedBatch[i]).degree;
      });

      //Adds additional space to C if it is too small
      size_t oldSize = C.size;
      C.size = std::max((uintE) C.size, deg[0] + 1);
      if (C.size >= C.capacity) {
        C.resize(0);
        par_for(oldSize, C.capacity, [&] (size_t i) {
          C.A[i] = nullptr;
        });
      }

      //--------------------------Updating Variables--------------------------//
      //Current number of vertices with degree >= |H| is |H| - |B| + |C[|H|]|
      uintE aboveH = this->hindex - bSize;
      if (this->hindex < C.size && C.A[this->hindex] != nullptr) aboveH += C.A[this->hindex]->size;

      //Number of vertices being added in batch with degree >= current |H|
      uintE added = pbbs::binary_search(deg, [&] (uintE x) { //Number of vertices that are above new hindex
        return x >= this->hindex;
      });

      aboveH += added;

      addToC(sortedBatch, deg); //Updates C

      //--------------------------Compute H-index--------------------------//
      //Begins constructing sum array where each number indicates the number of vertices lost in H if hindex increases
      pbbs::sequence<uintE> sum = pbbs::sequence<uintE>(sortedBatch.size() + 1);
      par_for(0, sortedBatch.size(), [&] (size_t i) {
        if (this->hindex + i < C.size && C.A[this->hindex + i] != nullptr) {
          sum[i] = C.A[this->hindex + i]->size;

        }
        else sum[i] = 0;
      });

      //Computes Prefix Sum (note scan_add_inplace is exclusive)
      pbbslib::scan_add_inplace(sum);
      par_for(0, sum.size(), [&] (size_t i) {
        if (this->hindex + i < C.size && C.A[this->hindex + i] != nullptr) {
          sum[i] += C.A[this->hindex + i]->size; //Compute inclusive sum
        }
      });

      //Creates the binary array described as M
      //Create M array
      pbbs::sequence<bool> M = pbbs::sequence<bool>(sum.size());
    
      par_for(0, M.size(), [&] (size_t i) {
        M[i] = aboveH - sum[i] >= this->hindex + i + 1;
      });

      //Find index of first 0
      pbbs::sequence<size_t> indexM = pbbs::sequence<size_t>(M.size());
      par_for(0, M.size(), [&] (size_t i) {
        indexM[i] = !M[i] ? i : UINT_E_MAX;
      });
      uintE hindexIncrease = pbbslib::reduce_min(indexM);
      if (hindexIncrease == UINT_E_MAX) hindexIncrease = sortedBatch.size();

      this->hindex += hindexIncrease;

      //--------------------------Updating Variables--------------------------//

      //Calculating B = hindex + C[hindex] - aboveH
      if (hindexIncrease == 0) {
        if (this->hindex >= C.size || C.A[this->hindex] == nullptr) bSize = 0;
        else {
          bSize = this->hindex + C.A[this->hindex]->size - aboveH;
        }
      }
      else {
        aboveH -= sum[hindexIncrease - 1];

        if (this->hindex >= C.size || C.A[this->hindex] == nullptr) bSize = 0;
        else {
          bSize = this->hindex + C.A[this->hindex]->size - aboveH;
        }
      }
    
      //Slight variation to hindex paper - C will store all vertices (instead of the vertices not in B)


      deg.clear();
      sum.clear();

      return this->hindex;
    }

    /**
     * Given a batch of vertices, removes all of them to HSet
     * Automatically resizes everything by itself
     *
     * @param batch, sequence of vertices to be added to HSet
     *     All vertices in batch must be in the graph and in HSet
     *     Can be empty
     *     CANNOT remove a vertex if it is not tracked by HSet (e.g. not present in C)
     * @param sorted, optional boolean determining if batch is already sorted by degree in nonascending order
     *     Set to false by default
     *     Saves time by not resorting batch if it's already sorted
     */
    uintE erase(sequence<uintE> batch, bool sorted = false) {

      pbbs::sequence<uintE> sortedBatch;
      if (!sorted) {
        sortedBatch = integer_sort(batch, [&] (uintE v) { return this->G->get_vertex(v).degree; }).rslice();
      }
      else {
        sortedBatch = batch;
      } 
      batch.clear();
      //Stores degrees and sorts in ascending order
      sequence<uintE> deg = pbbs::sequence<uintE>(sortedBatch.size());
      par_for(0, sortedBatch.size(), [&] (size_t i) {
        deg[i] = this->G->get_vertex(sortedBatch[i]).degree;
      });


      //--------------------------Updating Variables--------------------------//
      uintE aboveH = this->hindex - bSize;

      if (this->hindex < C.size && C.A[this->hindex] != nullptr) aboveH += C.A[this->hindex]->size;

      uintE removed = pbbs::binary_search(deg, [&] (uintE x) { //Number of vertices that are above new hindex
        return x >= this->hindex;
      });

      aboveH -= removed;

      removeFromC(sortedBatch, deg);


      //--------------------------Compute H-index--------------------------//
      //Compute prefix sum
      pbbs::sequence<uintE> sum = pbbs::sequence<uintE>(std::min(sortedBatch.size(), this->hindex));
      par_for(0, sum.size(), [&] (size_t i) {
        if (this->hindex - i - 1 < C.size && C.A[this->hindex - i - 1] != nullptr) { 
          sum[i] = C.A[this->hindex - i - 1]->size;
        }
        else sum[i] = 0;
      });

      pbbslib::scan_add_inplace(sum);

      par_for(0, sum.size(), [&] (size_t i) {
        if (this->hindex - i - 1 < C.size && C.A[this->hindex - i - 1] != nullptr) {
          sum[i] += C.A[this->hindex - i - 1]->size;
        }
      });

      //Create boolean M array
      pbbs::sequence<bool> M = pbbs::sequence<bool>(sum.size() + 1); //True means hindex works, need to find first true

      M[0] = (aboveH >= this->hindex);
      par_for(1, M.size(), [&] (size_t i) {
        M[i] = (aboveH + sum[i - 1] >= this->hindex - i);
      });

      //Find index of first true in M
      pbbs::sequence<size_t> indexM = pbbs::sequence<size_t>(M.size());
      par_for(0, M.size(), [&] (size_t i) {
       indexM[i] = M[i] ? i : UINT_E_MAX;
      });
      uintE hindexDecrease = pbbslib::reduce_min(indexM);
      if (hindexDecrease == UINT_E_MAX) hindexDecrease = sortedBatch.size();

      this->hindex -= hindexDecrease;

      //--------------------------Updating Variables--------------------------//
      //Calculating B = hindex + C[hindex] - aboveH
      if (hindexDecrease == 0) {
        if (this->hindex >= C.size || C.A[this->hindex] == nullptr) bSize = 0;
        else {
          bSize = this->hindex + C.A[this->hindex]->size - aboveH;
        }
      }
      else {
        aboveH += sum[hindexDecrease - 1];

        if (this->hindex >= C.size || C.A[this->hindex] == nullptr) bSize = 0;
        else {
          bSize = this->hindex + C.A[this->hindex]->size - aboveH;
        }
      }

      return this->hindex;
    }
  

    /**
     * Given a batch of vertices and their corresponding degrees, dynamically adds all of the vertices to C
     * Automatically resizes everything by itself
     * 
     * @param batch, sequence of vertices to be added to C
     *     Batch MUST be sorted by degree in nonascending order
     * @param deg, sequence of degrees mapping to the vertices in batch
     *     e.g. vertex at batch[i] must have degree of deg[i]
     */
    void addToC(pbbs::sequence<uintE> batch, pbbs::sequence<uintE> deg) {

      /*    
      //Subraction
      sequence<bool> difference = pbbs::sequence<bool>(batch.size() - 1);
    
      par_for(0, deg.size() - 1, [&] (size_t i) {
        difference[i] = (deg[i] != deg[i + 1]);  
      });
      */
      sequence<size_t> idx = pbbs::sequence<size_t>(batch.size());

      par_for(0, batch.size(), [&] (size_t i) {
        idx[i] = i;
      });

      auto f = [&] (size_t i) { return i == batch.size() - 1 || deg[i] != deg[i + 1]; };
      auto indices = pbbs::filter(idx, f);
      idx.clear();

      //Uses indices sequence to know the index of last entry for each clustered deg
      par_for(0, indices.size(), [&] (size_t i) {
        size_t start = (i == 0 ? 0 : indices[i - 1] + 1);
        size_t end = indices[i] + 1;

        pbbs::sequence<uintE> extra = batch.slice(start, end);

        if (C.A[deg[indices[i]]] == nullptr) {
          C.A[deg[indices[i]]] = new pbbslib::dyn_arr<uintE>(1);
        }
        C.A[deg[indices[i]]]->add(extra);

        extra.clear();
      });
      indices.clear();
    }

    /**
     * Given a batch of vertices and their corresponding degrees, dynamically removes all of the vertices from C
     * Automatically resizes everything by itself
     * 
     * @param batch, sequence of vertices to be removed from C
     *     Batch MUST be sorted by degree in nonascending order
     * @param deg, sequence of degrees mapping to the vertices in batch
     *     e.g. vertex at batch[i] must have degree of deg[i]
     */
    void removeFromC(pbbs::sequence<uintE> batch, pbbs::sequence<uintE> deg) {

      sequence<size_t> idx = pbbs::sequence<size_t>(batch.size());

      par_for(0, batch.size(), [&] (size_t i) {
        idx[i] = i;
      });

      auto f = [&] (size_t i) { return i == batch.size() - 1 || deg[i] != deg[i + 1]; };
      auto indices = pbbs::filter(idx, f);
      //Uses indices sequence to know the index of last entry for each clustered deg
      par_for(0, indices.size(), [&] (size_t i) {
        size_t start = (i == 0 ? 0 : indices[i - 1] + 1);
        size_t end = indices[i] + 1;
        pbbs::sequence<uintE> extra = batch.slice(start, end);
      
        auto remove = [&] (uintE element) {
          for (size_t i = 0; i < extra.size(); i++) {
            if (element == extra[i]) return false;
          }
          return true;
        };

        if (deg[indices[i]] != this->hindex) {

          auto filtered = pbbs::new_array_no_init<uintE>(C.A[deg[indices[i]]]->capacity);
          size_t removed = pbbslib::filter_seq(C.A[deg[indices[i]]]->A, filtered, C.A[deg[indices[i]]]->size, remove);

        
          if (C.A[deg[indices[i]]]->alloc) {
            pbbslib::free_array(C.A[deg[indices[i]]]->A);
          }

          C.A[deg[indices[i]]]->A = filtered;
          C.A[deg[indices[i]]]->size = removed;

        }
        else {
          auto filtered = pbbs::new_array_no_init<uintE>(C.A[deg[indices[i]]]->capacity);

          uintE removedFromB = pbbslib::filter_seq(C.A[deg[indices[i]]]->A, filtered, bSize, remove); 
          size_t removed = pbbslib::filter_seq(&(C.A[deg[indices[i]]]->A)[bSize], &(filtered)[removedFromB], C.A[deg[indices[i]]]->size - bSize, remove);
        
          if (C.A[deg[indices[i]]]->alloc) {
            pbbslib::free_array(C.A[deg[indices[i]]]->A);
          }
          C.A[deg[indices[i]]]->A = filtered;
          C.A[deg[indices[i]]]->size = removed + removedFromB;

          bSize = removedFromB;
        }


        if (C.A[deg[indices[i]]]->size == 0) {
          if (C.A[deg[indices[i]]]->alloc) {
            pbbslib::free_array(C.A[deg[indices[i]]]->A);
          }
          C.A[deg[indices[i]]] = nullptr;
        }      
        
        extra.clear();
      });

      indices.clear();
    }

    /**
     * Helper function that shrinks the size of C so it doesn't take too much memory 
     */
    void adjust() {
      /* Saves memory, more time - Sparse Graphs
      for (; C.size / 2 > 0; C.size--) {
        if (C.A[C.size - 1] != nullptr && C.A[C.size - 1]->size != 0) {
          break;
        }
      }
      */
    
      C.size = std::min(this->G->n, C.size); //More memory, saves time - Dense Graphs
      if (C.size < C.capacity / 4 && 2 * C.size > 1) {
        auto nA = pbbs::new_array_no_init<pbbslib::dyn_arr<uintE>*>(std::max(2 * C.size, (size_t) 1));
        par_for(0, C.size, [&] (size_t i) {
          nA[i] = C.A[i];
        });
        par_for(C.size, std::max(2 * C.size, (size_t) 1), [&] (size_t i) {
          nA[i] = nullptr;
        });
        if (C.alloc) {
         //pbbslib::free_array(C.A);
        }
        C.A = nA;
        C.capacity = std::max(2 * C.size, (size_t) 1);
      }
      //Alternative is to set C.size to min(number of vertices, C.size). Faster but could waste space on really sparse graphs
    }


  public:

    /**
     * Constructs HSetDynArr given a pointer to a dynamic graph
     *
     * @param _G, an unweighted dynamic_symmetric_graph, graph does not have to be empty
     */
    HSetDynArr(dynamic_symmetric_graph<dynamic_symmetric_vertex, pbbs::empty>* _G) : HSet(_G) {

      bSize = 0;

      C = pbbslib::dyn_arr<pbbslib::dyn_arr<uintE>*>(std::max(2 * this->G->n, (size_t) 1));

      par_for(0, C.capacity, [&] (size_t i) { 
        C.A[i] = nullptr;
      });
 
      //Add all starting graph elements into HSet
      if (_G->n > 0) {
        auto start = pbbs::sequence<uintE>(_G->n);
        par_for(0, start.size(), [&] (size_t i) {
          start[i] = i;
        });

        insert(start);
      }

    }

    /**
     * Given a batch of vertices, adds all of the new vertices to HSet in parallel
     *
     * @param vertices, sequence of vertices to be added, can contain existing vertices (will just be ignored)
     * @return the h-index after all the vertex insertions
     */
    uintE insertVertices(sequence<uintE> vertices) {
      this->G->batchAddVertices(vertices);
      insert(vertices);
      adjust();
      return this->hindex;
    }

    /**
     * Given a batch of vertices in the graph, deletes all of the existing ones from HSet in parallel
     *
     * @param vertices, sequence of vertices to be deleted, can contain vertices that don't exist yet (will just be ignored)
     * @return the h-index after all the vertex deletions
     */
    uintE eraseVertices(sequence<uintE> vertices) {
      erase(vertices);
      this->G->batchRemoveVertices(vertices);
      adjust();
      return this->hindex;
    }

    /**
     * Given a batch edges, inserts all of the new edges in parallel
     * Automatically adds any new vertices in the edge list
     *
     * @param edges, sequence of edges to be added
     *      Adding edge u, v also adds edge v, u since the graph is symmetric
     *      CANNOT contain duplicate  edges (use the getEdges() function from SubgraphCounting.cc to make sure)
     *      Can contain edges that already exist (will just be ignored)
     *      Edges can contain new vertices (will be added automatically)
     * @return the h-index after adding all the edges
     */
    uintE insertEdges(sequence<std::pair<uintE, uintE>> edges) {

      auto existEdges = pbbs::filter(edges, [&] (std::pair<uintE, uintE> e) { return !this->G->existEdge(e.first, e.second); } );

      //Get unique vertices
      sequence<uintE> vertices = sequence<uintE>(2 * existEdges.size());
      par_for(0, existEdges.size(), [&] (size_t i) {
       vertices[2 * i] = existEdges[i].first;
       vertices[(2 * i) + 1] = existEdges[i].second; 
      });

      pbbs::sequence<uintE> sortedV = integer_sort(vertices, [&] (uintE v) { return v; });
      vertices.clear();

      auto vTemp = sequence<uintE>(sortedV.size(), [&] (size_t i) {
        if(i == 0) {
          return sortedV[i];
        }
        else {
          if (sortedV[i] != sortedV[i - 1]) {
            return sortedV[i];
          }
          else {
            return UINT_E_MAX;
          }
        }
      });
      sortedV.clear();

      auto uniqueVertices = filter(vTemp, [&] (uintE i) { return (i != UINT_E_MAX); } );
      vTemp.clear();

      //Add vertices that don't exist yet
      pbbs::sequence<uintE> newVertices = filter(uniqueVertices, [&] (uintE v) { 
        if (v < this->G->existVertices.size) return !this->G->existVertices.A[v];
        else return true;
      });

      if (newVertices.size() != 0) insertVertices(newVertices);

      erase(uniqueVertices);
      this->G->batchAddEdges(existEdges);
      insert(uniqueVertices);

      uniqueVertices.clear();
      adjust();
      return this->hindex;
    }

    /**
     * Given a batch edges, deletes all of the new edges in parallel
     * Automatically removes any zero degree vertices from the graph after the deletion
     *
     * @param edges, sequence of edges to be erased
     *      Erasing edge u, v also erases edge v, u since the graph is symmetric
     *      CANNOT contain duplicate edges (use the getEdges() function from SubgraphCounting.cc to make sure)
     *      Can contain edges that don't exist (will be ignored)
     * @return the h-index after deleting all the edges
     */
    uintE eraseEdges(sequence<std::pair<uintE, uintE>> edges) {
    
      auto existEdges = pbbs::filter(edges, [&] (std::pair<uintE, uintE> e) { return this->G->existEdge(e.first, e.second); } );

      //Get unique vertices
      sequence<uintE> vertices = sequence<uintE>(2 * existEdges.size());
      par_for(0, existEdges.size(), [&] (size_t i) {
       vertices[2 * i] = existEdges[i].first;
       vertices[(2 * i) + 1] = existEdges[i].second; 
      });

      pbbs::sequence<uintE> sortedV = integer_sort(vertices, [&] (uintE v) { return v; });
      vertices.clear();

      auto vTemp = sequence<uintE>(sortedV.size(), [&] (size_t i) {
        if(i == 0) {
          return sortedV[i];
        }
        else {
          if (sortedV[i] != sortedV[i - 1]) {
            return sortedV[i];
          }
          else {
            return UINT_E_MAX;
          }
        }
      });
      sortedV.clear();

      auto uniqueVertices = filter(vTemp, [&] (uintE i) { return (i != UINT_E_MAX); } );
      vTemp.clear();

      auto edgesOrdered = pbbs::sequence<std::pair<uintE, uintE>>(existEdges.size(), [&] (size_t i) {
        if (existEdges[i].first > existEdges[i].second) return std::make_pair(existEdges[i].second, existEdges[i].first);
        else return existEdges[i];
      });


      erase(uniqueVertices);
      this->G->batchRemoveEdges(existEdges);
      insert(uniqueVertices);

      auto toRemove = pbbs::filter(uniqueVertices, [&] (uintE v) {
        if (v < this->G->existVertices.size) {
          return this->G->existVertices.A[v] && this->G->get_vertex(v).degree == 0;
        }
        else return false;
      });
      eraseVertices(toRemove);

      uniqueVertices.clear();
      adjust();
      return this->hindex;
    }

    /**
     * Returns a boolean that determines if a vertex belongs to H or not
     *
     * @param target, the vertex in question
           target does not have to be in the graph (function will return false)
     * @return whether or not target is in H
     */
    bool contains(uintE target) {
      if (target >= this->G->existVertices.size || !this->G->existVertices.A[target]) return false;
      uintE deg = this->G->get_vertex(target).degree;

      if (deg > this->hindex) return true;
      else if (deg == this->hindex) {
        for (size_t i = 0; i < bSize; i++) {
          if (C.A[this->hindex]->A[i] == target) return true;
        }
      }
      return false;
    }

    /**
     * Returns a sequence with all of the vertices in H
     *
     * @return sequence of all the vertices in H
     */
    pbbs::sequence<uintE> getH() {
      if (this->hindex == 0) {
        return pbbs::sequence<uintE>(0);
      }

      auto prefixSum = pbbs::sequence<uintE>(C.size - hindex - 1);
      par_for(0, prefixSum.size(), [&] (size_t i) {
        if (C.A[this->hindex + i + 1] != nullptr) {
          prefixSum[i] = C.A[this->hindex + i + 1]->size;
        }
        else prefixSum[i] = 0;
      });

      pbbslib::scan_add_inplace(prefixSum);

      par_for(0, prefixSum.size(), [&] (size_t i) {
        prefixSum[i] += bSize;
      });


      pbbs::sequence<uintE> hSeq = pbbs::sequence<uintE>(hindex);

      //Add everything in B
      if (C.A[this->hindex] != nullptr && bSize != 0) {
        par_for(0, bSize, [&] (size_t i) {
          hSeq[i] = C.A[this->hindex]->A[i];
        });
      }

      //Add rest of vertices
      par_for(0, prefixSum.size(), [&] (size_t i) {
        if (C.A[this->hindex + i + 1] != nullptr) {
          uintE offSet = prefixSum[i];
          par_for(0, C.A[this->hindex + i + 1]->size, [&] (size_t j) {
            hSeq[j + offSet] = C.A[this->hindex + i + 1]->A[j];
          });
        }
      });

      return hSeq;
    }

    /**
     * Frees all data structures used in HSet
     */
    void del() {
      par_for(0, C.size, [&] (size_t i) {
        if (C.A[i] != nullptr) pbbslib::free_array(C.A[i]->A);
      });
      pbbslib::free_array(C.A);
    }

};
