#pragma once

#include "hindex.h"

/**
 * HSet implementation that dynamically maintains H in parallel
 * Uses a dyn_arr to bucket the vertices
 *
 * H stores h vertices, where h is the largest number such
 * that there are at least h vertices with degree greater than or equal to h
 */
class HSetDynArr : public HSet {

  private:

    uintE bSize; //Number of elements with degree hindex and are inside H
    pbbslib::dyn_arr<pbbslib::dyn_arr<uintE>*> C; //Groups vertices by degree

    /**
     * Given a batch of vertices, adds all of them to HSet
     * Automatically resizes everything by itself
     *
     * @param batch, sequence of vertices to be added to HSet
     *     All vertices in batch must be in the graph
     *     Can be empty
     *     CANNOT add vertex if it is already tracked by HSet (present in C)
     * @param sorted, optional boolean
     *     Is batch already sorted by degree in nonascending order
     *     Set to false by default
     *     Saves time by not resorting batch if it's already sorted
     */
    uintE insert(sequence<uintE> batch, bool sorted = false) {

      //Exits if there are no vertices in batch
      if (batch.size() == 0) return hindex;

      //Sorts vertices by degree in nonascending order (if necessary)
      pbbs::sequence<uintE> sortedBatch;
      if (!sorted) {
        sortedBatch = integer_sort(batch, [&] (uintE v) {
          return G->get_vertex(v).degree;
        }).rslice();
      }
      else {
        sortedBatch = batch;
      }
      batch.clear();

      //Stores degrees of vertices in nonascending order
      sequence<uintE> deg = pbbs::sequence<uintE>(sortedBatch.size());
      par_for(0, sortedBatch.size(), [&] (size_t i) {
        deg[i] = deg[i] = this->G->get_vertex(sortedBatch[i]).degree;
      });

      //Adds additional space to C if it is too small
      //Makes sure it has room for the largest degree
      size_t oldSize = C.size;
      C.size = std::max((uintE) C.size, deg[0] + 1);
      if (C.size >= C.capacity) {
        C.resize(0);
        par_for(oldSize, C.capacity, [&] (size_t i) {
          C.A[i] = nullptr;
        });
      }

      //----------------------Updating Variables----------------------//
      //Current number of vertices with degree >= h is h - |B| + |C[h]|
      uintE aboveH = this->hindex - bSize;
      if (this->hindex < C.size && C.A[this->hindex] != nullptr) {
        aboveH += C.A[this->hindex]->size;
      }

      //Number of vertices being added in batch with degree >= current h
      uintE added = pbbs::binary_search(deg, [&] (uintE x) {
        return x >= this->hindex;
      });

      aboveH += added;

      addToC(sortedBatch, deg); //Updates C with new vertices

      //----------------------Compute H-index----------------------//
      //Begins constructing sum array
      //Indicates the number of vertices H loses if hindex increases

      //Stores the sizes of each entry of C
      pbbs::sequence<uintE> sum;
      sum = pbbs::sequence<uintE>(sortedBatch.size() + 1);

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

      //Creates the binary array M
      //True indicates that h can lose those vertices and still be valid
      pbbs::sequence<bool> M = pbbs::sequence<bool>(sum.size());
    
      //Still valid h when it increases by i + 1
      par_for(0, M.size(), [&] (size_t i) {
        M[i] = aboveH - sum[i] >= this->hindex + i + 1;
      });

      //Find index of first false
      pbbs::sequence<size_t> indexM = pbbs::sequence<size_t>(M.size());
      par_for(0, M.size(), [&] (size_t i) {
        indexM[i] = !M[i] ? i : UINT_E_MAX;
      });
      uintE hindexIncrease = pbbslib::reduce_min(indexM);
      if (hindexIncrease == UINT_E_MAX) hindexIncrease = sortedBatch.size();

      //Update h
      this->hindex += hindexIncrease;

      //----------------------Updating Variables----------------------//

      //If h does not increase, aboveH does not change
      if (hindexIncrease == 0) {
        //If C[h] is empty, B must also be empty
        if (this->hindex >= C.size || C.A[this->hindex] == nullptr) bSize = 0;
        else {
          //Calculate |B| = hindex + C[h] - aboveH
          bSize = this->hindex + C.A[this->hindex]->size - aboveH;
        }
      }
      else {
        //Update aboveH because h just increased and lost some vertices
        aboveH -= sum[hindexIncrease - 1];

        //If C[h] is empty, B must also be empty
        if (this->hindex >= C.size || C.A[this->hindex] == nullptr) bSize = 0;
        else {
          //Calculate |B| = hindex + C[h] - aboveH
          bSize = this->hindex + C.A[this->hindex]->size - aboveH;
        }
      }

      //Delete unused arrays before exiting
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
     *     CANNOT remove a vertex if it is not tracked by HSet
     * @param sorted, optional boolean
     *     Is batch already sorted by degree in nonascending order
     *     Set to false by default
     *     Saves time by not resorting batch if it's already sorted
     */
    uintE erase(sequence<uintE> batch, bool sorted = false) {

      //Sorts vertices by degree in nonascending order (if needed) 
      pbbs::sequence<uintE> sortedBatch;
      if (!sorted) {
        sortedBatch = integer_sort(batch, [&] (uintE v) {
          return this->G->get_vertex(v).degree;
        }).rslice();
      }
      else {
        sortedBatch = batch;
      } 
      batch.clear();

      //Stores degrees of vertices in nonascending order
      sequence<uintE> deg = pbbs::sequence<uintE>(sortedBatch.size());
      par_for(0, sortedBatch.size(), [&] (size_t i) {
        deg[i] = this->G->get_vertex(sortedBatch[i]).degree;
      });

      //----------------------Updating Variables----------------------//
      //Current number of vertices with degree >= h is h - |B| + |C[h]|
      uintE aboveH = this->hindex - bSize;

      if (this->hindex < C.size && C.A[this->hindex] != nullptr) {
        aboveH += C.A[this->hindex]->size;
      }

      //Number of vertices being added in batch with degree >= current h
      uintE removed = pbbs::binary_search(deg, [&] (uintE x) {
        return x >= this->hindex;
      });

      aboveH -= removed;

      removeFromC(sortedBatch, deg); //Updates C by removing vertices


      //----------------------Compute H-index----------------------//
      //Begins constructing sum array
      //Each number indicates the number of vertices H gains if h decreases

      //Stores the sizes of each entry of C
      //Makes sure sum array does not store sizes for negative degrees
      size_t size = std::min(sortedBatch.size(), this->hindex);
      pbbs::sequence<uintE> sum = pbbs::sequence<uintE>(size);
      par_for(0, sum.size(), [&] (size_t i) {
        if (this->hindex - i - 1 < C.size && 
            C.A[this->hindex - i - 1] != nullptr) { 

          sum[i] = C.A[this->hindex - i - 1]->size;
        }
        else sum[i] = 0;
      });

      //Computes Prefix Sum (note scan_add_inplace is exclusive)
      pbbslib::scan_add_inplace(sum);
      par_for(0, sum.size(), [&] (size_t i) {
        if (this->hindex - i - 1 < C.size &&
            C.A[this->hindex - i - 1] != nullptr) {

          sum[i] += C.A[this->hindex - i - 1]->size; //Compute inclusive sum
        }
      });

      //Create boolean M array
      //True indicates that h is valid after gaining those vertices
      pbbs::sequence<bool> M = pbbs::sequence<bool>(sum.size() + 1);

      //Is h valid after decreasing h by i
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

      //Update h
      this->hindex -= hindexDecrease;

      //----------------------Updating Variables----------------------//
      //If h does not increase, aboveH does not change
      if (hindexDecrease == 0) {
        //If C[h] is empty, B must also be empty
        if (this->hindex >= C.size || C.A[this->hindex] == nullptr) bSize = 0;
        else {
          //Calculating B = hindex + C[hindex] - aboveH
          bSize = this->hindex + C.A[this->hindex]->size - aboveH;
        }
      }
      else {
        //Update aboveH because h just increased and lost some vertices
        aboveH += sum[hindexDecrease - 1];

        //If C[h] is empty, B must also be empty
        if (this->hindex >= C.size || C.A[this->hindex] == nullptr) bSize = 0;
        else {
          //Calculate |B| = hindex + C[h] - aboveH
          bSize = this->hindex + C.A[this->hindex]->size - aboveH;
        }
      }

      return this->hindex;
    }
  

    /**
     * Given a batch of vertices and their corresponding degrees
     *   dynamically adds all of the vertices to C
     *
     * Automatically resizes everything by itself
     * 
     * @param batch, sequence of vertices to be added to C
     *     Batch MUST be sorted by degree in nonascending order
     * @param deg, sequence of degrees mapping to the vertices in batch
     *     e.g. vertex at batch[i] must have degree of deg[i]
     */
    void addToC(pbbs::sequence<uintE> batch, pbbs::sequence<uintE> deg) {

      //Stores array indices
      //Gives the ending index of a block of vertices with the same degree

      //Stores all indices (0, 1, ..., |batch|) in idx
      sequence<size_t> idx = pbbs::sequence<size_t>(batch.size());
      par_for(0, batch.size(), [&] (size_t i) {
        idx[i] = i;
      });

      //Filters so there only indices with degree different than the next 
      auto f = [&] (size_t i) { 
        return i == batch.size() - 1 || deg[i] != deg[i + 1];
      };
      auto indices = pbbs::filter(idx, f);
      idx.clear();

      //Uses indices sequence to know the index range of each block of
      //  of vertices with the same degree
      par_for(0, indices.size(), [&] (size_t i) {
        size_t start = (i == 0 ? 0 : indices[i - 1] + 1);
        size_t end = indices[i] + 1;

        //extra stores the vertices that should be added
        pbbs::sequence<uintE> extra = batch.slice(start, end);

        //If there isn't an entry for this degree, create one
        if (C.A[deg[indices[i]]] == nullptr) {
          C.A[deg[indices[i]]] = new pbbslib::dyn_arr<uintE>(1);
        }
        //Add all the new vertices
        C.A[deg[indices[i]]]->add(extra);

        extra.clear();
      });
      indices.clear();
    }

    /**
     * Given a batch of vertices and their corresponding degrees
     *   dynamically removes all of the vertices from C
     *
     * Automatically resizes everything by itself
     * 
     * @param batch, sequence of vertices to be removed from C
     *     Batch MUST be sorted by degree in nonascending order
     * @param deg, sequence of degrees mapping to the vertices in batch
     *     e.g. vertex at batch[i] must have degree of deg[i]
     */
    void removeFromC(pbbs::sequence<uintE> batch, pbbs::sequence<uintE> deg) {

      //Stores array indices
      //Gives the ending index of a block of vertices with the same degree

      //Stores all indices (0, 1, ..., |batch|) in idx
      sequence<size_t> idx = pbbs::sequence<size_t>(batch.size());

      par_for(0, batch.size(), [&] (size_t i) {
        idx[i] = i;
      });

      //Filters so there are only indices with degree different than the next
      auto f = [&] (size_t i) {
        return i == batch.size() - 1 || deg[i] != deg[i + 1];
      };

      auto indices = pbbs::filter(idx, f);
      //Uses indices sequence to know the index range of each block of
      //  of vertices with the same degree
      par_for(0, indices.size(), [&] (size_t i) {
        size_t start = (i == 0 ? 0 : indices[i - 1] + 1);
        size_t end = indices[i] + 1;
        pbbs::sequence<uintE> extra = batch.slice(start, end);
      
        //Filter out all vertices to be deleted
        auto remove = [&] (uintE element) {
          for (size_t i = 0; i < extra.size(); i++) {
            if (element == extra[i]) return false;
          }
          return true;
        };

        //Removing these vertices does not affect B
        if (deg[indices[i]] != this->hindex) {

          //Filter out elements in the entry of C
          auto filtered = pbbs::new_array_no_init<uintE>(
            C.A[deg[indices[i]]]->capacity
          );

          size_t removed = pbbslib::filter_seq(
            C.A[deg[indices[i]]]->A, 
            filtered, C.A[deg[indices[i]]]->size,
            remove
          );

          //Free the old entry
          if (C.A[deg[indices[i]]]->alloc) {
            pbbslib::free_array(C.A[deg[indices[i]]]->A);
          }

          //Update the new entry and its sizes
          C.A[deg[indices[i]]]->A = filtered;
          C.A[deg[indices[i]]]->size = removed;

        }

        //Removing these vertices possibly affects B
        else {
          auto filtered = pbbs::new_array_no_init<uintE>(
            C.A[deg[indices[i]]]->capacity
          );

          //Filter B, the first |B| elements of the entry
          uintE removedFromB = pbbslib::filter_seq(
            C.A[deg[indices[i]]]->A, filtered, bSize,
            remove
          );

          //Filter the rest of the entry
          size_t removed = pbbslib::filter_seq(
            &(C.A[deg[indices[i]]]->A)[bSize],
            &(filtered)[removedFromB],
            C.A[deg[indices[i]]]->size - bSize,
            remove
          );
        
          //Free the old entry
          if (C.A[deg[indices[i]]]->alloc) {
            pbbslib::free_array(C.A[deg[indices[i]]]->A);
          }
     
          //Update the new entry and its sizes (including the size of B)
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
     * Helper function that shrinks the size of C
     *   so it doesn't take too much memory 
     */
    void adjust() {
      /* Saves memory, uses more time - better for sparse graphs

      //Finds the largest degree in the graph and sets that to be the size
      for (; C.size / 2 > 0; C.size--) {
        if (C.A[C.size - 1] != nullptr && C.A[C.size - 1]->size != 0) {
          break;
        }
      }
      */

      //Uses more memory, but saves time - better for dense graphs

      //Assumes size (the largest degree) is at most the number of vertices
      C.size = std::min(this->G->n, C.size);

      //If it takes up less than the capacity
      if (C.size < C.capacity / 4 && 2 * C.size > 1) {
        auto nA = pbbs::new_array_no_init<pbbslib::dyn_arr<uintE>*>(
          std::max(2 * C.size, (size_t) 1)
        );

        par_for(0, C.size, [&] (size_t i) {
          nA[i] = C.A[i];
        });
        par_for(C.size, std::max(2 * C.size, (size_t) 1), [&] (size_t i) {
          nA[i] = nullptr;
        });
        if (C.alloc) {
          pbbslib::free_array(C.A);
        }
        C.A = nA;
        C.capacity = std::max(2 * C.size, (size_t) 1);
      }
    }

//-----------------------------------------------------------------------------
  public:

    /**
     * Constructs HSetDynArr given a pointer to a dynamic graph
     *
     * @param _G, an unweighted dynamic_symmetric_graph
     *     graph does not have to be empty
     */
    HSetDynArr(dynamic_symmetric_graph<dynamic_symmetric_vertex, pbbs::empty>* _G) : //TODO: SPLIT THE LINE

      HSet(_G) {

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
     * Adds all of the new vertices in a given batch to HSet in parallel
     *
     * @param vertices, sequence of vertices to be added
     *     can contain existing vertices (will be ignored)
     *     CANNOT contain duplicate vertices in batch
     * @return the h-index after all the vertex insertions
     */
    uintE insertVertices(sequence<uintE> vertices) {
      this->G->batchAddVertices(vertices);
      insert(vertices);
      adjust();
      return this->hindex;
    }

    /**
     * Deletes all of the existing vertices in batch from HSet in parallel
     *
     * @param vertices, sequence of vertices to be deleted
     *     can contain vertices that don't exist yet (will be ignored)
     *     CANNOT contain duplicate vertices in batch
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
     *      CANNOT contain duplicate  edges
     *          Filter with the getEdges() function from SubgraphCounting.cc
     *      Can contain edges that already exist (will just be ignored)
     *      Edges can contain new vertices (will be added automatically)
     * @return the h-index after adding all the edges
     */
    uintE insertEdges(sequence<std::pair<uintE, uintE>> edges) {

      //Filters edges so it only adds ones that don't exist yet
      auto existEdges = pbbs::filter(edges, [&] (std::pair<uintE, uintE> e) {
        return !this->G->existEdge(e.first, e.second);
      });

      //Gets all vertices of all edges (contains duplicates)
      sequence<uintE> vertices = sequence<uintE>(2 * existEdges.size());
      par_for(0, existEdges.size(), [&] (size_t i) {
       vertices[2 * i] = existEdges[i].first;
       vertices[(2 * i) + 1] = existEdges[i].second; 
      });

      //Sorts vertices
      pbbs::sequence<uintE> sortedV = integer_sort(vertices, [&] (uintE v) {
        return v;
      });
      vertices.clear();

      //Filters so there are only unique vertices
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

      auto uniqueVertices = filter(vTemp, [&] (uintE i) {
        return (i != UINT_E_MAX);
      });
      vTemp.clear();

      //Find vertices that don't exist yet
      pbbs::sequence<uintE> newVertices = filter(uniqueVertices, [&] (uintE v) { 
        if (v < this->G->existVertices.size) return !this->G->existVertices.A[v];
        else return true;
      });

      //Adds them to HSet and graph
      if (newVertices.size() != 0) insertVertices(newVertices);

      erase(uniqueVertices); //Removes from HSet
      this->G->batchAddEdges(existEdges); //Updates graph
      insert(uniqueVertices); //Re-adds back into HSet

      uniqueVertices.clear();
      adjust();
      return this->hindex;
    }

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
    uintE eraseEdges(sequence<std::pair<uintE, uintE>> edges) {
    
      //Filters edges so it only adds ones that don't exist yet
      auto existEdges = pbbs::filter(edges, [&] (std::pair<uintE, uintE> e) {
        return this->G->existEdge(e.first, e.second);
      });

      //Gets all vertices of all edges (contains duplicates)
      sequence<uintE> vertices = sequence<uintE>(2 * existEdges.size());
      par_for(0, existEdges.size(), [&] (size_t i) {
       vertices[2 * i] = existEdges[i].first;
       vertices[(2 * i) + 1] = existEdges[i].second; 
      });

      //Sorts vertices
      pbbs::sequence<uintE> sortedV = integer_sort(vertices, [&] (uintE v) {
        return v;
      });
      vertices.clear();

      //Filters so there are only unique vertices
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

      auto uniqueVertices = filter(vTemp, [&] (uintE i) {
        return (i != UINT_E_MAX);
      });
      vTemp.clear();

      erase(uniqueVertices); //Removes from HSet
      this->G->batchRemoveEdges(existEdges); //Updates edges
      insert(uniqueVertices); //Re-adds to HSet

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
     * @return whether or not target is in H
     */
    bool contains(uintE target) {
      //False if vertex does not exist
      if (target >= this->G->existVertices.size || !this->G->existVertices.A[target]) {
        return false;
      }

      uintE deg = this->G->get_vertex(target).degree;

      //If degree > h, it is definitely in H
      if (deg > this->hindex) return true;

      //If degree == h, check that it is in B
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
      //Return empty sequence if h is 0
      if (this->hindex == 0) {
        return pbbs::sequence<uintE>(0);
      }

      //Find the size of each entry in C
      auto prefixSum = pbbs::sequence<uintE>(C.size - hindex - 1);
      par_for(0, prefixSum.size(), [&] (size_t i) {
        if (C.A[this->hindex + i + 1] != nullptr) {
          prefixSum[i] = C.A[this->hindex + i + 1]->size;
        }
        else prefixSum[i] = 0;
      });

      //Take exclusive prefix sum to find where each block of
      //  same degree vertices start 
      pbbslib::scan_add_inplace(prefixSum);

      //Offset everything by |B| since that will come first
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

      //Add rest of vertices (using the prefix sum to find where they go)
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
     * Returns a sequence with all of the vertices in H
     *
     * @return sequence of all the vertices in H
     */
    pbbs::sequence<uintE> getP() {

      cout << hindex << " " << C.size << endl;

      //Return empty sequence if there are no vertices above 2 * h
      if (this->hindex == 0 || C.size < 2 * hindex) {
        return pbbs::sequence<uintE>(0);
      }

      //Find the size of each entry in C
      auto prefixSum = pbbs::sequence<uintE>(C.size - 2 * hindex);
      par_for(0, prefixSum.size(), [&] (size_t i) {
        if (C.A[2 * this->hindex + i] != nullptr) {
          prefixSum[i] = C.A[2 * this->hindex + i]->size;
        }
        else prefixSum[i] = 0;
      });

      for (int i = 0; i < prefixSum.size(); i++)  cout << prefixSum[i] << endl;

      //Take exclusive prefix sum to find where each block of
      //  same degree vertices start

      uintE total = prefixSum[prefixSum.size() - 1];
      pbbslib::scan_add_inplace(prefixSum);

      total += prefixSum[prefixSum.size() - 1];
      pbbs::sequence<uintE> pSeq = pbbs::sequence<uintE>(total); //Should be the last size of prefix sum
cout << "A" << endl;
      //Add rest of vertices (using the prefix sum to find where they go)
      //par_for(0, prefixSum.size(), [&] (size_t i) {
for (int i = 0; i < prefixSum.size(); i++) {
        if (2 * hindex + 1 < C.size && C.A[2 * this->hindex + i] != nullptr) {
          uintE offSet = prefixSum[i];
          //par_for(0, C.A[2 * this->hindex + i]->size, [&] (size_t j) {
for (int j = 0; j < C.A[2 * this->hindex + i]->size; j++) {
            pSeq[j + offSet] = C.A[this->hindex + i]->A[j];
          }//);
        }
      }//);
cout << "A" << endl;
      return pSeq;
    }

    /**
     * Frees all data structures used in HSet
     */
    void del() {
      //Frees each entry in C
      par_for(0, C.size, [&] (size_t i) {
        if (C.A[i] != nullptr) pbbslib::free_array(C.A[i]->A);
      });

      //Frees entire C array
      pbbslib::free_array(C.A);
    }

};
