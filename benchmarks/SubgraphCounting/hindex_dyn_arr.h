#pragma once

#include "hindex.h"

//---------------------------------------------------------------------------------//
//HSET DYN_ARR
//---------------------------------------------------------------------------------//

class HSetDynArr : public HSet {

  public:

    uintE bSize; //Number of elements with degree hindex and are inside H
    pbbslib::dyn_arr<pbbslib::dyn_arr<uintE>*> C; //Sorts all vertices by degree

    //Constructor - Takes in initial graph
    HSetDynArr(dynamic_symmetric_graph<dynamic_symmetric_vertex, pbbs::empty>* _G) : HSet(_G) {

      bSize = 0;

      C = pbbslib::dyn_arr<pbbslib::dyn_arr<uintE>*>(std::max(2 * this->G->n, INIT_C_SIZE));

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

    void initExtra(uintE n) {}

    //--------------------------INSERT--------------------------//
    //batch is a sequence of vertex ids (that can be used to get from graph)
    //sorted = true if batch is already sorted in descending order
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

    //--------------------------ERASE--------------------------//   
    //batch is a sequence of vertex ids (that can be used to get from graph)
    //sorted = true if batch is already sorted in descending order
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
  

    //--------------------------ADD TO C--------------------------//
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
          C.A[deg[indices[i]]] = new pbbslib::dyn_arr<uintE>(INIT_DEG_SIZE);
        }
        C.A[deg[indices[i]]]->add(extra);

        extra.clear();
      });
      indices.clear();
    }

    //--------------------------REMOVE FROM C--------------------------//
    void removeFromC(pbbs::sequence<uintE> batch, pbbs::sequence<uintE> deg) {
    
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

    void adjust() {
      /* Saves memory, more time - Sparse Graphs
      for (; C.size / 2 > 0; C.size--) {
        if (C.A[C.size - 1] != nullptr && C.A[C.size - 1]->size != 0) {
          cout << "---" << C.size << endl;
          break;
        }
      }
      */
    
      C.size = std::min(this->G->n, C.size); //More memory, saves time - Dense Graphs
      if (C.size < C.capacity / 4 && 2 * C.size > INIT_C_SIZE) {
        auto nA = pbbs::new_array_no_init<pbbslib::dyn_arr<uintE>*>(std::max(2 * C.size, INIT_C_SIZE));
        par_for(0, C.size, [&] (size_t i) {
          nA[i] = C.A[i];
        });
        par_for(C.size, C.capacity, [&] (size_t i) {
          nA[i] = nullptr;
        });
        if (C.alloc) {
         pbbslib::free_array(C.A);
        }
        C.A = nA;
        C.capacity = std::max(2 * C.size, INIT_C_SIZE);
      }
      //Alternative is to set C.size to min(number of vertices, C.size). Faster but could waste space on really sparse graphs
    }
  

    uintE insertVertices(sequence<uintE> vertices) {
      this->G->batchAddVertices(vertices);
      insert(vertices);
      adjust();
      return this->hindex;
    }

    uintE eraseVertices(sequence<uintE> vertices) {
      this->G->batchRemoveVertices(vertices);
      erase(vertices);
      adjust();
      return this->hindex;
    }

    //Insert edges once e.g. u--v inserts v--u as well
    uintE insertEdges(sequence<std::pair<uintE, uintE>> edges) {

      //Get unique vertices
      sequence<uintE> vertices = sequence<uintE>(2 * edges.size());
      par_for(0, edges.size(), [&] (size_t i) {
       vertices[2 * i] = edges[i].first;
       vertices[(2 * i) + 1] = edges[i].second; 
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
      this->G->batchAddEdges(edges);
      insert(uniqueVertices);

      uniqueVertices.clear();
      adjust();
      return this->hindex;
    }

    uintE eraseEdges(sequence<std::pair<uintE, uintE>> edges) {
    
      //Get unique vertices
      sequence<uintE> vertices = sequence<uintE>(2 * edges.size());
      par_for(0, edges.size(), [&] (size_t i) {
       vertices[2 * i] = edges[i].first;
       vertices[(2 * i) + 1] = edges[i].second; 
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

      auto edgesOrdered = pbbs::sequence<std::pair<uintE, uintE>>(edges.size(), [&] (size_t i) {
        if (edges[i].first > edges[i].second) return std::make_pair(edges[i].second, edges[i].first);
        else return edges[i];
      });


      erase(uniqueVertices);
      this->G->batchRemoveEdges(edges);
      insert(uniqueVertices);

      uniqueVertices.clear();
      adjust();
      return this->hindex;
    }

    //Contains can just check the degree of a vertex and see if it is greater than |H|
    //If vertex's degree is the hindex, we check if it is contained in B
    bool contains(uintE target) {
      uintE deg = this->G->get_vertex(target).degree;

      if (deg > this->hindex) return true;
      else if (deg == this->hindex) {
        for (size_t i = 0; i < bSize; i++) {
          if (C.A[this->hindex]->A[i] == target) return true;
        }
      }
      return false;
    }

    pbbs::sequence<uintE> getH() {
      if (this->hindex == 0) {
        return pbbs::sequence<uintE>(0);
      }

      auto H = make_sparse_table<uintE, pbbs::empty,hash_uintE>(2 * this->hindex + 1 ,std::make_tuple(UINT_E_MAX, pbbs::empty()),hash_uintE());
      par_for(this->hindex + 1, C.size, [&] (size_t i) {
        if (C.A[i] != nullptr) {
          par_for(0, C.A[i]->size, [&] (size_t j) {
            H.insert(std::make_tuple(C.A[i]->A[j], pbbs::empty()));
          });
        }
      });

      if (C.A[this->hindex] != nullptr && bSize != 0) {
        par_for(0, bSize, [&] (size_t i) {
          H.insert(std::make_tuple(C.A[this->hindex]->A[i], pbbs::empty()));
        });
      }

      auto hSeq = H.entries();
      pbbs::sequence<uintE> result = pbbs::sequence<uintE>(this->hindex);

      par_for(0, result.size(), [&] (size_t i) {
        result[i] = std::get<0>(hSeq[i]);
      });

      return result;
    
    }

};
