#pragma once

#include "hindex.h"

//--------------------------------------------------------------------------//
//HSET THRESHOLD/HYBRID
//--------------------------------------------------------------------------//

class HSetThreshold : public HSet {

  public:
    uintE bSize;
    pbbslib::dyn_arr<pbbslib::dyn_arr<uintE>*> lowC;
    sparse_table<uintE, pbbslib::dyn_arr<uintE>*, hash_uintE> highC;

    uintE threshold;

    uintE highCStored;

    //Constructor
    //HSet(Graph& _G) { //Don't pass by reference
    HSetThreshold(dynamic_symmetric_graph<dynamic_symmetric_vertex, pbbs::empty>* _G, uintE a) : HSet(_G) {

      bSize = 0;

      threshold = a;

      highCStored = 0;

      lowC = pbbslib::dyn_arr<pbbslib::dyn_arr<uintE>*>(threshold);
      highC = make_sparse_table<uintE, pbbslib::dyn_arr<uintE>*, hash_uintE>
        (std::max(2 * _G->n, INIT_C_SIZE), std::make_tuple(UINT_E_MAX, nullptr),hash_uintE());

      par_for(0, lowC.capacity, [&] (size_t i) { 
        lowC.A[i] = nullptr;
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

    //--------------------------INSERT--------------------------//
    //batch is a sequence of vertex ids (that can be used to get from graph)
    //sorted = true if batch is already sorted in descending order
    uintE insert(sequence<uintE> batch, bool sorted = false) {

      if (batch.size() == 0) return this->hindex;

      pbbs::sequence<uintE> sortedBatch;
      if (!sorted) {
        sortedBatch = integer_sort(batch, [&] (uintE v) { return this->G->get_vertex(v).degree; }).rslice();
      }

      else {
        sortedBatch = batch;
      }
      batch.clear();

      //Stores degrees and sorts in descending order
      sequence<uintE> deg = pbbs::sequence<uintE>(sortedBatch.size());
      par_for(0, sortedBatch.size(), [&] (size_t i) {
        deg[i] = deg[i] = this->G->get_vertex(sortedBatch[i]).degree;
      });

      //Adds additional space to C if it is too small
      /*
      size_t highestLowC = pbbs::binary_search(deg, [&] (uintE x) {
        return x >= threshold;
      });
      if (deg[highestLowC] <= deg.size()) {
        lowC.size = std::min(std::max((uintE) lowC.size, deg[highestLowC] + 1), threshold);
      }
      */

      auto needNewEntry = pbbs::sequence<uintE>(deg.size());
      par_for(0, deg.size(), [&] (size_t i) {
        if ((i == 0 || deg[i] != deg[i - 1]) && deg[i] >= threshold) {
          if (highC.contains(deg[i])) needNewEntry[i] = 0;
          else needNewEntry[i] = 1;
        }
        else needNewEntry[i] = 0;
      });

      highCStored += pbbslib::reduce_add(needNewEntry);

      if (highCStored >= highC.m) {

        auto entries = pbbs::filter(highC.entries(), [&] (std::tuple<uintE, pbbslib::dyn_arr<uintE>*> element) {
          return std::get<1>(element) != nullptr;
        });

        if (highC.alloc) pbbslib::free_array(highC.table);
        
        highC = make_sparse_table<uintE, pbbslib::dyn_arr<uintE>*, hash_uintE>
          (std::max((size_t) 2 * highCStored, INIT_C_SIZE), std::make_tuple(UINT_E_MAX, nullptr),hash_uintE());

        par_for(0, entries.size(), [&] (size_t i) {
          highC.insert(entries[i]);
        });

        highCStored = entries.size();

        needNewEntry = pbbs::sequence<uintE>(deg.size());
        par_for(0, deg.size(), [&] (size_t i) {
          if ((i == 0 || deg[i] != deg[i - 1]) && deg[i] >= threshold) {
            if (highC.contains(deg[i])) needNewEntry[i] = 0;
            else needNewEntry[i] = 1;
          }
          else needNewEntry[i] = 0;
        });
        highCStored += pbbslib::reduce_add(needNewEntry);
      }

      //--------------------------Updating Variables--------------------------//
      uintE aboveH = this->hindex - bSize;
      if (getC(this->hindex) != nullptr) aboveH += getC(this->hindex)->size;

      uintE added = pbbs::binary_search(deg, [&] (uintE x) { //Number of vertices that are above new hindex
        return x >= this->hindex;
      });

      aboveH += added;
      
      addToC(sortedBatch, deg);

      //--------------------------Compute H-index--------------------------//
      pbbs::sequence<uintE> sum = pbbs::sequence<uintE>(sortedBatch.size() + 1);
      par_for(0, sortedBatch.size(), [&] (size_t i) {
        if (getC(this->hindex + i) != nullptr) {
          sum[i] = getC(this->hindex + i)->size;

        }
        else sum[i] = 0;
      });

      //Computes Prefix Sum (note scan_add_inplace is exclusive)
      pbbslib::scan_add_inplace(sum);
      par_for(0, sum.size(), [&] (size_t i) {
        if (getC(this->hindex + i) != nullptr) {
          sum[i] += getC(this->hindex + i)->size; //Compute inclusive sum
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
        if (getC(this->hindex) == nullptr) bSize = 0;
        else {
          bSize = this->hindex + getC(this->hindex)->size - aboveH;
        }
      }
      else {
        aboveH -= sum[hindexIncrease - 1];

        if (getC(this->hindex) == nullptr) bSize = 0;
        else {
          bSize = this->hindex + getC(this->hindex)->size - aboveH;
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

      if (getC(this->hindex) != nullptr) aboveH += getC(this->hindex)->size;

      uintE removed = pbbs::binary_search(deg, [&] (uintE x) { //Number of vertices that are above new hindex
        return x >= this->hindex;
      });

      aboveH -= removed;

      removeFromC(sortedBatch, deg);

      //--------------------------Compute H-index--------------------------//
      //Compute prefix sum
      pbbs::sequence<uintE> sum = pbbs::sequence<uintE>(std::min(sortedBatch.size(), this->hindex));
      par_for(0, sum.size(), [&] (size_t i) {
        if (getC(this->hindex - i - 1) != nullptr) { 
          sum[i] = getC(this->hindex - i - 1)->size;
        }
        else sum[i] = 0;
      });

      pbbslib::scan_add_inplace(sum);

      par_for(0, sum.size(), [&] (size_t i) {
        if (getC(this->hindex - i - 1) != nullptr) {
          sum[i] += getC(this->hindex - i - 1)->size;
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
        if (getC(this->hindex) == nullptr) bSize = 0;
        else {
          bSize = this->hindex + getC(this->hindex)->size - aboveH;
        }
      }
      else {
        aboveH += sum[hindexDecrease - 1];

        if (getC(this->hindex) == nullptr) bSize = 0;
        else {
          bSize = this->hindex + getC(this->hindex)->size - aboveH;
        }
      }

      return this->hindex;
    }


    //--------------------------ADD TO C--------------------------//
    void addToC(pbbs::sequence<uintE> batch, pbbs::sequence<uintE> deg) {

      
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

        if (getC(deg[indices[i]]) == nullptr) {
          if (deg[indices[i]] >= threshold) {
            if (highC.contains(deg[indices[i]])) {
              highC.change(deg[indices[i]], new pbbslib::dyn_arr<uintE>(INIT_DEG_SIZE));
            }
            else {
              highC.insert(std::make_tuple(deg[indices[i]], new pbbslib::dyn_arr<uintE>(INIT_DEG_SIZE)));
            }
          }
          else {
            lowC.A[deg[indices[i]]] = new pbbslib::dyn_arr<uintE>(INIT_DEG_SIZE);
          }     
        }

        getC(deg[indices[i]])->add(extra);

        //Make sure not wasting space
        if (getC(deg[indices[i]])->size <= getC(deg[indices[i]])->capacity / 4 && 2 * getC(deg[indices[i]])->size > INIT_DEG_SIZE) {

          auto nA = pbbs::new_array_no_init<uintE>(std::max(2 * getC(deg[indices[i]])->size, INIT_DEG_SIZE));
          par_for(0, getC(deg[indices[i]])->size, [&] (size_t j) {
            nA[j] = getC(deg[indices[i]])->A[j];
          });
   
        
          if (getC(deg[indices[i]])->alloc) {
            pbbslib::free_array(getC(deg[indices[i]])->A);
          }

          getC(deg[indices[i]])->A = nA;
          getC(deg[indices[i]])->capacity = 2 * getC(deg[indices[i]])->size;
          
        }

        extra.clear();
      });
      indices.clear();
    }

    //--------------------------REMOVE FROM C--------------------------//
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

          auto filtered = pbbs::new_array_no_init<uintE>(getC(deg[indices[i]])->capacity);
          size_t removed = pbbslib::filter_seq(getC(deg[indices[i]])->A, filtered, getC(deg[indices[i]])->size, remove);

        
          if (getC(deg[indices[i]])->alloc) {
            pbbslib::free_array(getC(deg[indices[i]])->A);
          }

          getC(deg[indices[i]])->A = filtered;
          getC(deg[indices[i]])->size = removed;

        }
        else {
          auto filtered = pbbs::new_array_no_init<uintE>(getC(deg[indices[i]])->capacity);

          uintE removedFromB = pbbslib::filter_seq(getC(deg[indices[i]])->A, filtered, bSize, remove); 
          size_t removed = pbbslib::filter_seq(&(getC(deg[indices[i]])->A)[bSize], &(filtered)[removedFromB], getC(deg[indices[i]])->size - bSize, remove);
        
          if (getC(deg[indices[i]])->alloc) {
            pbbslib::free_array(getC(deg[indices[i]])->A);
          }
          getC(deg[indices[i]])->A = filtered;
          getC(deg[indices[i]])->size = removed + removedFromB;

          bSize = removedFromB;
        }


        if (getC(deg[indices[i]])->size == 0) {
          if (getC(deg[indices[i]])->alloc) {
            pbbslib::free_array(getC(deg[indices[i]])->A);
          }
          
          if (deg[indices[i]] >= threshold) highC.change(deg[indices[i]], nullptr);
          else lowC.A[deg[indices[i]]] = nullptr;
        }      
        
        extra.clear();
      });

      indices.clear();
    }
  

    pbbslib::dyn_arr<uintE>* getC(uintE deg) {
      if (deg >= threshold) {
        return highC.find(deg, nullptr);
      }
      else return lowC.A[deg];
    }


    void adjust() {
      if (highCStored <= highC.m / 4 && 2 * highCStored >= INIT_C_SIZE) {

        auto entries = pbbs::filter(highC.entries(), [&] (std::tuple<uintE, pbbslib::dyn_arr<uintE>*> element) {
          return std::get<1>(element) != nullptr;
        });

        if (highC.alloc) {
          pbbslib::free_array(highC.table);
        }
        highC = make_sparse_table<uintE, pbbslib::dyn_arr<uintE>*, hash_uintE>
          (2 * highCStored, std::make_tuple(UINT_E_MAX, nullptr),hash_uintE());

        par_for(0, entries.size(), [&] (size_t i) {
          highC.insert(entries[i]);
        });

        highCStored = entries.size();
        
      }
    }



    bool contains(uintE target) {
      uintE deg = this->G->get_vertex(target).degree;

      if (deg > this->hindex) return true;
      else if (deg == this->hindex) {
        for (size_t i = 0; i < bSize; i++) {
          if (getC(deg)->A[i] == target) return true;
        }
      }
      return false;
    }

    uintE insertVertices(pbbs::sequence<uintE> vertices) { 
      this->G->batchAddVertices(vertices);
      insert(vertices);
      adjust();
      return this->hindex;
    }

    uintE eraseVertices(pbbs::sequence<uintE> vertices) { 
      this->G->batchRemoveVertices(vertices);
      erase(vertices);
      adjust();
      return this->hindex;
    }
    uintE insertEdges(pbbs::sequence<std::pair<uintE, uintE>> edges) { 

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

    uintE eraseEdges(pbbs::sequence<std::pair<uintE, uintE>> edges) { 

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

      erase(uniqueVertices);
      this->G->batchRemoveEdges(edges);
      insert(uniqueVertices);

      uniqueVertices.clear();
      adjust();
      return this->hindex;
    }

    
    pbbs::sequence<uintE> getH() {
      if (this->hindex == 0) {
        return pbbs::sequence<uintE>(0);
      }

      auto prefixSum = pbbs::sequence<uintE>(this->G->n - hindex - 1);
      par_for(0, prefixSum.size(), [&] (size_t i) {
        if (getC(this->hindex + i + 1) != nullptr) {
          prefixSum[i] = getC(this->hindex + i + 1)->size;
        }
        else prefixSum[i] = 0;
      });

      pbbslib::scan_add_inplace(prefixSum);

      par_for(0, prefixSum.size(), [&] (size_t i) {
        prefixSum[i] += bSize;
      });


      pbbs::sequence<uintE> hSeq = pbbs::sequence<uintE>(hindex);

      //Add everything in B
      if (getC(this->hindex) != nullptr && bSize != 0) {
        par_for(0, bSize, [&] (size_t i) {
          hSeq[i] = getC(this->hindex)->A[i];
        });
      }

      //Add rest of vertices
      par_for(0, this->G->n, [&] (size_t i) {
        if (getC(this->hindex + i + 1) != nullptr) {
          uintE offSet = prefixSum[i];
          par_for(0, getC(this->hindex + i + 1)->size, [&] (size_t j) {
            hSeq[j + offSet] = getC(this->hindex + i + 1)->A[j];
          });
        }
      });

      return hSeq;

      /*
      if (this->hindex == 0) {
        return pbbs::sequence<uintE>(0);
      }

      auto H = make_sparse_table<uintE, pbbs::empty,hash_uintE>(2 * this->hindex + 1, std::make_tuple(UINT_E_MAX, pbbs::empty()),hash_uintE());

      if (threshold > 0 && this->hindex + 1 < threshold) {
        par_for(this->hindex + 1, threshold, [&] (size_t i) {
          if (lowC.A[i] != nullptr) {
            par_for(0, lowC.A[i]->size, [&] (size_t j) {
              H.insert(std::make_tuple(lowC.A[i]->A[j], pbbs::empty()));
            });
          }
        });
      }

      auto highCEntries = highC.entries();
      par_for(0, highCEntries.size(), [&] (size_t i) {
        if (std::get<0>(highCEntries[i]) > hindex && std::get<1>(highCEntries[i]) != nullptr) {
          par_for(0, std::get<1>(highCEntries[i])->size, [&] (size_t j) {
            H.insert(std::make_tuple(std::get<1>(highCEntries[i])->A[j], pbbs::empty()));
          });
        }
      });

      if (bSize != 0 && getC(this->hindex) != nullptr) {
        par_for(0, bSize, [&] (size_t i) {
          H.insert(std::make_tuple(getC(this->hindex)->A[i], pbbs::empty()));
        });
      }

      auto hSeq = H.entries();
      pbbs::sequence<uintE> result = pbbs::sequence<uintE>(this->hindex);

      par_for(0, result.size(), [&] (size_t i) {
        result[i] = std::get<0>(hSeq[i]);
      });

      return result;
      */
    }
    

};
