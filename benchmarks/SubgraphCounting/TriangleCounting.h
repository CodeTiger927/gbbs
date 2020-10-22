
/**
 * Wrapper for hashing pairs
 */
struct hash_pair {
  inline size_t operator () (const std::pair<uintE,uintE> & a) {return ((pbbs::hash64(a.first) * 3) ^ (pbbs::hash64(a.second) >> 32));}
};

/**
 * TriangleCounting stores the graph and HSet, and maintains the number of Triangles.
 */
struct TriangleCounting {
  HSet* hset;
  sparse_table<std::pair<uintE,uintE>,uintE,hash_pair> wedges;
  uintE stored;

  uintE total;

  TriangleCounting(HSet* _hset) {
    hset = _hset;
    // Actual number of wedges pair (u,v)
    total = 0;
    // Number of stored wedges pair (u,v). Some might be 0
    stored = 0;

    wedges = make_sparse_table<std::pair<uintE,uintE>,uintE,hash_pair>(hset -> G -> n,std::make_tuple(std::make_pair(UINT_E_MAX,UINT_E_MAX),0),hash_pair());

    initialize(hset -> G -> n);
  }

  void del() {
    wedges.del();
  }

  /**
  * Converts a sequence to dynamic array.
  * 
  * @param s Input Sequence
  * @return dynamic version of the sequence
  */

  pbbslib::dyn_arr<std::pair<uintE, uintE>> seq2da(sequence<std::pair<uintE,uintE>> s) {
      pbbslib::dyn_arr<std::pair<uintE,uintE>> res = pbbslib::dyn_arr<std::pair<uintE,uintE>>(s.size());

      par_for(0,s.size(),[&](size_t j) {res.A[j] = s[j];});
      res.size = s.size();

      return res;
  }

  /**
  * Concatenates a dynamic array of dynamic arrays
  * 
  * @param s A dynamic array of dynamic arrays
  * @return Concatenated result of the dynamic arrays
  */
  pbbslib::dyn_arr<std::pair<uintE,uintE>> concatSeq(pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>> & s) {
    sequence<uintE> sizes = sequence<uintE>(s.size + 1,[&](size_t i) {
      if(i == s.size) {
        return (size_t)0;
      }
      return s.A[i].capacity;});
    pbbslib::scan_add_inplace(sizes);

    pbbslib::dyn_arr<std::pair<uintE,uintE>> res = pbbslib::dyn_arr<std::pair<uintE,uintE>>(sizes[sizes.size() - 1]);
    par_for(0,s.size,[&](size_t i) {
      par_for(0,(s.A[i]).size,[&](size_t j) {
        res.A[sizes[i] + j] = s.A[i].A[j];
      }); 
    });
    res.size = res.capacity;

    s.clear();
    return res;
  }

  /**
  * Add the second dynamic array to the tail of the first one.
  * 
  * @param f The first dynamic array
  * @param s The second dyamic array
  * @return combined version of the two dynamic arrays.
  */
  pbbslib::dyn_arr<std::pair<uintE,uintE>> concat2Seq(pbbslib::dyn_arr<std::pair<uintE,uintE>> & f,pbbslib::dyn_arr<std::pair<uintE,uintE>> & s) {

    pbbslib::dyn_arr<std::pair<uintE,uintE>> res = pbbslib::dyn_arr<std::pair<uintE,uintE>>(f.size + s.size);
    res.size = f.size + s.size;

    par_for(0,res.size,[&](size_t i) {
      if(i < f.size) {
        res.A[i] = f.A[i];
      }else{
        res.A[i] = s.A[i - f.size];
      }
    });
    return res;
  }

  /**
  * Resizes the Wedges to an apporpriate capacity
  * 
  * @param amount The new expected number of wedges pairs
  */
  void resizeWedges(uintE amount) {
    amount = std::max(amount,(uintE) 4);
    if(amount << 1 <= wedges.m && amount << 2 >= wedges.m) return;
    
    auto entries = wedges.entries();
    wedges.del();
    pbbs::sequence<uintE> added = pbbs::sequence<uintE>(entries.size());

    wedges = make_sparse_table<std::pair<uintE,uintE>,uintE,hash_pair>(2 * amount,std::make_tuple(std::make_pair(UINT_E_MAX,UINT_E_MAX),0),hash_pair());
    par_for(0, entries.size(),[&](size_t i) { 
      if(std::get<1>(entries[i])) {
        wedges.insert(entries[i]);
        added[i] = 1;
      }
      else added[i] = 0;
    });
    stored = pbbslib::reduce_add(added);
    added.clear();
    entries.clear();
  }

  /**
  * Prints the elements of HSet.
  */
  void printH() {
    for(int i = 0;i < hset -> getH().size();++i) {
      std::cout << hset -> getH()[i] << " ";
    }
    std::cout << endl;
  }

  /**
  * Output the count for each pair of wedges pair (i,j)
  */
  void printAllWedges() {
    for(int i = 0; i < hset -> G -> n;++i) {
      for(int j = i + 1;j < hset -> G -> n;++j) {
        cout << i << " - " << j << " : " << wedges.find(std::make_pair(i,j),0) << endl;
      }
    }
    cout << "-----------------------------" << endl;
  }

  /**
  * Adds c to wedges pair (a,b)
  * 
  * @param a the first index of the pair
  * @param b the second index of the pair
  * @param c the amount added, has to be non-negative
  */
  void addToWedge(uintE a,uintE b,uintE c) {
    uintE res = wedges.find(std::make_pair(a,b),0);
    if(res == 0) {
      wedges.insert(std::make_tuple(std::make_pair(a,b),c));
      wedges.change(std::make_pair(a,b),c);
    }else{
      wedges.change(std::make_pair(a,b),res + c);
    }
  }

  /**
  * Removes c to wedges pair (a,b)
  * 
  * @param a the first index of the pair
  * @param b the second index of the pair
  * @param c the amount removed, has to be non-negative
  */
  void remFrWedge(uintE a,uintE b,uintE c) {
    uintE res = wedges.find(std::make_pair(a,b),0);
    if(res == 0) {
      wedges.insert(std::make_tuple(std::make_pair(a,b),res - c));
    }else{
      wedges.change(std::make_pair(a,b),res - c);
    }
  }

  /**
  * Given a sequence of wedges, add them to the wedges pair
  * 
  * @param seq the sequence of pairs of the wedges that need to be removed
  */
  void addWedges(sequence<std::pair<uintE,uintE>> & seq) {
    resizeWedges(stored + seq.size());
    stored += seq.size();

    par_for(0,seq.size(),[&](size_t i) {
      if(i != 0) {
        if(seq[i] != seq[i - 1]) {
          // wedges[seq[i - 1].first][seq[i - 1].second] += i;
          addToWedge(seq[i - 1].first,seq[i - 1].second,i);
        }
      }
    });

    if(seq.size() != 0) addToWedge(seq[seq.size() - 1].first,seq[seq.size() - 1].second,seq.size());

    par_for(0,seq.size(),[&](size_t i) {
      if(i != 0) {
        if(seq[i] != seq[i - 1]) {
          // wedges[seq[i].first][seq[i].second] -= i;
          remFrWedge(seq[i].first,seq[i].second,i);
        }
      }
    });

  }

  /**
  * Given a sequence of wedges, remove them from the wedges pair
  * 
  * @param seq the sequence of pairs of the wedges that need to be removed
  */
  void removeWedges(sequence<std::pair<uintE,uintE>> & seq) {
    par_for(0,seq.size(),[&](size_t i) {
      if(i != 0) {
        if(seq[i] != seq[i - 1]) {
          addToWedge(seq[i].first,seq[i].second,i);
        }
      }
    });

    par_for(0,seq.size(),[&](size_t i) {
      if(i != 0) {
        if(seq[i] != seq[i - 1]) {
          remFrWedge(seq[i - 1].first,seq[i - 1].second,i);
        }
      }
    });

    if(seq.size() != 0) remFrWedge(seq[seq.size() - 1].first,seq[seq.size() - 1].second,seq.size());
  }


  /**
  * After modifying the H-Sets, this function adjusts the wedges so that everything is according to the new H Set.
  * 
  * @param originalH the original HSet elements in sequence
  * @param hs the new HSet elements in sequence
  * @param allEdges all the modified edges stored in a sparse table for O(1) lookup
  * @param originalHSet the original HSet elements stored in a sparse tbale for O(1) lookup
  */
  void adjustHSetWedges(sequence<uintE>& originalH, sequence<uintE>& hs,sparse_table<std::pair<uintE,uintE>, bool, hash_pair>& allEdges, sparse_table<uintE,bool,hash_uintE>& originalHSet) {

    // Wedges that need to be added due to being removed form HSet
    pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>> newWedges = pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>>(originalH.size());
    newWedges.size = originalH.size();
    par_for(0,originalH.size(),[&](size_t i) {
      // If it was in H-Set but not anymore.
      if(!hset -> contains(originalH[i])) {
        // Need to add wedges back
        // I need to find the nodes not in special set first
        auto es = hset -> G -> v_data.A[originalH[i]].neighbors.entries();
        sequence<uintE> specialSet = filter(sequence<uintE>(es.size(),[&](size_t j) {
          if(!std::get<1>(es[j]) || allEdges.find(std::make_pair(originalH[i],std::get<0>(es[j])),false)) {
            return UINT_E_MAX;
          }
          return std::get<0>(es[j]);
        }),[&](uintE u) {return u != UINT_E_MAX;});
        newWedges.A[i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(specialSet.size() * (specialSet.size() - 1));
        par_for(0,specialSet.size(),[&](size_t j) {
          newWedges.A[i].size = specialSet.size() * (specialSet.size() - 1);
          par_for(0,specialSet.size(),[&](size_t k) {
            if(k < j) {
              // cout << specialSet[j] << " " << specialSet[k] << endl;
              newWedges.A[i].A[j * (specialSet.size() - 1) + k] = std::make_pair(specialSet[j],specialSet[k]);
            }else if(k > j) {
              // cout << specialSet[j] << " " << specialSet[k] << endl;
              newWedges.A[i].A[j * (specialSet.size() - 1) + k - 1] = std::make_pair(specialSet[j],specialSet[k]);
            }
          });
        });

        es.clear();
        specialSet.clear();
      }else{
        newWedges.A[i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);
      }
    });


    // Wedges that need to be removed due to being added to HSet
    pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>> removedWedges = pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>>(hs.size());
    removedWedges.size = hs.size();
    par_for(0,hs.size(),[&](size_t i) {
      // If it wasn't in H-Set but now is.
      if(!(originalHSet.find(hs[i],false))) {
        // Need to remove the wedges
        // I need to find the nodes not in special set first
        auto es = hset -> G -> v_data.A[hs[i]].neighbors.entries();
        sequence<uintE> specialSet = pbbs::filter(sequence<uintE>(es.size(),[&](size_t j) {
          if(!std::get<1>(es[j]) || allEdges.find(std::make_pair(hs[i],std::get<0>(es[j])),false)) {
            return UINT_E_MAX;
          }
          return std::get<0>(es[j]);
        }),[&](uintE u) {return u != UINT_E_MAX;});


        removedWedges.A[i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(specialSet.size() * (specialSet.size() - 1));
        par_for(0,specialSet.size(),[&](size_t j) {
          removedWedges.A[i].size = specialSet.size() * (specialSet.size() - 1);
          par_for(0,specialSet.size(),[&](size_t k) {
            if(k < j) {
              removedWedges.A[i].A[j * (specialSet.size() - 1) + k] = std::make_pair(specialSet[j],specialSet[k]);
            }else if(k > j) {
              removedWedges.A[i].A[j * (specialSet.size() - 1) + k - 1] = std::make_pair(specialSet[j],specialSet[k]);
            }
          });
        });
      }else{
        removedWedges.A[i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);
      }
    });

    sequence<std::pair<uintE,uintE>> allNewWedges = merge_sort(concatSeq(newWedges).to_seq(),[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a > b;});
    addWedges(allNewWedges);

    sequence<std::pair<uintE,uintE>> allRemovedWedges = merge_sort(concatSeq(removedWedges).to_seq(),[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a > b;});
    removeWedges(allRemovedWedges);

    newWedges.del();
    allNewWedges.clear();
    removedWedges.del();
    allRemovedWedges.clear();
  }

  /**
  * Initializes as a graph of n vertices(or with id of at most n - 1)
  *
  * @param n the expected max id of the vertices in the grpah.
  */
  void initialize(size_t n) {
    resizeWedges(n * hset -> hindex);
    sequence<uintE> nodes = sequence<uintE>(n,[&](size_t i) {return i;});
    hset -> insertVertices(nodes);
    nodes.clear();
  }

  /**
  * Adds edges to the graph, which will also modify the H-Set and update the triangle counts
  *
  * @param edges The sequence of edges that need to be added
  */
  void addEdges(sequence<std::pair<uintE,uintE>> edges) {
    edges = pbbs::filter(edges, [&] (std::pair<uintE, uintE> e) { return !hset -> G -> existEdge(e.first, e.second); } );

    // STEP 1 - Adjust the H Set

    // AllEdges in the form of a sparse table, allowing O(1) lookup if an edge exists.
    auto allEdges = make_sparse_table<std::pair<uintE,uintE>,bool,hash_pair>(2 * edges.size() + 1,std::make_tuple(std::make_pair(UINT_E_MAX,UINT_E_MAX),false),hash_pair());
    par_for(0,edges.size(),[&](size_t i) {
      allEdges.insert(std::make_tuple(edges[i],true));
      allEdges.insert(std::make_tuple(std::make_pair(edges[i].second,edges[i].first),true));
    });

    // The original HSet elements
    sequence<uintE> originalH = hset -> getH();

    hset -> insertEdges(edges);

    // The new HSet elements
    sequence<uintE> hs = hset -> getH();

    // The original HSet elements in a sparse table, allowing O(1) lookup if an element was originally in the hset.
    auto originalHSet = make_sparse_table<uintE,bool,hash_uintE>(originalH.size() * 2 + 1,std::make_tuple(UINT_E_MAX,false),hash_uintE());

    par_for(0,originalH.size(),[&](size_t i) {
      originalHSet.insert(std::make_tuple(originalH[i],true));
    });  

    // Adjusts the wedges now that we have updated the HSet
    adjustHSetWedges(originalH,hs,allEdges,originalHSet);

    // triangles added by wedges
    pbbslib::dyn_arr<uintE> wedgeTriangles = pbbslib::dyn_arr<uintE>(edges.size());
    // triangles added not by h-set but by multiple added edges
    pbbslib::dyn_arr<uintE> multiTriangles = pbbslib::dyn_arr<uintE>(edges.size());

    par_for(0,edges.size(),[&](size_t i) {
      wedgeTriangles.A[i] = 0;
      multiTriangles.A[i] = 0;
    });

    // newly formed wedges

    // STEP 2 - Find all the triangles/wedges by either wedges or entities formed solely by the newly added edges
    pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>> newWedges = pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>>(2 * edges.size());

    par_for(0,edges.size(),[&](size_t i) {
      uintE u = std::get<1>(edges[i]);
      uintE v = std::get<0>(edges[i]);

      newWedges.A[2 * i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);
      newWedges.A[2 * i + 1] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);

      // Case 1, added by wedges
      wedgeTriangles.A[i] = wedges.find(std::make_pair(v,u),0);

      if(!hset -> contains(v)) {
        auto edgesV = hset -> G -> v_data.A[v].neighbors.entries();

        pbbslib::dyn_arr<uintE> multiTrianglesV = pbbslib::dyn_arr<uintE>(edgesV.size());

        pbbslib::dyn_arr<std::pair<uintE,uintE>> newWedgesV = pbbslib::dyn_arr<std::pair<uintE,uintE>>(2 * edgesV.size());
        par_for(0,edgesV.size(),[&](size_t i) {
          multiTrianglesV.A[i] = 0;
          uintE next = std::get<0>(edgesV[i]);
          if(next == u || !std::get<1>(edgesV[i])) {
            newWedgesV.A[2 * i] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
            newWedgesV.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
            return;
          }
          if(allEdges.find(std::make_pair(v,next),false)) {
            // Both are added edges

            newWedgesV.A[2 * i] = std::make_pair(u,next);
            newWedgesV.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
          }else{
            newWedgesV.A[2 * i] = std::make_pair(u,next);
            newWedgesV.A[2 * i + 1] = std::make_pair(next,u);
          }

          if(hset -> G -> existEdge(v,next) && hset -> G -> existEdge(next,u)) {
            // A triangle!
            if(allEdges.find(std::make_pair(v,next),false) && !(allEdges.find(std::make_pair(u,next),false))) {
              // (u,v) and (v,next) are added edges but not (next,u).
              // Possible cases are 2 and 3
              if(!hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)) {
                // Case 2, 2 added edges and no h-set
                multiTrianglesV.A[i] = 3;
              }else if(hset -> contains(v) && !hset -> contains(u) && !hset -> contains(next)){
                // Case 3, 2 added edges and joint node is in h-set
                multiTrianglesV.A[i] = 6;
              }
            }else if(allEdges.find(std::make_pair(u,next),false) && !(allEdges.find(std::make_pair(v,next),false))) {
              // (u,v) and (u,next) are added edges but not (v,next).
              // same with above
              if(!hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)) {
                // Case 2, 2 added edges and no h-set
                multiTrianglesV.A[i] = 3;
              }else if(hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)){
                // Case 3, 2 added edges and joint node is in h-set
                multiTrianglesV.A[i] = 6;
              }
            }else if(allEdges.find(std::make_pair(u,next),false) && allEdges.find(std::make_pair(v,next),false)) {
              if(!hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)) {
                // Case 8, 3 added edges and no nodes in h-set
                multiTrianglesV.A[i] = 2;
              }
            }
          }
        });
        newWedgesV.size = newWedgesV.capacity;

        newWedges.A[2 * i] = seq2da(filter(newWedgesV.to_seq(),[&](std::pair<uintE,uintE> cur) {
          if((cur.first) == UINT_E_MAX && (cur.second) == UINT_E_MAX) return false;
          return true;
        }));
        multiTrianglesV.size = multiTrianglesV.capacity;
        multiTriangles.A[i] = pbbs::reduce(multiTrianglesV.to_seq(),pbbs::addm<uintE>());

        multiTrianglesV.del();
        newWedgesV.del();
      }

      if(!hset -> contains(u)) {
        auto edgesU = hset -> G -> v_data.A[u].neighbors.entries();


        pbbslib::dyn_arr<uintE> multiTrianglesU = pbbslib::dyn_arr<uintE>(edgesU.size());

        // newly formed wedges (one existing, one added)
        pbbslib::dyn_arr<std::pair<uintE,uintE>> newWedgesU = pbbslib::dyn_arr<std::pair<uintE,uintE>>(2 * edgesU.size());

        par_for(0,edgesU.size(),[&](size_t i) {
          multiTrianglesU.A[i] = 0;
          uintE next = std::get<0>(edgesU[i]);
          if(next == v || !std::get<1>(edgesU[i])) {

            newWedgesU.A[2 * i] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
            newWedgesU.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
            return;
          }
          if(allEdges.find(std::make_pair(u,next),false)) {
            // Both are added edges
            newWedgesU.A[2 * i] = std::make_pair(v,next);
            newWedgesU.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
          }else{
            newWedgesU.A[2 * i] = std::make_pair(v,next);
            newWedgesU.A[2 * i + 1] = std::make_pair(next,v);
          }

          if(hset -> G -> existEdge(next,u) && hset -> G -> existEdge(next,v)) {
            // A triangle!
            if(allEdges.find(std::make_pair(v,next),false) && !(allEdges.find(std::make_pair(u,next),false))) {
              // (u,v) and (v,next) are added edges but not (next,u).
              // Possible cases are 2 and 3
              if(!hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)) {
                // Case 2, 2 added edges and no h-set
                multiTrianglesU.A[i] = 3;
              }else if(hset -> contains(v) && !hset -> contains(u) && !hset -> contains(next)){
                // Case 3, 2 added edges and joint node is in h-set
                multiTrianglesU.A[i] = 6;
              }
            }else if(allEdges.find(std::make_pair(u,next),false) && !(allEdges.find(std::make_pair(v,next),false))) {
              // (u,v) and (u,next) are added edges but not (v,next)/
              // same with above

              if(!hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)) {
                // Case 2, 2 added edges and no h-set
                multiTrianglesU.A[i] = 3;
              }else if(hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)){
                // Case 3, 2 added edges and joint node is in h-set
                multiTrianglesU.A[i] = 6;
              }
            }else if(allEdges.find(std::make_pair(u,next),false) && allEdges.find(std::make_pair(v,next),false)) {
              if(!hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)) {
                // Case 8, 3 added edges and no nodes in h-set
                multiTrianglesU.A[i] = 2;
              }
            }
          }
        });

        newWedgesU.size = newWedgesU.capacity;

        newWedges.A[2 * i + 1] = seq2da(filter(newWedgesU.to_seq(),[&](std::pair<uintE,uintE> cur) {
          if((cur.first) == UINT_E_MAX && (cur.second) == UINT_E_MAX) return false;
          return true;
        }));
        multiTrianglesU.size = multiTrianglesU.capacity;
        multiTriangles.A[i] += pbbs::reduce(multiTrianglesU.to_seq(),pbbs::addm<uintE>());

        multiTrianglesU.del();
        newWedgesU.del();
      }
    });


    // STEP 3 - Find triangles formed with the H Set
    pbbslib::dyn_arr<uintE> hsetTriangles = pbbslib::dyn_arr<uintE>(edges.size());
    par_for(0,edges.size(),[&](size_t i) {
      hsetTriangles.A[i] = 0;
      pbbslib::dyn_arr<uintE> uvhsetTriangles = pbbslib::dyn_arr<uintE>(hs.size());
      int u = edges[i].first;
      int v = edges[i].second;
      par_for(0,hs.size(),[&](size_t j) {
        uvhsetTriangles.A[j] = 0;
        int next = hs[j];
        if(u == next || v == next) return;
        if(hset -> G -> existEdge(u,next) && hset -> G -> existEdge(v,next)) {
          // Triangle!
          // I decided to this instead, since if i do case work, there would be like 26 if statements
          int counter = 1;
          if(allEdges.find(std::make_pair(next,u),false) && hset -> contains(v)) {
            counter++;
          }
          if(allEdges.find(std::make_pair(v,next),false) && hset -> contains(u)) {
            counter++;
          }
          uvhsetTriangles.A[j] = 6 / counter;
        }
      });
      uvhsetTriangles.size = uvhsetTriangles.capacity;
      hsetTriangles.A[i] = pbbs::reduce(uvhsetTriangles.to_seq(),pbbs::addm<uintE>());
    });

    // STEP 4, update the counts and the wedges

    newWedges.size = newWedges.capacity;
    wedgeTriangles.size = wedgeTriangles.capacity;
    total += pbbs::reduce(wedgeTriangles.to_seq(),pbbs::addm<uintE>());
    wedgeTriangles.del();

    multiTriangles.size = multiTriangles.capacity;
    total += pbbs::reduce(multiTriangles.to_seq(),pbbs::addm<uintE>()) / 12;
    multiTriangles.del();

    hsetTriangles.size = hsetTriangles.capacity;
    total += pbbs::reduce(hsetTriangles.to_seq(),pbbs::addm<uintE>()) / 6;
    hsetTriangles.del();

    // STEP 4

    sequence<std::pair<uintE,uintE>> allNewWedges = merge_sort(concatSeq(newWedges).to_seq(),[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a > b;});
    addWedges(allNewWedges);
    newWedges.del();
    allNewWedges.clear();
    allEdges.del();
    originalH.clear();
    hs.clear();
    originalHSet.del();
  }


  /**
  * Removes edges from the graph, which will also modify the H-Set and update the triangle counts
  *
  * @param edges The sequence of edges that need to be added
  */
  void removeEdges(sequence<std::pair<uintE,uintE>> edges) {

    edges = pbbs::filter(edges, [&] (std::pair<uintE, uintE> e) { return hset -> G -> existEdge(e.first, e.second); } );

    auto allEdges = make_sparse_table<std::pair<uintE,uintE>,bool,hash_pair>
        (2 * edges.size() + 1,std::make_tuple(std::make_pair(UINT_E_MAX,UINT_E_MAX),false),hash_pair());
    par_for(0,edges.size(),[&](size_t i) {
      allEdges.insert(std::make_tuple(edges[i],true));
      allEdges.insert(std::make_tuple(std::make_pair(edges[i].second,edges[i].first),true));
    });

    sequence<uintE> originalH = hset -> getH();

    // triangles removed by wedges
    pbbslib::dyn_arr<uintE> wedgeTriangles = pbbslib::dyn_arr<uintE>(edges.size());
    // triangles removed not by h-set
    pbbslib::dyn_arr<uintE> multiTriangles = pbbslib::dyn_arr<uintE>(edges.size());

    par_for(0,edges.size(),[&](size_t i) {
      wedgeTriangles.A[i] = 0;
      multiTriangles.A[i] = 0;
    });

    // wedges that need to be removed
    // STEP 1 - Finds all the triangles/wedges erased due to either the wedges or by the removed edges themselves
    pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>> removedWedges = pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>>(2 * edges.size());

    par_for(0,edges.size(),[&](size_t i) {
      uintE u = std::get<1>(edges[i]);
      uintE v = std::get<0>(edges[i]);

      removedWedges.A[2 * i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);
      removedWedges.A[2 * i + 1] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);

      // Case 1, removed by wedges
      wedgeTriangles.A[i] = wedges.find(std::make_pair(v,u), 0);

      if(!hset -> contains(v)) {
        auto edgesV = hset -> G -> v_data.A[v].neighbors.entries();

        pbbslib::dyn_arr<uintE> multiTrianglesV = pbbslib::dyn_arr<uintE>(edgesV.size());

        pbbslib::dyn_arr<std::pair<uintE,uintE>> removedWedgesV = pbbslib::dyn_arr<std::pair<uintE,uintE>>(2 * edgesV.size());
        par_for(0,edgesV.size(),[&](size_t i) {
          multiTrianglesV.A[i] = 0;
          uintE next = std::get<0>(edgesV[i]);
          if(next == u || !std::get<1>(edgesV[i])) {
            removedWedgesV.A[2 * i] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
            removedWedgesV.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
            return;
          }
          if(allEdges.find(std::make_pair(v,next),false)) {
            // Both are removed edges
            removedWedgesV.A[2 * i] = std::make_pair(u,next);
            removedWedgesV.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
          }else{
            removedWedgesV.A[2 * i] = std::make_pair(u,next);
            removedWedgesV.A[2 * i + 1] = std::make_pair(next,u);
          }

          if(hset -> G -> existEdge(v,next) && hset -> G -> existEdge(next,u)) {
            // A triangle!
            if(allEdges.find(std::make_pair(v,next),false) && !(allEdges.find(std::make_pair(u,next),false))) {
              // (u,v) and (v,next) are added edges but not (next,u).
              // Possible cases are 2 and 3
              if(!hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)) {
                // Case 2, 2 added edges and no h-set
                multiTrianglesV.A[i] = 3;
              }else if(hset -> contains(v) && !hset -> contains(u) && !hset -> contains(next)){
                // Case 3, 2 added edges and joint node is in h-set
                multiTrianglesV.A[i] = 6;
              }
            }else if(allEdges.find(std::make_pair(u,next),false) && !(allEdges.find(std::make_pair(v,next),false))) {
              // (u,v) and (u,next) are added edges but not (v,next).
              // same with above
              if(!hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)) {
                // Case 2, 2 added edges and no h-set
                multiTrianglesV.A[i] = 3;
              }else if(hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)){
                // Case 3, 2 added edges and joint node is in h-set
                multiTrianglesV.A[i] = 6;
              }
            }else if(allEdges.find(std::make_pair(u,next),false) && allEdges.find(std::make_pair(v,next),false)) {
              if(!hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)) {
                // Case 8, 3 added edges and no nodes in h-set
                multiTrianglesV.A[i] = 4;
              }
            }
          }
        });
        removedWedgesV.size = removedWedgesV.capacity;

        removedWedges.A[2 * i] = seq2da(filter(removedWedgesV.to_seq(),[&](std::pair<uintE,uintE> cur) {
          if((cur.first) == UINT_E_MAX && (cur.second) == UINT_E_MAX) return false;
          return true;
        }));
        multiTrianglesV.size = multiTrianglesV.capacity;
        multiTriangles.A[i] = pbbs::reduce(multiTrianglesV.to_seq(),pbbs::addm<uintE>());

        multiTrianglesV.del();
        removedWedgesV.del();
      }

      if(!hset -> contains(u)) {
        auto edgesU = hset -> G -> v_data.A[u].neighbors.entries();


        pbbslib::dyn_arr<uintE> multiTrianglesU = pbbslib::dyn_arr<uintE>(edgesU.size());

        // newly formed wedges (one existing, one added)
        pbbslib::dyn_arr<std::pair<uintE,uintE>> removedWedgesU = pbbslib::dyn_arr<std::pair<uintE,uintE>>(2 * edgesU.size());

        par_for(0,edgesU.size(),[&](size_t i) {
          multiTrianglesU.A[i] = 0;
          uintE next = std::get<0>(edgesU[i]);
          if(next == v || !std::get<1>(edgesU[i])) {
            removedWedgesU.A[2 * i] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
            removedWedgesU.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
            return;
          }
          if(allEdges.find(std::make_pair(u,next),false)) {
            // Both are added edges
            removedWedgesU.A[2 * i] = std::make_pair(v,next);
            removedWedgesU.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
          }else{
            removedWedgesU.A[2 * i] = std::make_pair(v,next);
            removedWedgesU.A[2 * i + 1] = std::make_pair(next,v);
          }

          if(hset -> G -> existEdge(v,next) && hset -> G -> existEdge(next,u)) {
            // A triangle!
            if(allEdges.find(std::make_pair(v,next),false) && !(allEdges.find(std::make_pair(u,next),false))) {
              // (u,v) and (v,next) are added edges but not (next,u).
              // Possible cases are 2 and 3
              if(!hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)) {
                // Case 2, 2 added edges and no h-set
                multiTrianglesU.A[i] = 3;
              }else if(hset -> contains(v) && !hset -> contains(u) && !hset -> contains(next)){
                // Case 3, 2 added edges and joint node is in h-set
                multiTrianglesU.A[i] = 6;
              }
            }else if(allEdges.find(std::make_pair(u,next),false) && !(allEdges.find(std::make_pair(v,next),false))) {
              // (u,v) and (u,next) are added edges but not (v,next).
              // same with above
              if(!hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)) {
                // Case 2, 2 added edges and no h-set
                multiTrianglesU.A[i] = 3;
              }else if(hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)){
                // Case 3, 2 added edges and joint node is in h-set
                multiTrianglesU.A[i] = 6;
              }
            }else if(allEdges.find(std::make_pair(u,next),false) && allEdges.find(std::make_pair(v,next),false)) {
              if(!hset -> contains(u) && !hset -> contains(v) && !hset -> contains(next)) {
                // Case 8, 3 added edges and no nodes in h-set
                multiTrianglesU.A[i] = 4;
              }
            }
          }
        });

        removedWedgesU.size = removedWedgesU.capacity;

        removedWedges.A[2 * i + 1] = seq2da(filter(removedWedgesU.to_seq(),[&](std::pair<uintE,uintE> cur) {
          if((cur.first) == UINT_E_MAX && (cur.second) == UINT_E_MAX) return false;
          return true;
        }));
        multiTrianglesU.size = multiTrianglesU.capacity;
        multiTriangles.A[i] += pbbs::reduce(multiTrianglesU.to_seq(),pbbs::addm<uintE>());

        multiTrianglesU.del();
        removedWedgesU.del();
      }
    });


    // STEP 2 - triangles removed with h-set
    pbbslib::dyn_arr<long long> hsetTriangles = pbbslib::dyn_arr<long long>(edges.size());
    par_for(0,edges.size(),[&](size_t i) {
      hsetTriangles.A[i] = 0;
      pbbslib::dyn_arr<long long> uvhsetTriangles = pbbslib::dyn_arr<long long>(originalH.size());
      int u = edges[i].first;
      int v = edges[i].second;
      par_for(0,originalH.size(),[&](size_t j) {
        uvhsetTriangles.A[j] = 0;
        int next = originalH[j];
        if(hset -> G -> existEdge(u,next) && hset -> G -> existEdge(v,next)) {
          // Triangle!
          // I decided to this instead, since if i do case work, there would be like 26 if statements
          // TODO: What is the meaning of these counters?
          int counter = 1;
          int counter2 = -1;
          if(allEdges.find(std::make_pair(u,next),false) && hset -> contains(v)) {
            counter++;
          }
          if(allEdges.find(std::make_pair(v,next),false) && hset -> contains(u)) {
            counter++;
          }
          if(allEdges.find(std::make_pair(u,next),false) && !hset -> contains(v)) {
            counter2++;
          }
          if(allEdges.find(std::make_pair(v,next),false) && !hset -> contains(u)) {
            counter2++;
          }
          uvhsetTriangles.A[j] = 6 * counter2 / counter;
        }
      });
      uvhsetTriangles.size = uvhsetTriangles.capacity;
      hsetTriangles.A[i] = pbbs::reduce(uvhsetTriangles.to_seq(),pbbs::addm<long long>());
    });

    // STEP 3 - Updates all the counts

    removedWedges.size = removedWedges.capacity;
    wedgeTriangles.size = wedgeTriangles.capacity;
    total -= pbbs::reduce(wedgeTriangles.to_seq(),pbbs::addm<uintE>());
    wedgeTriangles.del();
    multiTriangles.size = multiTriangles.capacity;
    total += pbbs::reduce(multiTriangles.to_seq(),pbbs::addm<uintE>()) / 12;
    multiTriangles.del();
    hsetTriangles.size = hsetTriangles.capacity;
    total += pbbs::reduce(hsetTriangles.to_seq(),pbbs::addm<long long>()) / 6;
    hsetTriangles.del();

    // STEP 4 - Updates the HSet and the wedges

    hset -> eraseEdges(edges);

    sequence<uintE> hs = hset -> getH();

    auto originalHSet = make_sparse_table<uintE,bool,hash_uintE>(originalH.size() * 2 + 1,std::make_tuple(UINT_E_MAX,false),hash_uintE());

    par_for(0,originalH.size(),[&](size_t i) {
      originalHSet.insert(std::make_tuple(originalH[i],true));
    });

    sequence<std::pair<uintE,uintE>> allRemovedWedges = merge_sort(concatSeq(removedWedges).to_seq(),[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a > b;});
    removeWedges(allRemovedWedges);
    adjustHSetWedges(originalH,hs,allEdges,originalHSet); 
    removedWedges.del();
    allRemovedWedges.clear();
    allEdges.del();
    originalH.clear();
    hs.clear();
    originalHSet.del();
  }
};
