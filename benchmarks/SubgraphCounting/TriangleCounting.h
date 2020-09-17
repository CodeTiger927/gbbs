
struct hash_pair {
  inline size_t operator () (const std::pair<uintE,uintE> & a) {return ((pbbs::hash64(a.first) * 3) ^ (pbbs::hash64(a.second) >> 32));}
};

// TODO: You shouldn't have global variables like this. You can make it a function
// if you'd like.
// Converts sequences to dynamic array
auto seq2da = [&](sequence<std::pair<uintE,uintE>> s) -> pbbslib::dyn_arr<std::pair<uintE,uintE>> {pbbslib::dyn_arr<std::pair<uintE,uintE>> res = pbbslib::dyn_arr<std::pair<uintE,uintE>>(s.size()); par_for(0,s.size(),[&](size_t j) {res.A[j] = s[j];}); res.size = s.size(); return res;};

struct TriangleCounting {
// Here is a lot of the things that need to be manually fixed
uintE N = 1005;

dynamic_symmetric_graph<dynamic_symmetric_vertex,pbbs::empty> dsg = createEmptyDynamicSymmetricGraph<dynamic_symmetric_vertex,pbbs::empty>();
HSetAlex hset = HSetAlex(&dsg);
// uintE wedges[1005][1005];
sparse_table<std::pair<uintE,uintE>,uintE,hash_pair> wedges = make_sparse_table<std::pair<uintE,uintE>,uintE,hash_pair>(N,std::make_tuple(std::make_pair(UINT_E_MAX,UINT_E_MAX),0),hash_pair());
uintE total = 0;

TriangleCounting() {

}

// Concatenates a lot of dynamic dynamic arrays into one
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
  return res;
}

// Concatenates two dynamic arrays into one
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


// Function for debugging. Prints all the elements in H Set
void printH() {
  for(int i = 0;i < hset.allH().size();++i) {
    std::cout << hset.allH()[i] << " ";
  }
  std::cout << endl;
}

// Function for debugging. Prints the counts for all the wedges
void printAllWedges() {
  for(int i = 0; i < N;++i) {
    for(int j = i + 1;j < N;++j) {
      cout << i << " - " << j << " : " << wedges.find(std::make_pair(i,j),0) << endl;
    }
  }
  cout << "-----------------------------" << endl;
}

// Resizes Wedges sparse table, and ensures that the sparsity is always at least 2, so that find functions takes less than O(2) expected value
void resizeWedges(uintE amount) {
  amount = std::max(amount,(uintE)4);
  if(amount << 1 <= wedges.m && amount << 2 >= wedges.m) return;
  auto entries = wedges.entries();
  wedges = make_sparse_table<std::pair<uintE,uintE>,uintE,hash_pair>(2 * amount,std::make_tuple(std::make_pair(UINT_E_MAX,UINT_E_MAX),0),hash_pair());
  par_for(0,entries.size(),[&](size_t i) {if(std::get<1>(entries[i])) wedges.insert(entries[i]);});
}

// Adds c to wedges (a,b)
void addToWedge(uintE a,uintE b,uintE c) {
  uintE res = wedges.find(std::make_pair(a,b),0);
  if(res == 0) {
    wedges.insert(std::make_tuple(std::make_pair(a,b),c));
    wedges.change(std::make_pair(a,b),c);
  }else{
    wedges.change(std::make_pair(a,b),res + c);
  }
}

// Subtracts c from wedges (a,b)
void remFrWedge(uintE a,uintE b,uintE c) {
  uintE res = wedges.find(std::make_pair(a,b),0);
  if(res == 0) {
    wedges.insert(std::make_tuple(std::make_pair(a,b),res - c));
  }else{
    wedges.change(std::make_pair(a,b),res - c);
  }
}

// Given a series of wedges, add them
void addWedges(sequence<std::pair<uintE,uintE>> seq) {
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

// Given a series of wedges, remove them
void removeWedges(sequence<std::pair<uintE,uintE>> seq) {
  par_for(0,seq.size(),[&](size_t i) {
    if(i != 0) {
      if(seq[i] != seq[i - 1]) {
        addToWedge(seq[i].first,seq[i].second,i);
        // wedges[seq[i].first][seq[i].second] += i;
      }
    }
  });

  par_for(0,seq.size(),[&](size_t i) {
    if(i != 0) {
      if(seq[i] != seq[i - 1]) {
        remFrWedge(seq[i - 1].first,seq[i - 1].second,i);
        // wedges[seq[i - 1].first][seq[i - 1].second] -= i;
      }
    }
  });

  if(seq.size() != 0) remFrWedge(seq[seq.size() - 1].first,seq[seq.size() - 1].second,seq.size());
}


// After modifying the H-Sets, this function adjusts the wedges so that everything is according to the new H Set.
void adjustHSetWedges(sequence<uintE>& originalH,sequence<uintE>& hs,sparse_table<std::pair<uintE,uintE>,bool,hash_pair>& allEdges,sparse_table<uintE,bool,hash_uintE>& originalHSet) {

  // Wedges that need to be added due to being removed form HSet
  pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>> dW = pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>>(originalH.size());
  dW.size = originalH.size();
  par_for(0,originalH.size(),[&](size_t i) {
    // If it was in H-Set but not anymore.
    if(!hset.inH(originalH[i])) {
      // Need to add wedges back
      // I need to find the nodes not in special set first
      auto es = dsg.v_data.A[originalH[i]].neighbors.entries();
      sequence<uintE> specialSet = filter(sequence<uintE>(es.size(),[&](size_t j) {
        if(!std::get<1>(es[j]) || allEdges.find(std::make_pair(originalH[i],std::get<0>(es[j])),false)) {
          return UINT_E_MAX;
        }
        return std::get<0>(es[j]);
      }),[&](uintE u) {return u != UINT_E_MAX;});
      dW.A[i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(specialSet.size() * (specialSet.size() - 1));
      par_for(0,specialSet.size(),[&](size_t j) {
        dW.A[i].size = specialSet.size() * (specialSet.size() - 1);
        par_for(0,specialSet.size(),[&](size_t k) {
          if(k < j) {
            // cout << specialSet[j] << " " << specialSet[k] << endl;
            dW.A[i].A[j * (specialSet.size() - 1) + k] = std::make_pair(specialSet[j],specialSet[k]);
          }else if(k > j) {
            // cout << specialSet[j] << " " << specialSet[k] << endl;
            dW.A[i].A[j * (specialSet.size() - 1) + k - 1] = std::make_pair(specialSet[j],specialSet[k]);
          }
        });
      });
    }else{
      dW.A[i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);
    }
  });


  // Wedges that need to be removed due to being added to HSet
  pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>> rW = pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>>(hs.size());
  rW.size = hs.size();
  par_for(0,hs.size(),[&](size_t i) {
    // If it wasn't in H-Set but now is.
    if(!(originalHSet.find(hs[i],false))) {
      // Need to remove the wedges
      // I need to find the nodes not in special set first
      auto es = dsg.v_data.A[hs[i]].neighbors.entries();
      sequence<uintE> specialSet = filter(sequence<uintE>(es.size(),[&](size_t j) {
        if(!std::get<1>(es[j]) || allEdges.find(std::make_pair(hs[i],std::get<0>(es[j])),false)) {
          return UINT_E_MAX;
        }
        return std::get<0>(es[j]);
      }),[&](uintE u) {return u != UINT_E_MAX;});


      rW.A[i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(specialSet.size() * (specialSet.size() - 1));
      par_for(0,specialSet.size(),[&](size_t j) {
        rW.A[i].size = specialSet.size() * (specialSet.size() - 1);
        par_for(0,specialSet.size(),[&](size_t k) {
          if(k < j) {
            rW.A[i].A[j * (specialSet.size() - 1) + k] = std::make_pair(specialSet[j],specialSet[k]);
          }else if(k > j) {
            rW.A[i].A[j * (specialSet.size() - 1) + k - 1] = std::make_pair(specialSet[j],specialSet[k]);
          }
        });
      });
    }else{
      rW.A[i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);
    }
  });

  pbbslib::dyn_arr<std::pair<uintE,uintE>> fdW = concatSeq(dW);
  sequence<std::pair<uintE,uintE>> all2 = merge_sort(fdW.to_seq(),[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a > b;});
  addWedges(all2);

  pbbslib::dyn_arr<std::pair<uintE,uintE>> finalRemove = concatSeq(rW);
  sequence<std::pair<uintE,uintE>> allR = merge_sort(finalRemove.to_seq(),[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a > b;});
  removeWedges(allR);
}

// Initializes as a graph of n vertices(or with id of at most n - 1)
void initialize(size_t n) {
  hset.resizeV(n);
  resizeWedges(n * n);
  sequence<uintE> nodes = sequence<uintE>(n,[&](size_t i) {return i;});
  sequence<std::pair<uintE,uintE>> hPairs = sequence<std::pair<uintE,uintE>>(n,[&](size_t i) {return std::make_pair(i,0);});
  hset.insertVertices(nodes);
  hset.insert(hPairs);
}

// Adds edges to the graph, which will also modify the H-Set and update the triangle counts
void addEdges(sequence<std::pair<uintE,uintE>> edges) {
  // STEP 1 - Adjust the H Set
  auto allEdges = make_sparse_table<std::pair<uintE,uintE>,bool,hash_pair>(2 * edges.size() + 1,std::make_tuple(std::make_pair(UINT_E_MAX,UINT_E_MAX),false),hash_pair());
  par_for(0,edges.size(),[&](size_t i) {
    allEdges.insert(std::make_tuple(edges[i],true));
    allEdges.insert(std::make_tuple(std::make_pair(edges[i].second,edges[i].first),true));
  });
  sequence<uintE> originalH = hset.allH();


  hset.insertEdges(edges);

  sequence<uintE> hs = hset.allH();

  auto originalHSet = make_sparse_table<uintE,bool,hash_uintE>(originalH.size() * 2 + 1,std::make_tuple(UINT_E_MAX,false),hash_uintE());

  par_for(0,originalH.size(),[&](size_t i) {
    originalHSet.insert(std::make_tuple(originalH[i],true));
  });  

  adjustHSetWedges(originalH,hs,allEdges,originalHSet);

  // triangles added by wedges
  pbbslib::dyn_arr<uintE> wT = pbbslib::dyn_arr<uintE>(edges.size());
  // triangles added not by h-set
  pbbslib::dyn_arr<uintE> aT = pbbslib::dyn_arr<uintE>(edges.size());

  par_for(0,edges.size(),[&](size_t i) {
    wT.A[i] = 0;
    aT.A[i] = 0;
  });

  // newly formed wedges

  // STEP 2 - Find all the triangles/wedges by either wedges or entities formed solely by the newly added edges
  pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>> tW = pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>>(2 * edges.size());

  par_for(0,edges.size(),[&](size_t i) {
    uintE u = std::get<1>(edges[i]);
    uintE v = std::get<0>(edges[i]);

    tW.A[2 * i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);
    tW.A[2 * i + 1] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);

    // Case 1, added by wedges
    wT.A[i] = wedges.find(std::make_pair(v,u),0);

    if(!hset.inH(v)) {
      auto edgesV = dsg.v_data.A[v].neighbors.entries();

      pbbslib::dyn_arr<uintE> aTv = pbbslib::dyn_arr<uintE>(edgesV.size());

      pbbslib::dyn_arr<std::pair<uintE,uintE>> tWv = pbbslib::dyn_arr<std::pair<uintE,uintE>>(2 * edgesV.size());
      par_for(0,edgesV.size(),[&](size_t i) {
        aTv.A[i] = 0;
        uintE next = std::get<0>(edgesV[i]);
        if(next == u || !std::get<1>(edgesV[i])) {
          tWv.A[2 * i] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
          tWv.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
          return;
        }
        if(allEdges.find(std::make_pair(v,next),false)) {
          // Both are added edges

          tWv.A[2 * i] = std::make_pair(u,next);
          tWv.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
        }else{
          tWv.A[2 * i] = std::make_pair(u,next);
          tWv.A[2 * i + 1] = std::make_pair(next,u);
        }

        if(dsg.existEdge(v,next) && dsg.existEdge(next,u)) {
          // A triangle!
          if(allEdges.find(std::make_pair(v,next),false) && !(allEdges.find(std::make_pair(u,next),false))) {
            // (u,v) and (v,next) are added edges but not (next,u).
            // Possible cases are 2 and 3
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 2, 2 added edges and no h-set
              aTv.A[i] = 3;
            }else if(hset.inH(v) && !hset.inH(u) && !hset.inH(next)){
              // Case 3, 2 added edges and joint node is in h-set
              aTv.A[i] = 6;
            }
          }else if(allEdges.find(std::make_pair(u,next),false) && !(allEdges.find(std::make_pair(v,next),false))) {
            // (u,v) and (u,next) are added edges but not (v,next).
            // same with above
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 2, 2 added edges and no h-set
              aTv.A[i] = 3;
            }else if(hset.inH(u) && !hset.inH(v) && !hset.inH(next)){
              // Case 3, 2 added edges and joint node is in h-set
              aTv.A[i] = 6;
            }
          }else if(allEdges.find(std::make_pair(u,next),false) && allEdges.find(std::make_pair(v,next),false)) {
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 8, 3 added edges and no nodes in h-set
              aTv.A[i] = 2;
            }
          }
        }
      });
      tWv.size = tWv.capacity;

      tW.A[2 * i] = seq2da(filter(tWv.to_seq(),[&](std::pair<uintE,uintE> cur) {if((cur.first) == UINT_E_MAX && (cur.second) == UINT_E_MAX) {return false;} return true;}));
      aTv.size = aTv.capacity;
      aT.A[i] = pbbs::reduce(aTv.to_seq(),pbbs::addm<uintE>());
    }

    if(!hset.inH(u)) {
      auto edgesU = dsg.v_data.A[u].neighbors.entries();


      pbbslib::dyn_arr<uintE> aTu = pbbslib::dyn_arr<uintE>(edgesU.size());

      // newly formed wedges (one existing, one added)
      pbbslib::dyn_arr<std::pair<uintE,uintE>> tWu = pbbslib::dyn_arr<std::pair<uintE,uintE>>(2 * edgesU.size());

      par_for(0,edgesU.size(),[&](size_t i) {
        aTu.A[i] = 0;
        uintE next = std::get<0>(edgesU[i]);
        if(next == v || !std::get<1>(edgesU[i])) {

          tWu.A[2 * i] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
          tWu.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
          return;
        }
        if(allEdges.find(std::make_pair(u,next),false)) {
          // Both are added edges
          tWu.A[2 * i] = std::make_pair(v,next);
          tWu.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
        }else{
          tWu.A[2 * i] = std::make_pair(v,next);
          tWu.A[2 * i + 1] = std::make_pair(next,v);
        }

        if(dsg.existEdge(next,u) && dsg.existEdge(next,v)) {
          // A triangle!
          if(allEdges.find(std::make_pair(v,next),false) && !(allEdges.find(std::make_pair(u,next),false))) {
            // (u,v) and (v,next) are added edges but not (next,u).
            // Possible cases are 2 and 3
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 2, 2 added edges and no h-set
              aTu.A[i] = 3;
            }else if(hset.inH(v) && !hset.inH(u) && !hset.inH(next)){
              // Case 3, 2 added edges and joint node is in h-set
              aTu.A[i] = 6;
            }
          }else if(allEdges.find(std::make_pair(u,next),false) && !(allEdges.find(std::make_pair(v,next),false))) {
            // (u,v) and (u,next) are added edges but not (v,next)/
            // same with above

            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 2, 2 added edges and no h-set
              aTu.A[i] = 3;
            }else if(hset.inH(u) && !hset.inH(v) && !hset.inH(next)){
              // Case 3, 2 added edges and joint node is in h-set
              aTu.A[i] = 6;
            }
          }else if(allEdges.find(std::make_pair(u,next),false) && allEdges.find(std::make_pair(v,next),false)) {
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 8, 3 added edges and no nodes in h-set
              aTu.A[i] = 2;
            }
          }
        }
      });

      tWu.size = tWu.capacity;

      tW.A[2 * i + 1] = seq2da(filter(tWu.to_seq(),[&](std::pair<uintE,uintE> cur) {if((cur.first) == UINT_E_MAX && (cur.second) == UINT_E_MAX) {return false;} return true;}));
      aTu.size = aTu.capacity;
      aT.A[i] += pbbs::reduce(aTu.to_seq(),pbbs::addm<uintE>());
    }
  });


  // STEP 3 - Find triangles formed with the H Set
  // triangles formed with h-set
  pbbslib::dyn_arr<uintE> hT = pbbslib::dyn_arr<uintE>(edges.size());
  par_for(0,edges.size(),[&](size_t i) {
    hT.A[i] = 0;
    pbbslib::dyn_arr<uintE> uvhT = pbbslib::dyn_arr<uintE>(hs.size());
    int u = edges[i].first;
    int v = edges[i].second;
    par_for(0,hs.size(),[&](size_t j) {
      uvhT.A[j] = 0;
      int next = hs[j];
      if(u == next || v == next) return;
      if(dsg.existEdge(u,next) && dsg.existEdge(v,next)) {
        // Triangle!
        // I decided to this instead, since if i do case work, there would be like 26 if statements
        int counter = 1;
        if(allEdges.find(std::make_pair(next,u),false) && hset.inH(v)) {
          counter++;
        }
        if(allEdges.find(std::make_pair(v,next),false) && hset.inH(u)) {
          counter++;
        }
        uvhT.A[j] = 6 / counter;
      }
    });
    uvhT.size = uvhT.capacity;
    hT.A[i] = pbbs::reduce(uvhT.to_seq(),pbbs::addm<uintE>());
  });

  // STEP 4, update the counts and the wedges

  tW.size = tW.capacity;
  wT.size = wT.capacity;
  total += pbbs::reduce(wT.to_seq(),pbbs::addm<uintE>());

  aT.size = aT.capacity;
  total += pbbs::reduce(aT.to_seq(),pbbs::addm<uintE>()) / 12;

  hT.size = hT.capacity;
  total += pbbs::reduce(hT.to_seq(),pbbs::addm<uintE>()) / 6;

  // STEP 4

  pbbslib::dyn_arr<std::pair<uintE,uintE>> finalAdd = concatSeq(tW);

  sequence<std::pair<uintE,uintE>> all = merge_sort(finalAdd.to_seq(),[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a > b;});
  addWedges(all);

}

// Removes edges from the graph, which will also modify the H-Set and update the triangle counts
void removeEdges(sequence<std::pair<uintE,uintE>> edges) {
  auto allEdges = make_sparse_table<std::pair<uintE,uintE>,bool,hash_pair>(2 * edges.size() + 1,std::make_tuple(std::make_pair(UINT_E_MAX,UINT_E_MAX),false),hash_pair());
  par_for(0,edges.size(),[&](size_t i) {
    allEdges.insert(std::make_tuple(edges[i],true));
    allEdges.insert(std::make_tuple(std::make_pair(edges[i].second,edges[i].first),true));
  });

  sequence<uintE> originalH = hset.allH();

  // triangles removed by wedges
  pbbslib::dyn_arr<uintE> wT = pbbslib::dyn_arr<uintE>(edges.size());
  // triangles removed not by h-set
  pbbslib::dyn_arr<uintE> aT = pbbslib::dyn_arr<uintE>(edges.size());

  par_for(0,edges.size(),[&](size_t i) {
    wT.A[i] = 0;
    aT.A[i] = 0;
  });

  // wedges that need to be removed
  // STEP 1 - Finds all the triangles/wedges erased due to either the wedges or by the removed edges themselves
  pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>> tW = pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>>(2 * edges.size());

  par_for(0,edges.size(),[&](size_t i) {
    uintE u = std::get<1>(edges[i]);
    uintE v = std::get<0>(edges[i]);

    tW.A[2 * i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);
    tW.A[2 * i + 1] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);

    // Case 1, removed by wedges
    wT.A[i] = wedges.find(std::make_pair(v,u),0);

    if(!hset.inH(v)) {
      auto edgesV = dsg.v_data.A[v].neighbors.entries();

      pbbslib::dyn_arr<uintE> aTv = pbbslib::dyn_arr<uintE>(edgesV.size());

      pbbslib::dyn_arr<std::pair<uintE,uintE>> tWv = pbbslib::dyn_arr<std::pair<uintE,uintE>>(2 * edgesV.size());
      par_for(0,edgesV.size(),[&](size_t i) {
        aTv.A[i] = 0;
        uintE next = std::get<0>(edgesV[i]);
        if(next == u || !std::get<1>(edgesV[i])) {
          tWv.A[2 * i] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
          tWv.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
          return;
        }
        if(allEdges.find(std::make_pair(v,next),false)) {
          // Both are removed edges
          tWv.A[2 * i] = std::make_pair(u,next);
          tWv.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
        }else{
          tWv.A[2 * i] = std::make_pair(u,next);
          tWv.A[2 * i + 1] = std::make_pair(next,u);
        }

        if(dsg.existEdge(v,next) && dsg.existEdge(next,u)) {
          // A triangle!
          if(allEdges.find(std::make_pair(v,next),false) && !(allEdges.find(std::make_pair(u,next),false))) {
            // (u,v) and (v,next) are added edges but not (next,u).
            // Possible cases are 2 and 3
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 2, 2 added edges and no h-set
              aTv.A[i] = 3;
            }else if(hset.inH(v) && !hset.inH(u) && !hset.inH(next)){
              // Case 3, 2 added edges and joint node is in h-set
              aTv.A[i] = 6;
            }
          }else if(allEdges.find(std::make_pair(u,next),false) && !(allEdges.find(std::make_pair(v,next),false))) {
            // (u,v) and (u,next) are added edges but not (v,next).
            // same with above
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 2, 2 added edges and no h-set
              aTv.A[i] = 3;
            }else if(hset.inH(u) && !hset.inH(v) && !hset.inH(next)){
              // Case 3, 2 added edges and joint node is in h-set
              aTv.A[i] = 6;
            }
          }else if(allEdges.find(std::make_pair(u,next),false) && allEdges.find(std::make_pair(v,next),false)) {
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 8, 3 added edges and no nodes in h-set
              aTv.A[i] = 4;
            }
          }
        }
      });
      tWv.size = tWv.capacity;

      tW.A[2 * i] = seq2da(filter(tWv.to_seq(),[&](std::pair<uintE,uintE> cur) {if((cur.first) == UINT_E_MAX && (cur.second) == UINT_E_MAX) {return false;} return true;}));
      aTv.size = aTv.capacity;
      aT.A[i] = pbbs::reduce(aTv.to_seq(),pbbs::addm<uintE>());
    }

    if(!hset.inH(u)) {
      auto edgesU = dsg.v_data.A[u].neighbors.entries();


      pbbslib::dyn_arr<uintE> aTu = pbbslib::dyn_arr<uintE>(edgesU.size());

      // newly formed wedges (one existing, one added)
      pbbslib::dyn_arr<std::pair<uintE,uintE>> tWu = pbbslib::dyn_arr<std::pair<uintE,uintE>>(2 * edgesU.size());

      par_for(0,edgesU.size(),[&](size_t i) {
        aTu.A[i] = 0;
        uintE next = std::get<0>(edgesU[i]);
        if(next == v || !std::get<1>(edgesU[i])) {
          tWu.A[2 * i] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
          tWu.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
          return;
        }
        if(allEdges.find(std::make_pair(u,next),false)) {
          // Both are added edges
          tWu.A[2 * i] = std::make_pair(v,next);
          tWu.A[2 * i + 1] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
        }else{
          tWu.A[2 * i] = std::make_pair(v,next);
          tWu.A[2 * i + 1] = std::make_pair(next,v);
        }

        if(dsg.existEdge(v,next) && dsg.existEdge(next,u)) {
          // A triangle!
          if(allEdges.find(std::make_pair(v,next),false) && !(allEdges.find(std::make_pair(u,next),false))) {
            // (u,v) and (v,next) are added edges but not (next,u).
            // Possible cases are 2 and 3
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 2, 2 added edges and no h-set
              aTu.A[i] = 3;
            }else if(hset.inH(v) && !hset.inH(u) && !hset.inH(next)){
              // Case 3, 2 added edges and joint node is in h-set
              aTu.A[i] = 6;
            }
          }else if(allEdges.find(std::make_pair(u,next),false) && !(allEdges.find(std::make_pair(v,next),false))) {
            // (u,v) and (u,next) are added edges but not (v,next).
            // same with above
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 2, 2 added edges and no h-set
              aTu.A[i] = 3;
            }else if(hset.inH(u) && !hset.inH(v) && !hset.inH(next)){
              // Case 3, 2 added edges and joint node is in h-set
              aTu.A[i] = 6;
            }
          }else if(allEdges.find(std::make_pair(u,next),false) && allEdges.find(std::make_pair(v,next),false)) {
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 8, 3 added edges and no nodes in h-set
              aTu.A[i] = 4;
            }
          }
        }
      });

      tWu.size = tWu.capacity;

      tW.A[2 * i + 1] = seq2da(filter(tWu.to_seq(),[&](std::pair<uintE,uintE> cur) {if((cur.first) == UINT_E_MAX && (cur.second) == UINT_E_MAX) {return false;} return true;}));
      aTu.size = aTu.capacity;
      aT.A[i] += pbbs::reduce(aTu.to_seq(),pbbs::addm<uintE>());
    }
  });


  // STEP 2 - triangles removed with h-set
  pbbslib::dyn_arr<long long> hT = pbbslib::dyn_arr<long long>(edges.size());
  par_for(0,edges.size(),[&](size_t i) {
    hT.A[i] = 0;
    pbbslib::dyn_arr<long long> uvhT = pbbslib::dyn_arr<long long>(originalH.size());
    int u = edges[i].first;
    int v = edges[i].second;
    par_for(0,originalH.size(),[&](size_t j) {
      uvhT.A[j] = 0;
      int next = originalH[j];
      if(dsg.existEdge(u,next) && dsg.existEdge(v,next)) {
        // Triangle!
        // I decided to this instead, since if i do case work, there would be like 26 if statements
        // TODO: What is the meaning of these counters?
        int counter = 1;
        int counter2 = -1;
        if(allEdges.find(std::make_pair(u,next),false) && hset.inH(v)) {
          counter++;
        }
        if(allEdges.find(std::make_pair(v,next),false) && hset.inH(u)) {
          counter++;
        }
        if(allEdges.find(std::make_pair(u,next),false) && !hset.inH(v)) {
          counter2++;
        }
        if(allEdges.find(std::make_pair(v,next),false) && !hset.inH(u)) {
          counter2++;
        }
        uvhT.A[j] = 6 * counter2 / counter;
      }
    });
    uvhT.size = uvhT.capacity;
    hT.A[i] = pbbs::reduce(uvhT.to_seq(),pbbs::addm<long long>());
  });

  // STEP 3 - Updates all the counts

  tW.size = tW.capacity;
  wT.size = wT.capacity;
  total -= pbbs::reduce(wT.to_seq(),pbbs::addm<uintE>());
  aT.size = aT.capacity;
  total += pbbs::reduce(aT.to_seq(),pbbs::addm<uintE>()) / 12;
  hT.size = hT.capacity;
  total += pbbs::reduce(hT.to_seq(),pbbs::addm<long long>()) / 6;

  // STEP 4 - Updates the HSet and the wedges

  hset.eraseEdges(edges);

  sequence<uintE> hs = hset.allH();

  auto originalHSet = make_sparse_table<uintE,bool,hash_uintE>(originalH.size() * 2 + 1,std::make_tuple(UINT_E_MAX,false),hash_uintE());

  par_for(0,originalH.size(),[&](size_t i) {
    originalHSet.insert(std::make_tuple(originalH[i],true));
  });

  pbbslib::dyn_arr<std::pair<uintE,uintE>> finalRemove = concatSeq(tW);
  sequence<std::pair<uintE,uintE>> allR = merge_sort(finalRemove.to_seq(),[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a > b;});
  removeWedges(allR);
  adjustHSetWedges(originalH,hs,allEdges,originalHSet); 
}
};

