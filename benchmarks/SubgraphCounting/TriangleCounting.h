

uintE N = 10;

dynamic_symmetric_graph<dynamic_symmetric_vertex,uintE> dsg = createEmptyDynamicSymmetricGraph<dynamic_symmetric_vertex,uintE>();
HSetAlex hset = HSetAlex();
uintE wedges[10][10];
uintE total = 0;

auto seq2da = [&](sequence<std::pair<uintE,uintE>> s) -> pbbslib::dyn_arr<std::pair<uintE,uintE>> {pbbslib::dyn_arr<std::pair<uintE,uintE>> res = pbbslib::dyn_arr<std::pair<uintE,uintE>>(s.size()); par_for(0,s.size(),[&](size_t j) {res.A[j] = s[j];}); res.size = s.size(); return res;};

struct hash_pair {
  inline size_t operator () (const std::pair<uintE,uintE> & a) {return ((pbbs::hash64(a.first) * 3) ^ (pbbs::hash64(a.second) >> 32));}
};


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


void initialize(size_t n) {
  hset.resizeV(n);
  sequence<uintE> nodes = sequence<uintE>(n,[&](size_t i) {return i;});
  sequence<std::pair<uintE,uintE>> hPairs = sequence<std::pair<uintE,uintE>>(n,[&](size_t i) {return std::make_pair(i,0);});
  dsg.batchAddVertices(nodes);
  hset.insert(hPairs);
}

void addEdges(sequence<std::pair<uintE,uintE>> edges) {
  auto allEdges = make_sparse_table<std::pair<uintE,uintE>,bool,hash_pair>(2 * edges.size() + 1,std::make_tuple(std::make_pair(UINT_E_MAX,UINT_E_MAX),false),hash_pair());
  par_for(0,edges.size(),[&](size_t i) {
    allEdges.insert(std::make_tuple(edges[i],true));
    allEdges.insert(std::make_tuple(std::make_pair(edges[i].second,edges[i].first),true));
  });

  auto aN = merge_sort(sequence<uintE>(2 * edges.size(),[&](size_t i) {if(i % 2) {
    return edges[i / 2].first;
  }else{
    return edges[i / 2].second;
  }}),[&](uintE a,uintE b){return a > b;});


  auto allNodes = sequence<uintE>(2 * edges.size(),[&](size_t i) {
    if(i == 0) {
      return aN[i];
    }else{
      if(aN[i] != aN[i - 1]) {
        return aN[i];
      }else{
        return UINT_E_MAX;
      }
    }
  });

  allNodes = filter(allNodes,[&](uintE i) {return (i != UINT_E_MAX);});


  // Update H-Set first or work might become > h^2. 

  // Find all the nodes that weren't in H-Set but now are. Also nodes that were in H-Set and now aren't
  
  dsg.batchAddEdges(edges);

  sequence<uintE> originalH = hset.allH();


  auto tmp = sequence<std::pair<uintE,uintE>>(allNodes.size(),[&](size_t i) {return std::make_pair(allNodes[i],dsg.v_data.A[allNodes[i]].degree);});
  // par_for(0,allNodes.size(),[&](size_t i) {cout << allNodes[i] << "   " << dsg.v_data.A[allNodes[i]].degree << endl;});

  hset.modify(sequence<std::pair<uintE,uintE>>(allNodes.size(),[&](size_t i) {return std::make_pair(allNodes[i],dsg.v_data.A[allNodes[i]].degree);}));

  sequence<uintE> hs = hset.allH();

  auto originalHSet = make_sparse_table<uintE,bool,hash_uintE>(originalH.size() * 2 + 1,std::make_tuple(UINT_E_MAX,false),hash_uintE());

  par_for(0,originalH.size(),[&](size_t i) {
    originalHSet.insert(std::make_tuple(originalH[i],true));
  });


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





  // triangles added by wedges
  pbbslib::dyn_arr<uintE> wT = pbbslib::dyn_arr<uintE>(edges.size());
  // triangles added not by h-set
  pbbslib::dyn_arr<uintE> aT = pbbslib::dyn_arr<uintE>(edges.size());

  par_for(0,edges.size(),[&](size_t i) {
    wT.A[i] = 0;
    aT.A[i] = 0;
  });


  // newly formed wedges
  pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>> tW = pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>>(2 * edges.size());

  par_for(0,edges.size(),[&](size_t i) {
    uintE u = std::get<1>(edges[i]);
    uintE v = std::get<0>(edges[i]);

    tW.A[2 * i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);
    tW.A[2 * i + 1] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);

    // Case 1, added by wedges
    wT.A[i] = wedges[v][u];

    if(!hset.inH(v)) {
      auto edgesV = dsg.v_data.A[v].neighbors.entries();

      pbbslib::dyn_arr<uintE> aTv = pbbslib::dyn_arr<uintE>(edgesV.size());

      pbbslib::dyn_arr<std::pair<uintE,uintE>> tWv = pbbslib::dyn_arr<std::pair<uintE,uintE>>(2 * edgesV.size());
      par_for(0,edgesV.size(),[&](size_t i) {
        aTv.A[i] = 0;
        uintE next = std::get<0>(edgesV[i]);
        if(next == u) {
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
            // (u,v) and (u,next) are added edges but not (v,next)/
            // same with above
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 2, 2 added edges and no h-set
              aTv.A[i] = 3;
            }else if(hset.inH(v) && !hset.inH(u) && !hset.inH(next)){
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
        if(next == v) {

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

        if(dsg.existEdge(next,u)) {
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
            }else if(hset.inH(v) && !hset.inH(u) && !hset.inH(next)){
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
      if(dsg.existEdge(u,next) && dsg.existEdge(v,next)) {
        // Triangle!
        // I decided to this instead, since if i do case work, there would be like 26 if statements
        int counter = 1;
        if(allEdges.find(std::make_pair(u,next),false) && hset.inH(v)) {
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

  tW.size = tW.capacity;
  // cout << "lalalala" << endl;
  wT.size = wT.capacity;
  // cout << pbbs::reduce(wT.to_seq(),pbbs::addm<uintE>()) << endl;
  total += pbbs::reduce(wT.to_seq(),pbbs::addm<uintE>());

  aT.size = aT.capacity;
  total += pbbs::reduce(aT.to_seq(),pbbs::addm<uintE>()) / 12;

  hT.size = hT.capacity;
  total += pbbs::reduce(hT.to_seq(),pbbs::addm<uintE>()) / 6;


  // for(int i = 0;i < dW.size;++i) {
  //   cout << dW.A[i].size << endl;
  // }


  pbbslib::dyn_arr<std::pair<uintE,uintE>> finalAdd = concatSeq(tW);


  pbbslib::dyn_arr<std::pair<uintE,uintE>> fdW = concatSeq(dW);

  finalAdd = concat2Seq(finalAdd,fdW);

  sequence<std::pair<uintE,uintE>> all = merge_sort(finalAdd.to_seq(),[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a > b;});
  

  par_for(0,all.size(),[&](size_t i) {
    if(i != 0) {
      if(all[i] != all[i - 1]) {
        wedges[all[i - 1].first][all[i - 1].second] += i;
      }
    }
  });

  par_for(0,all.size(),[&](size_t i) {
    if(i != 0) {
      if(all[i] != all[i - 1]) {
        wedges[all[i].first][all[i].second] -= i;
      }
    }
  });

 if(all.size() != 0) wedges[all[all.size() - 1].first][all[all.size() - 1].second] += all.size();


  pbbslib::dyn_arr<std::pair<uintE,uintE>> finalRemove = concatSeq(rW);
  sequence<std::pair<uintE,uintE>> allR = merge_sort(finalRemove.to_seq(),[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a > b;});

  par_for(0,allR.size(),[&](size_t i) {
    if(i != 0) {
      if(allR[i] != allR[i - 1]) {
        wedges[allR[i - 1].first][allR[i - 1].second] -= i;
      }
    }
  });

  par_for(0,allR.size(),[&](size_t i) {
    if(i != 0) {
      if(allR[i] != allR[i - 1]) {
        wedges[allR[i].first][allR[i].second] += i;
      }
    }
  });
  if(allR.size() != 0) wedges[allR[allR.size() - 1].first][allR[allR.size() - 1].second] -= allR.size();

}

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
  pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>> tW = pbbslib::dyn_arr<pbbslib::dyn_arr<std::pair<uintE,uintE>>>(2 * edges.size());


  par_for(0,edges.size(),[&](size_t i) {
    uintE u = std::get<1>(edges[i]);
    uintE v = std::get<0>(edges[i]);

    tW.A[2 * i] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);
    tW.A[2 * i + 1] = pbbslib::dyn_arr<std::pair<uintE,uintE>>(0);

    // Case 1, removed by wedges
    wT.A[i] = wedges[v][u];

    if(!hset.inH(v)) {
      auto edgesV = dsg.v_data.A[v].neighbors.entries();

      pbbslib::dyn_arr<uintE> aTv = pbbslib::dyn_arr<uintE>(edgesV.size());

      pbbslib::dyn_arr<std::pair<uintE,uintE>> tWv = pbbslib::dyn_arr<std::pair<uintE,uintE>>(2 * edgesV.size());
      par_for(0,edgesV.size(),[&](size_t i) {
        aTv.A[i] = 0;
        uintE next = std::get<0>(edgesV[i]);
        if(next == u) {
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
            // (u,v) and (u,next) are added edges but not (v,next)/
            // same with above
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 2, 2 added edges and no h-set
              aTv.A[i] = 3;
            }else if(hset.inH(v) && !hset.inH(u) && !hset.inH(next)){
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
        if(next == v) {

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
            // (u,v) and (u,next) are added edges but not (v,next)/
            // same with above
            if(!hset.inH(u) && !hset.inH(v) && !hset.inH(next)) {
              // Case 2, 2 added edges and no h-set
              aTu.A[i] = 3;
            }else if(hset.inH(v) && !hset.inH(u) && !hset.inH(next)){
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



  // triangles removed with h-set
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
  tW.size = tW.capacity;
  // cout << "lalalala" << endl;
  wT.size = wT.capacity;
  total -= pbbs::reduce(wT.to_seq(),pbbs::addm<uintE>());
  aT.size = aT.capacity;
  total += pbbs::reduce(aT.to_seq(),pbbs::addm<uintE>()) / 12;

  hT.size = hT.capacity;
  total += pbbs::reduce(hT.to_seq(),pbbs::addm<long long>()) / 6;

  auto aN = merge_sort(sequence<uintE>(2 * edges.size(),[&](size_t i) {if(i % 2) {
    return edges[i / 2].first;
  }else{
    return edges[i / 2].second;
  }}),[&](uintE a,uintE b){return a > b;});

  auto allNodes = sequence<uintE>(2 * edges.size(),[&](size_t i) {
    if(i == 0) {
      return aN[i];
    }else{
      if(aN[i] != aN[i - 1]) {
        return aN[i];
      }else{
        return UINT_E_MAX;
      }
    }
  });

  allNodes = filter(allNodes,[&](uintE i) {return (i != UINT_E_MAX);});
  dsg.batchRemoveEdges(edges);

  auto tmp = sequence<std::pair<uintE,uintE>>(allNodes.size(),[&](size_t i) {return std::make_pair(allNodes[i],dsg.v_data.A[allNodes[i]].degree);});
  // par_for(0,allNodes.size(),[&](size_t i) {cout << allNodes[i] << "   " << dsg.v_data.A[allNodes[i]].degree << endl;});

  hset.modify(sequence<std::pair<uintE,uintE>>(allNodes.size(),[&](size_t i) {return std::make_pair(allNodes[i],dsg.v_data.A[allNodes[i]].degree);}));

  sequence<uintE> hs = hset.allH();

  auto originalHSet = make_sparse_table<uintE,bool,hash_uintE>(originalH.size() * 2 + 1,std::make_tuple(UINT_E_MAX,false),hash_uintE());

  par_for(0,originalH.size(),[&](size_t i) {
    originalHSet.insert(std::make_tuple(originalH[i],true));
  });


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
            dW.A[i].A[j * (specialSet.size() - 1) + k] = std::make_pair(specialSet[j],specialSet[k]);
          }else if(k > j) {
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


  // Added edges due to being removed from h-set
  pbbslib::dyn_arr<std::pair<uintE,uintE>> fdW = concatSeq(dW);


  sequence<std::pair<uintE,uintE>> all = merge_sort(fdW.to_seq(),[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a > b;});

  par_for(0,all.size(),[&](size_t i) {
    if(i != 0) {
      if(all[i] != all[i - 1]) {
        wedges[all[i - 1].first][all[i - 1].second] += i;
      }
    }
  });

  par_for(0,all.size(),[&](size_t i) {
    if(i != 0) {
      if(all[i] != all[i - 1]) {
        wedges[all[i].first][all[i].second] -= i;
      }
    }
  });

  if(all.size() != 0) wedges[all[all.size() - 1].first][all[all.size() - 1].second] += all.size();

  pbbslib::dyn_arr<std::pair<uintE,uintE>> frW = concatSeq(rW);
  pbbslib::dyn_arr<std::pair<uintE,uintE>> finalRemove = concatSeq(tW);
  finalRemove = concat2Seq(finalRemove,frW);

  sequence<std::pair<uintE,uintE>> allR = merge_sort(finalRemove.to_seq(),[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a > b;});

  par_for(0,allR.size(),[&](size_t i) {
    if(i != 0) {
      if(allR[i] != allR[i - 1]) {
        wedges[allR[i - 1].first][allR[i - 1].second] -= i;
      }
    }
  });

  par_for(0,allR.size(),[&](size_t i) {
    if(i != 0) {
      if(allR[i] != allR[i - 1]) {
        wedges[allR[i].first][allR[i].second] += i;
      }
    }
  });

  if(allR.size() != 0) wedges[allR[allR.size() - 1].first][allR[allR.size() - 1].second] -= allR.size();

}

// // Basic concept and not parallel
// uintE connectEdge(uintE v,uintE u) {
//   if(dsg.existEdge(u,v)) return total;
//   hset.connect(u,v);

//   auto entries = hset.P.entries();
//   total += wedges[v][u];
//   par_for(0,entries.size(),1,[&](size_t i){
//     if(dsg.existEdge(std::get<0>(entries[0]),u) && dsg.existEdge(std::get<0>(entries[0]),v)) ++total;
//   });
//   if(!checkTF(hset.P,v)) {
//     par_for(0,dsg.v_data[v].entries.size,1,[&](size_t i) {
//       ++wedges[u][std::get<0>(dsg.v_data[v].entries[i])];
//       ++wedges[std::get<0>(dsg.v_data[v].entries[i])][u];
//     });
//   }

//   if(!checkTF(hset.P,u)) {
//     par_for(0,dsg.v_data[u].entries.size,1,[&](size_t i) {
//       ++wedges[v][std::get<0>(dsg.v_data[u].entries[i])];
//       ++wedges[std::get<0>(dsg.v_data[u].entries[i])][v];
//     });
//   }
//   dsg.addEdge(u,v);
//   return total;
// }