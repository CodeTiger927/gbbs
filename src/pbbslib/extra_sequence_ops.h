// This file is a bridge connecting the "lib interface" gbbs exports and the
// interfact that the current pbbslib exports. We would like to support both
// C++11 users, and the current (C++17) implementation of the lib. Using this
// bridge will hopefully simplify having two separate implementations of the lib
// interface.

#pragma once

#include "pbbslib/seq.h"
#include "pbbslib/sequence_ops.h"
#include "pbbslib/macros.h"
#include "pbbslib/monoid.h"

// C++17 bridge
namespace pbbs {

  // used so second template argument can be inferred
  template <class T, class F>
  inline delayed_sequence<T,F> make_sequence (size_t n, F f) {
    return delayed_sequence<T,F>(n,f);
  } 

  // Scans the input sequence using the addm monoid.
  template <class In_Seq>
  inline auto scan_add_inplace(In_Seq const& In) -> typename In_Seq::value_type {
    using T = typename In_Seq::value_type;
    return pbbslib::scan_inplace(s, pbbslib::addm<T>());
  }

}

namespace pbbs {

  template <class Seq>
  inline auto reduce_max(Seq I, flags fl = no_flag) -> typename Seq::T {
    using T = typename Seq::T;
    return reduce(I, maxm<T>(), fl);
  }

  template <class Seq>
  inline auto reduce_min(Seq I, flags fl = no_flag) -> typename Seq::T {
    using T = typename Seq::T;
    return reduce(I, minm<T>(), fl);
  }

  template <class Seq>
  inline auto reduce_xor(Seq I, flags fl = no_flag) -> typename Seq::T {
    using T = typename Seq::T;
    return reduce(I, xorm<T>(), fl);
  }


  template <class Idx_Type, class D, class F>
  inline sequence<std::tuple<Idx_Type, D> > pack_index_and_data(
      F& f, size_t size, flags fl = no_flag) {
    auto identity = [&](size_t i) {
      return std::make_tuple((Idx_Type)i, std::get<1>(f(i)));
    };
    auto flgs_f = [&](size_t i) { return std::get<0>(f(i)); };
    auto flgs_in =
        make_sequence<bool>(size, flgs_f);
    return pack(make_sequence<std::tuple<Idx_Type, D> >(size, identity), flgs_in,
                fl);
  }

  template <class T, class Pred>
  inline size_t filter_seq(T* in, T* out, size_t n, Pred p) {
    size_t k = 0;
    for (size_t i = 0; i < n; i++)
      if (p(in[i])) out[k++] = in[i];
    return k;
  }

  // Faster for a small number in output (about 40% or less)
  // Destroys the input.   Does not need a bool array.
  template <class T, class PRED>
  inline size_t filterf(T* In, T* Out, size_t n, PRED p) {
    size_t b = _F_BSIZE;
    if (n < b) return filter_seq(In, Out, n, p);
    size_t l = num_blocks(n, b);
    size_t* Sums = new_array_no_init<size_t>(l + 1);
    par_for(0, l, 1, [&] (size_t i) {
      size_t s = i * b;
      size_t e = std::min(s + b, n);
      size_t k = s;
      for (size_t j = s; j < e; j++) {
        if (p(In[j])) In[k++] = In[j];
      }
      Sums[i] = k - s;
    });
    auto isums = sequence<size_t>(Sums, l);
    size_t m = scan_add_inplace(isums);
    Sums[l] = m;
    par_for(0, l, 1, [&] (size_t i) {
      T* I = In + i * b;
      T* O = Out + Sums[i];
      for (size_t j = 0; j < Sums[i + 1] - Sums[i]; j++) {
        O[j] = I[j];
      }
    });
    pbbslib::free_array(Sums);
    return m;
  }


  // Faster for a small number in output (about 40% or less)
  // Destroys the input.   Does not need a bool array.
  template <class T, class PRED, class OUT>
  inline size_t filterf(T* In, size_t n, PRED p, OUT out, size_t out_off) {
    size_t b = _F_BSIZE;
    if (n < b) {
      size_t k = out_off;
      for (size_t i = 0; i < n; i++) {
        if (p(In[i])) out(k++, In[i]);
      }
      return k - out_off;
    }
    size_t l = num_blocks(n, b);
    size_t* Sums = new_array_no_init<size_t>(l + 1);
    par_for(0, l, 1, [&] (size_t i) {
      size_t s = i * b;
      size_t e = std::min(s + b, n);
      size_t k = s;
      for (size_t j = s; j < e; j++) {
        if (p(In[j])) In[k++] = In[j];
      }
      Sums[i] = k - s;
    });
    auto isums = sequence<size_t>(Sums, l);
    size_t m = scan_add(isums, isums);
    Sums[l] = m;
    par_for(0, l, 1, [&] (size_t i) {
      T* I = In + i * b;
      size_t si = out_off + Sums[i];
      for (size_t j = 0; j < Sums[i + 1] - Sums[i]; j++) {
        out(si + j, I[j]);
      }
    });
    pbbslib::free_array(Sums);
    return m;
  }

  template <class T, class PRED>
  inline size_t filterf_and_clear(T* In, T* Out, size_t n, PRED p, T& empty) {
    size_t b = _F_BSIZE;
    if (n < b) {
      size_t ret = filter_seq(In, Out, n, p);
      for (size_t i=0; i<n; i++) {
        if (p(In[i])) {
          In[i] = empty;
        }
      }
      return ret;
    }
    size_t l = num_blocks(n, b);
    b = num_blocks(n, l);
    size_t* Sums = new_array_no_init<size_t>(l + 1);

    par_for(0, l, 1, [&] (size_t i) {
      size_t s = i * b;
      size_t e = std::min(s + b, n);
      size_t k = s;
      for (size_t j = s; j < e; j++) {
        if (p(In[j])) {
          In[k] = In[j];
          if (k != j) {
            In[j] = empty;
          }
          k++;
        }
      }
      Sums[i] = k - s;
    });
    auto isums = sequence<size_t>(Sums, l);
    size_t m = scan_add(isums, isums);
    Sums[l] = m;
    par_for(0, l, 1, [&] (size_t i) {
      T* I = In + (i * b);
      size_t i_off = Sums[i];
      size_t num_i = Sums[i+1] - i_off;
      T* O = Out + i_off;
      for (size_t j = 0; j < num_i; j++) {
        O[j] = I[j];
        I[j] = empty;
      }
    });
    pbbslib::free_array(Sums);
    return m;
  }

}
