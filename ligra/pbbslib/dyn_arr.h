// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010-2016 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// A simple dynamic array implementation.
//
// Note: Have only tested this for simple datatypes, in particular, it uses newA
// (i.e. malloc) to resize, which will cause cryptic segfaults if E requires a
// constructor call to properly initialize memory. see pbbslib::new_array(..)
#pragma once

#include "ligra/bridge.h"

namespace pbbslib {

  constexpr size_t kDynArrMinBktSize = 2000;

  template <class E>
  struct dyn_arr {
    E* A;
    size_t size;
    size_t capacity;
    bool alloc;

    dyn_arr() : A(NULL), size(0), capacity(0), alloc(false) {}
    dyn_arr(size_t s) : size(0), capacity(s), alloc(true) { A = pbbslib::new_array_no_init<E>(s); }
    dyn_arr(E* _A, long _size, long _capacity, bool _alloc)
        : A(_A), size(_size), capacity(_capacity), alloc(_alloc) {}

    //Copy constructor
    dyn_arr(const dyn_arr& array) {
      A = pbbslib::new_array_no_init<E>(array.capacity);
      par_for(0, array.size, [&] (size_t i) { A[i] = array.A[i]; });
      size = array.size;
      capacity = array.capacity;
    }

    //Copy assignment operator
    dyn_arr& operator=(const dyn_arr& array) {
      A = pbbslib::new_array_no_init<E>(array.capacity);
      par_for(0, array.size, [&] (size_t i) { A[i] = array.A[i]; });
      size = array.size;
      capacity = array.capacity;
      return *this;
    }

    void del() {
      if (alloc) {
        pbbslib::free_array(A);
        alloc = false;
      }
    }

    pbbs::sequence<E> to_seq() {
      assert(A);
      auto ret = pbbs::sequence<E>(A, size);
      size = 0;
      A = nullptr;
      return std::move(ret);
    }

    void clear() { size = 0; }

    inline void resize(size_t n) {
      if (n + size > capacity) {
        size_t new_capacity = std::max(2 * (n + size), (size_t)kDynArrMinBktSize);
        E* nA = pbbslib::new_array_no_init<E>(new_capacity);
        par_for(0, size, 2000, [&] (size_t i) { nA[i] = A[i]; });
        if (alloc) {
          pbbslib::free_array(A);
        }
        A = nA;
        capacity = new_capacity;
        alloc = true;
      }
    }

    inline void insert(E val, size_t pos) { A[size + pos] = val; }

    inline void push_back(E val) {
      A[size] = val;
      size++;
    }

    // TODO: What do you mean by "does not maintain order"?
    // TODO: You have to be careful about using this kind of stuff when you
    // want to remove an element from dyn_arr, because it will leave
    // gaps in your array that you have to deal with in your code.
    // This is lazy deletion. Probably, it's better to have a wrapper
    // around dyn_arr that's for dynamic insertion/deletion,
    // that will keep an "empty" value (as sparse table does), and
    // then lazily delete by putting in empty values. It should then
    // try and keep track of how many deletions create holes in your 
    // array, and at some bound go in and compress everything 
    // to be in order again. When it retrieves elements,
    // it'll go to some upper limit count that is given based off of 
    // how many holes we've allowed through deletions, and it'll
    // filter out empty elements.
    // TODO: Also, you may want to put in a batch erase so it's parallel.
    //Erase does not maintain order (does it have to?)
    inline void erase(E val) {
      size_t index = indexOf(val);

      if (index < size - 1) {
        // TODO: Why would you replace with the last element?
        A[index] = A[size - 1]; //Replaces value with last element
      }
      size--; //Decreases size
    }

    //Similar to sparse_table implementation
    inline size_t indexOf(E val) {
      // TODO: Don't use int. Use size_t. We'll exceed int, and note
      // that size is of type size_t -- always good practice to match types.
      for (int i = 0; i < size; i++) {
        if (A[i] == val) return i;
      }
      // TODO: Generally, our vertex neighbor lists will be sorted.
      // If this is what you're looking up into, you should not have
      // this kind of indexOf function -- instead, use a binary search.
      // (We have binary search implementations).
      return -1;
    }

    inline void add(E val) {
      if (size >= capacity) {
        resize(1); //Should double size
      }
      push_back(val);
    }

    template <class F>
    void map(F f) {
      par_for(0, size, 2000, [&] (size_t i) { f(A[i]); });
    }

    template <class F>
    inline void copyIn(F& f, size_t n) {
      resize(n);
      par_for(0, n, 2000, [&] (size_t i) { A[size + i] = f[i]; });
      size += n;
    }

    template <class F>
    inline void copyInF(F f, size_t n) {
      resize(n);
      par_for(0, n, 2000, [&] (size_t i) { A[size + i] = f(i); });
      size += n;
    }
  };
}; // namespace pbbslib
