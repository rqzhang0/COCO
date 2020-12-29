// shared_ptr.h
// Ben Gamari
// August 2009

#ifndef _GWU_SHARED_PTR_H
#define _GWU_SHARED_PTR_H

namespace qcd {

template<typename T>
struct shared_ptr {
private:
  T* ptr;
  int* ref_count;

public:
  shared_ptr(T* p) : ptr(p), ref_count(new int(1)) { }

  shared_ptr(const shared_ptr& s) : ptr(s.ptr), ref_count(s.ref_count) {
    add_ref();
  }

  shared_ptr& operator=(const shared_ptr& s) {
    remove_ref();
    this->ptr = s.ptr;
    this->ref_count = s.ref_count;
    add_ref();
    return *this;
  }

  ~shared_ptr() {
    remove_ref();
  }

  void add_ref() {
    (*ref_count)++;
  }

  void remove_ref() {
    (*ref_count)--;
    if (*ref_count == 0) {
      delete ptr;
      delete ref_count;
    }
  }

  T* operator->() const {
    return ptr;
  }

  T* get() const {
    return ptr;
  }

  T& operator*() const {
    return *ptr;
  }
};

} // qcd namespace

#endif


