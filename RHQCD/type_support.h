// type_support.h
// A. Alexandru
// July 2014

#pragma once

namespace qcd{

template<bool, typename T = void> 
  struct enable_if {};

template<typename T>
  struct enable_if<true, T> 
  {
    typedef T type;
  };

template <typename T>
  struct is_numeric {const static bool val = false;};

template <> struct is_numeric<double> {const static bool val = true;};
template <> struct is_numeric<float > {const static bool val = true;};
template <> struct is_numeric<int   > {const static bool val = true;};

template <typename T1, typename T2>
  struct max_type {};

template <typename T> struct max_type<T, T>       {typedef T res;};
template <typename T> struct max_type<double, T>  {typedef double res;};
template <typename T> struct max_type<T, double>  {typedef double res;};
template <>           struct max_type<float, int> {typedef float res;};
template <>           struct max_type<int, float> {typedef float res;};

template <typename T>
  struct is_basictype {const static bool val = is_numeric<T>::val;};

template <> struct is_basictype<bool> {const static bool val=true;};
template <> struct is_basictype<char> {const static bool val=true;};
} // qcd namespace
