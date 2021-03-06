#ifndef PRINTABLE_H
#define PRINTABLE_H

#include <iostream>
#include <set>
#include <list>
#include <map>
#include <vector>
// #include <tr1/unordered_map>
// #include <tr1/unordered_set>
#include "boost/tr1/unordered_map.hpp"
#include "boost/tr1/unordered_set.hpp"


/* prints a list<T> */
template <class T>
std::ostream & operator<<(std::ostream & out, std::list<T> & s) {
  out << "[";
  for(typename std::list<T>::iterator item_it = s.begin(); item_it != s.end(); ++item_it) {
    out << *item_it << " ";
  }
  out << "]";
  return out;
}      

/* prints a unordered_set<T> */
template <class T>
std::ostream & operator<<(std::ostream & out, std::tr1::unordered_set<T> & s) {
  out << "{";
  for(typename std::tr1::unordered_set<T>::iterator item_it = s.begin(); item_it != s.end(); ++item_it) {
    out << *item_it << " ";
  }
  out << "}";
  return out;
}      

/* prints a pair<T1, T2> */
template <class T1, class T2>
std::ostream & operator<<(std::ostream & out, std::pair<T1,T2> & s) {
  out << "(" << s.first  << " " << s.second  << ")";
  return out;
}      

/* prints a vector<T> */
template <class T>
std::ostream & operator<<(std::ostream & out, std::vector<T> & s) {
  out << "{";
  for(typename std::vector<T>::iterator item_it = s.begin(); item_it != s.end(); ++item_it) {
    out << *item_it << " ";
  }
  out << "}";
  return out;
}      

/* prints a set<T> */
template <class T>
std::ostream & operator<<(std::ostream & out, std::set<T> & s) {
  out << "{";
  for(typename std::set<T>::iterator item_it = s.begin(); item_it != s.end(); ++item_it) {
    out << *item_it << " ";
  }
  out << "}";
  return out;
}      

/* prints a map<TKey,TVal> */
template <class TKey, class TVal>
std::ostream & operator<<(std::ostream & out, std::map<TKey,TVal> & m) {
  out << "{";
  for(typename std::map<TKey,TVal>::iterator item_it = m.begin(); item_it != m.end(); ++item_it) {
    out << item_it->first << ": " << item_it->second << ", ";
  }
  out << "}";
  return out;
}      

/* prints a unordered_map<TKey,TVal> */
template <class TKey, class TVal>
std::ostream & operator<<(std::ostream & out, std::tr1::unordered_map<TKey,TVal> & mm) {
  std::map<TKey,TVal> m;
  m.insert(mm.begin(),mm.end());
  out << m;
  return out;
}      

#endif
