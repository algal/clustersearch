#include <iostream>
#include <map>
#include <iterator>
#include "boost/iterator/transform_iterator.hpp"

using std::map;
typedef std::string Key;
typedef std::string Val;

map<Key,Val>::key_type get_key(map<Key,Val>::value_type aPair) {
  return aPair.first;
}

typedef map<Key,Val>::key_type (*get_key_t)(map<Key,Val>::value_type);
typedef map<Key,Val>::iterator map_iterator;
typedef boost::transform_iterator<get_key_t, map_iterator> mapkey_iterator;

int main() {
  map<Key,Val> m;
  m["a"]="A";
  m["b"]="B";
  m["c"]="C";

  // iterate over the map's (key,val) pairs as usual
  for(map_iterator i = m.begin(); i != m.end(); i++) {
    std::cout << i->first << " " << i->second << std::endl;
  }

  // iterate over the keys using the transformed iterators
  mapkey_iterator keybegin(m.begin(), get_key);
  mapkey_iterator keyend(m.end(), get_key);
  for(mapkey_iterator i = keybegin; i != keyend; i++) {
    std::cout << *i << std::endl;
  }
}
