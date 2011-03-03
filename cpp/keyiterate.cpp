#include <iostream>
#include <map>
#include <iterator>

#include "boost/iterator/transform_iterator.hpp"

using std::map;
using std::string;
using std::cout;
using std::endl;

typedef string geno;
typedef string pheno;

map<geno,pheno>::key_type get_key(map<geno,pheno>::value_type aPair) {
  return aPair.first;
}

typedef map<geno,pheno>::key_type (*get_key_t)(map<geno,pheno>::value_type);
typedef map<geno,pheno>::iterator map_iterator;
typedef boost::transform_iterator<get_key_t, map_iterator> mapkey_iterator;

int main() {
  map<geno,pheno> m;
  m["a"]="A";
  m["b"]="B";
  m["c"]="C";

  mapkey_iterator keybegin(m.begin(), get_key);
  mapkey_iterator keyend(m.end(), get_key);

  for(map_iterator i = m.begin(); i != m.end(); i++) {
    cout << i->first << " " << i->second << endl;
  }

  for(mapkey_iterator i = keybegin; i != keyend; i++) {
    cout << *i << endl;
  }
}
