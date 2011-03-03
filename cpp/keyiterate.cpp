#include <iostream>
#include <list>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>

#include "boost/iterator/transform_iterator.hpp"

using std::map;
using std::string;
using std::list;
using std::cout;
using std::endl;
using std::set;

typedef string geno;
typedef string pheno;

map<geno,pheno>::key_type get_key(map<geno,pheno>::value_type aPair) {
  return aPair.first;
}

//typedef boost::function<map<geno,pheno>::key_type, map<geno,pheno>::value_type> func_type;
//typedef map<geno,pheno>::key_type (*func_t)(map<geno,pheno>::value_type);
typedef map<geno,pheno>::iterator map_iterator;
typedef boost::transform_iterator<get_key, map_iterator> mapkey_iterator;

int main() {
  map<geno,pheno> m;
  m["a"]="A";
  m["b"]="B";
  m["c"]="C";

  mapkey_iterator keybegin(m.begin(), func);
  mapkey_iterator keyend(m.end(), func);

  typedef map<geno,pheno>::value_type;

  typedef boost::function;
  typedef map::<geno,pheno>::key_type (*F)(map<geno,pheno>::value_type);
  typedef boost::transform_iterator<F,map_iterator> mapkey_iterator;

}
