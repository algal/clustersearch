#include <iostream>
#include <list>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>


using std::map;
using std::string;
using std::list;
using std::cout;
using std::endl;
using std::set;
#using std::inserter

typedef string geno;
typedef string pheno;

// key iterators on maps


typedef map<geno,pheno>::iterator map_iterator;
typedef map<geno,pheno>::value_type

typedef boost::function
typedef map::<geno,pheno>::key_type (*F)(map<geno,pheno>::value_type);
typedef boost::transform_iterator<F,map_iterator> mapkey_iterator;


pheno range[] = {"white", "red", "blue"};

set<geno> mut(geno g) {
  return set<geno>;
}

set<geno> keys(map<geno> m) {
  set<geno> result;
  
}


map<geno,pheno> search(const geno& root) {
  const pheno SPECIAL_COLOR = range[0];

  map<geno,pheno> observed;
  observed[geno]=SPECIAL_COLOR;

  list<geno> to_traverse;
  to_traverse.push_back(root);

  geno cursor;
  while( !to_traverse.empty() ) {
    cursor = to_traverse.front();
    to_traverse.pop_front();
    set<geno> all_neighbors(mut(cursor));
    set<geno> all_observed(keys(observed));
    set<geno> new_neighbors;
    set_difference(all_neighbors.begin(),all_neighbors.end()
		   all_observed.begin(),all_observed.end(),
		   std::inserter(new_neighbors, new_neighbors.end()));

    if(!new_neighbors.empty()) {
      
    


  }

}

