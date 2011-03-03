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
//using std::inserter

typedef string geno;
typedef string pheno;

pheno colors[] = {"white", "red", "blue"};
/**
FIXME: implement
*/
set<geno> mut(geno g) { return set<geno>(); }

pheno colorOf(geno g) { return colors[0]; }

// iterators over map keys
typedef map<geno,pheno>::iterator map_iterator;
typedef map<geno,pheno>::key_type (*get_key_t)(map<geno,pheno>::value_type);
typedef boost::transform_iterator<get_key_t, map_iterator> key_iterator;
map<geno,pheno>::key_type get_key(map<geno,pheno>::value_type aPair) { return aPair.first; }

/**
   Searches breadth-first from root node, recording color of every
   node observed, but traversing onward only from nodes with the color

   Nodes x and y are connected if (x in mut(y)). Each node x has
   color colorOf(x).

   Returns a std::map of the observed nodes and their colors.
*/
map<geno,pheno> search(const geno& root) {
  const pheno SPECIAL_COLOR = colors[0];

  map<geno,pheno> observed;
  observed[root]=SPECIAL_COLOR;

  list<geno> to_traverse;
  to_traverse.push_back(root);

  geno cursor;
  while( !to_traverse.empty() ) {
    cursor = to_traverse.front();
    to_traverse.pop_front();
    set<geno> all_neighbors(mut(cursor));
    list<geno> new_neighbors;
    set_difference(all_neighbors.begin(),all_neighbors.end(),
		   key_iterator(observed.begin(), get_key), key_iterator(observed.end(), get_key),
		   std::inserter(new_neighbors, new_neighbors.end()));

    if(!new_neighbors.empty()) {
      map<geno,pheno> newnodes;
      for(list<geno>::iterator g = new_neighbors.begin(); g != new_neighbors.end(); ++g) {
	newnodes[*g] = colorOf(*g);
	if (newnodes[*g] == SPECIAL_COLOR) 
	  to_traverse.push_back(*g);
      }
      observed.insert(newnodes.begin(),newnodes.end());
    }
  }
  return observed;
}
      

int main()
{
  std::cout << "Hello, world" << std::endl;
}
