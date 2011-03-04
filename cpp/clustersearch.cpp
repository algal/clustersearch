#include <iostream>
#include <set>
#include <list>
#include <map>
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

string alphabet = "ABC";

/**
   Generates mutants of a genotype.
*/
set<geno> mut(geno g) { 
  set<geno> result;

  list<set<geno>::iterator> ends;

  const size_t geno_length = g.length();
  for(int pos = geno_length - 1; pos != -1; --pos) { // go backwards to generate in-order
    string::iterator alphabet_end;
    for(string::iterator alternative = alphabet.begin(), alphabet_end = alphabet.end();
	alternative != alphabet_end; ++alternative) {
      if (*alternative != g[pos]) {
	string mutant(g);
	mutant[pos] = *alternative;
	result.insert(result.end(),mutant);
	//	cout << "generated mutant: " << mutant << endl;
      }
    }
  }

  for(list<set<geno>::iterator>::iterator it = ends.begin(); it != ends.end(); ++it) {
    if(*it != ends.front()) 
      cout << "nonequality detected" << endl;
  }
    
  return result;
}


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
      map<geno,pheno> newly_observed;
      for(list<geno>::iterator g = new_neighbors.begin(); g != new_neighbors.end(); ++g) {
	newly_observed[*g] = colorOf(*g);
	if (newly_observed[*g] == SPECIAL_COLOR) 
	  to_traverse.push_back(*g);
      }
      observed.insert(newly_observed.begin(),newly_observed.end());
    }
  }
  return observed;
}

/* prints a set<string> */
std::ostream & operator<< (std::ostream & out, set<string> & s) {
  out << "{";
  for(set<string>::iterator item_it = s.begin(); item_it != s.end(); ++item_it) {
    out << *item_it << " ";
  }
  out << "}";
  return out;
}      

int main()
{
  geno g = "AAAA";
  cout << "mutants of " << g << " are " << endl;
  set<geno> s(mut(g));
  cout << s << endl;
  return 0;
}
