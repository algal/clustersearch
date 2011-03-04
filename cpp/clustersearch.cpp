#include <iostream>
#include <set>
#include <list>
#include <map>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <cstdlib>
#include "boost/iterator/transform_iterator.hpp"

using std::map;
using std::string;
using std::list;
using std::cout;
using std::endl;
using std::set;

typedef string geno;
typedef unsigned int pheno;

// don't touch these globals, which define how all functions work
string alphabet = "ABC"; // must be in lexicographical order
size_t numOfColors = 3;
const pheno CLUSTER_COLOR = 0;

/** Initialize the alphabet to a different size */
void initialize_alphabet_size(unsigned int length) {
  const string max_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  if (length > 52) {
    std::cerr << "ERROR: maximum alphabet size is 52. defaulting to 52" << endl;
    length=52;
  }
  alphabet = max_alphabet.substr(0,length);
}

/** Initialize the number of possible phenotypes */
void initialize_numOfColors(const unsigned int num) { 
  numOfColors = num; 
}

/** creates string with all chars equal, implicitly defining length of
    all string in the genotype space */
string createRootString(size_t length) { 
  return string(length, alphabet[0]); 
}
  

/* prints a list<T> */
template <class T>
std::ostream & operator<<(std::ostream & out, list<T> & s) {
  out << "[";
  for(typename list<T>::iterator item_it = s.begin(); item_it != s.end(); ++item_it) {
    out << *item_it << " ";
  }
  out << "]";
  return out;
}      

/* prints a set<T> */
template <class T>
std::ostream & operator<<(std::ostream & out, set<T> & s) {
  out << "{";
  for(typename set<T>::iterator item_it = s.begin(); item_it != s.end(); ++item_it) {
    out << *item_it << " ";
  }
  out << "}";
  return out;
}      

/* prints a map<TKey,TVal> */
template <class TKey, class TVal>
std::ostream & operator<<(std::ostream & out, map<TKey,TVal> & m) {
  out << "{";
  for(typename map<TKey,TVal>::iterator item_it = m.begin(); item_it != m.end(); ++item_it) {
    out << item_it->first << ": " << item_it->second << ", ";
  }
  out << "}";
  return out;
}      

/**
   Returns mutants in lexicographical order.
*/
set<string> mut(const string g) { 
  set<string> result;
  const size_t geno_length = g.length(); // go right-to-left to generate in-order
  for(int pos = geno_length - 1; pos != -1; --pos) { 
    string::iterator alphabet_end;
    for(string::iterator alternative = alphabet.begin(), alphabet_end = alphabet.end();
	alternative != alphabet_end; ++alternative) {
      if (*alternative != g[pos]) {
	string mutant(g);
	mutant[pos] = *alternative;
	result.insert(result.end(),mutant);
	//	cout << "generated mutant " << mutant << endl;
      }
    }
  }
  return result;
}


/* generates random int phenotype */
inline
pheno colorOf(const geno g) {
  return rand() % numOfColors;
}

// iterators over map keys
typedef map<geno,pheno>::iterator map_iterator;
typedef map<geno,pheno>::key_type (*get_key_t)(map<geno,pheno>::value_type);
typedef boost::transform_iterator<get_key_t, map_iterator> key_iterator;
map<geno,pheno>::key_type get_key(const map<geno,pheno>::value_type aPair) { return aPair.first; }

/**
   Searches breadth-first from root node, recording color of every
   node observed, but traversing onward only from nodes with the color

   Nodes x and y are connected if (x in mut(y)). Each node x has
   color colorOf(x).

   Returns a std::map of the observed nodes and their colors.
*/
map<geno,pheno> search(const geno& root) {
  map<geno,pheno> observed;
  observed[root]=CLUSTER_COLOR;

  list<geno> to_traverse;
  to_traverse.push_back(root);

  geno cursor;
  while( !to_traverse.empty() ) {
    // process head of queue of nodes till to be traversed
    cursor = to_traverse.front();
    to_traverse.pop_front();
    cout << "Traversing from node " << cursor << endl;
    // compute the neighbors not previously observed
    set<geno> all_neighbors(mut(cursor));
    list<geno> new_neighbors;
    set_difference(all_neighbors.begin(),all_neighbors.end(),
		   key_iterator(observed.begin(), get_key), key_iterator(observed.end(), get_key),
		   std::inserter(new_neighbors, new_neighbors.end()));
    cout << "\tFound " << cursor << " had neighbors " << all_neighbors << " of which the previously unobserved were " << new_neighbors << endl;
    // if there are some new nodes to observe ...
    if(!new_neighbors.empty()) {
      map<geno,pheno> newly_observed;
      for(list<geno>::iterator g = new_neighbors.begin(); g != new_neighbors.end(); ++g) {
	// "discover" the colors 
	newly_observed[*g] = colorOf(*g);
	// plan to visit only the special ones later
	if (newly_observed[*g] == CLUSTER_COLOR) {
	  to_traverse.push_back(*g);
	}
      }
      cout << "\tObserved these nodes to be colored: " << newly_observed << endl;
      // add them to the db of observed nodes
      observed.insert(newly_observed.begin(),newly_observed.end());
    }
  }
  return observed;
}

int main()
{
  srand(time(NULL)); // seed the random number generator
  geno g = createRootString(4);
  initialize_alphabet_size(10);
  initialize_numOfColors(3);


  /*
  cout << "mutants of " << g << " are: ";
  set<geno> s(mut(g));
  cout << s << endl;
  */

  cout << "search starting at " << g << endl;
  map<geno,pheno> m = search(g);
  cout << m << endl;

  // for(int i =0; i < 10; ++i) {
  //   map<geno,pheno> m = search(g);
  //   //    mut(g);
  // }

  return 0;
}
