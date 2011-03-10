#include <iostream>
#include <set>
#include <list>
#include <map>
#include <vector>
#include <iterator>
#include <ctime>
#include <cstdlib>
#include "boost/tr1/unordered_map.hpp"
#include "boost/tr1/unordered_set.hpp"

#include "printable.hpp"

using std::string;
using std::list;
using std::cout;
using std::endl;
using std::set;
using std::vector;

using std::tr1::unordered_set;
using std::tr1::unordered_map;

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


/**
   Returns mutants
*/
vector<string> mut(const string g) { 
  vector<string> result;
  result.reserve((alphabet.size() -1) * g.length() );
  const size_t geno_length = g.length(); // go right-to-left to generate in-order
  for(int pos = geno_length - 1; pos != -1; --pos) { 
    string::iterator alphabet_end;
    for(string::iterator alternative = alphabet.begin(), alphabet_end = alphabet.end();
	alternative != alphabet_end; ++alternative) {
      if (*alternative != g[pos]) {
	string mutant(g);
	mutant[pos] = *alternative;
	result.push_back(mutant);
	//	cout << "generated mutant " << mutant << endl;
      }
    }
  }
  return result;
}


/* generates random int phenotype */
inline
pheno colorOf(const geno & g) {
  return rand() % numOfColors;
}

/* s1 - keys(m)
   
   @param[in] s1 a mathematical set of Ts
   @param[in] m an unordered_map with keys of T

 */
template <class T, class TVal>
vector<T> set_difference(const vector<T> & s1, 
			 const unordered_map<T,TVal> & m) {
  vector<T> result;
  for(typename vector<T>::const_iterator it_add = s1.begin(); it_add != s1.end(); ++it_add) {
    if (m.find(*it_add) == m.end() ) {
      result.push_back(*it_add);
    }
  }
  return result;
}

/**
   Searches breadth-first from ROOT, recording the color of every node
   observed, but traversing onward only from nodes with the same color.

   Nodes x and y are connected if x is in mut(y). Each node x has
   color colorOf(x).

   Returns a std::map of the observed nodes and their colors.
*/
unordered_map<geno,pheno> search(const geno& root) {
  unordered_map<geno,pheno> observed;
  observed[root]=CLUSTER_COLOR;

  list<geno> to_traverse;
  to_traverse.push_back(root);

  geno cursor;
  while( !to_traverse.empty() ) {
    // process head of queue of nodes till to be traversed
    cursor = to_traverse.front();
    to_traverse.pop_front();
//    cout << "Traversing from node " << cursor << endl;
    // compute the neighbors not previously observed
    vector<geno> all_neighbors(mut(cursor));
    vector<geno> new_neighbors(set_difference(all_neighbors, observed));
//    cout << "\tFound " << cursor << " had neighbors " << all_neighbors << " of which the previously unobserved were " << new_neighbors << endl;
    // if there are some new nodes to observe ...
    if(!new_neighbors.empty()) {
      unordered_map<geno,pheno> newly_observed;
      for(vector<geno>::iterator g = new_neighbors.begin(); 
	  g != new_neighbors.end(); ++g) {
	// "discover" the colors and record the visit
	if ((observed[*g] = colorOf(*g)) == CLUSTER_COLOR) {
	  // and  plan to visit only the special ones later
	  to_traverse.push_back(*g);
	}
      }
    }
  }
  return observed;
}

int main()
{
  //  srand(time(NULL)); // seed the random number generator
  srand(0); // seed the random number generator
  geno g = createRootString(4);
  initialize_alphabet_size(10);
  initialize_numOfColors(3);

  /*
  cout << "mutants of " << g << " are: ";
  set<geno> s(mut(g));
  cout << s << endl;
  */

  // display one search
  if (false) {
    cout << "search starting at " << g << endl;
    unordered_map<geno,pheno> mm(search(g));
    std::map<geno,pheno> m;
    m.insert(mm.begin(),mm.end());
    cout << m << endl;
  }

  if (true) {
    // benchmark 10 random searches
    for(int i =0; i < 20; ++i) {
      search(g);
    }
  }
  return 0;
}
