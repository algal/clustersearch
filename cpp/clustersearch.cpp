#include <iostream>
#include <set>
#include <list>
#include <map>
#include <vector>
#include <iterator>
#include <ctime>
#include <cstdlib>
#include <tr1/unordered_map>
#include <tr1/unordered_set>

#include "printable.hpp"

using std::string;
using std::list;
using std::cout;
using std::endl;
using std::set;
using std::vector;

using std::tr1::unordered_map;
using std::tr1::unordered_set;

typedef string geno;
typedef unsigned char pheno;

// don't touch these globals, which define how all functions work
const pheno CLUSTER_COLOR = 0;
pheno numOfColors = 3;
string alphabet = "ABC"; // must be in lexicographical order
unsigned int length=0;

/** Initialize the alphabet to a different size */
void initialize_alphabet_size(unsigned int length) {
  const string max_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  if (length > 52) {
    std::cerr << "ERROR: maximum alphabet size is 52. defaulting to 52" << endl;
    ::length = 52;
  }
  ::alphabet = max_alphabet.substr(0,length);
}

/** Initialize the number of possible phenotypes */
void initialize_numOfColors(const unsigned int num) { 
  ::numOfColors = num; 
}

/** creates string with all chars equal, implicitly defining length of
    all string in the genotype space */
void initialize_length(size_t length) { 
  ::length = length;
}

/** creates string with all chars equal, implicitly defining length of
    all string in the genotype space */
string createRootString() { 
  return string(::length, alphabet[0]); 
}

/**
   Returns mutants
*/
vector<string> mut(const string g) { 
  vector<string> result;
  result.reserve((alphabet.size() -1) * ::length );
  for(int pos = ::length - 1; pos != -1; --pos) { 
    string::iterator alphabet_end;
    for(string::iterator alternative = alphabet.begin(), alphabet_end = alphabet.end();
	alternative != alphabet_end; ++alternative) {
      if (*alternative != g[pos]) {
	string mutant(g);
	mutant.reserve(::length); // may improve perf
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

/**
   Do a search run.
   
   Returns dictionary over all points in a cluster (i.e., neutral
   network) and its mutational neighborhood, 'picked' from random
   string graph with string LENGTH, built from an alphabet of
   ALPHABETSIZE, where every node is painted one color out of
   NUMOFCOLORS.
   
   'Picks' this cluster by doing a breadth-first search, assigning
   colors randomly as it progresses.
   
   NOTE:
   - modifies global variables alphabet and colors
   - alphabetsize and numOfColors must be < 26 

*/
unordered_map<geno,pheno> doRun(const unsigned int length, const unsigned int alphabetsize, const unsigned int numOfColors) {
  initialize_alphabet_size(alphabetsize);
  initialize_numOfColors(numOfColors);
  initialize_length(length);

  const geno origin = createRootString();
  return search(origin);
}

extern "C" struct cluster_measures {
  unsigned int cluster_size;
  unsigned int perimeter_size;
  unsigned int colors;
  unsigned int exits_size;
};

/**
   Calculates measures of a cluster.

   Relies on that the search only observes the cluster and its perimeter.
 */
cluster_measures calculate_measures_from_run(const unordered_map<geno,pheno> & m) {
  cluster_measures results;

  unordered_set<pheno> perimeter_colors;
  unsigned int exits_size=0;
  unsigned int perimeter_size =0;

  // for every item in the perimeter ...
  unordered_map<geno,pheno>::const_iterator end(m.end());
  for(unordered_map<geno,pheno>::const_iterator it = m.begin(); it != end; ++it) {
    const geno g = it->first;
    const pheno p = it->second;
    if(p != CLUSTER_COLOR) {
      // ... tally it, add its pheno to a unique set
      ++perimeter_size;
      perimeter_colors.insert(p);

      const vector<geno> mutants( mut(g) );
      vector<geno>::const_iterator end_mut( mutants.end() );
      for(vector<geno>::const_iterator it_mut = mutants.begin(); it_mut != end_mut; ++it_mut) {
	if(m.at(*it_mut) == CLUSTER_COLOR) {
	  ++exits_size;
	}
      }
    }
  }

  results.perimeter_size = perimeter_size;
  results.cluster_size = m.size() - results.perimeter_size;
  results.colors = perimeter_colors.size();
  results.exits_size = exits_size;

  return results;
}

extern "C" cluster_measures calculate_measures(const unsigned int length, const unsigned int alphabetsize, const unsigned int numOfColors) {
  srand(0);
  return calculate_measures_from_run(doRun(length,alphabetsize,numOfColors));
}



extern "C" size_t cluster_size(const unsigned int length, const unsigned int alphabetsize, const unsigned int numOfColors) {
  srand(0);
  unordered_map<geno,pheno> m(doRun(length,alphabetsize,numOfColors));
  return m.size();
}

int main()
{
  //  srand(time(NULL)); // seed the random number generator
  srand(0); // seed the random number generator

  /*
  cout << "mutants of " << g << " are: ";
  set<geno> s(mut(g));
  cout << s << endl;
  */

  // display one search
  if (false) {
    unordered_map<geno,pheno> mm(doRun(11,2,2));
    std::map<geno,pheno> m;
    m.insert(mm.begin(),mm.end());
    cout << m << endl;
  }

  if (true) {
    // benchmark 10 random searches
    for(int i =0; i < 1000; ++i) {
      doRun(11,2,2);
    }
  }
  return 0;
}
