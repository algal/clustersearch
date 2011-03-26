#include <iostream>
#include <list>
#include <vector>
#include <iterator>
#include <ctime>
#include <cstdlib>
// #include <tr1/unordered_map>
// #include <tr1/unordered_set>
#include "boost/tr1/unordered_map.hpp"
#include "boost/tr1/unordered_set.hpp"

#include "boost/accumulators/accumulators.hpp"
#include "boost/accumulators/statistics/stats.hpp"
#include "boost/accumulators/statistics/mean.hpp"
#include "boost/accumulators/statistics/moment.hpp"

#include "printable.hpp"

//#define DEBUG

#ifdef DEBUG
#define TRACE(arg) (arg)
#else
#define TRACE(arg) ((void)0)
#endif

using std::string;
using std::list;
using std::cout;
using std::endl;
using std::vector;

using std::tr1::unordered_map;
using std::tr1::unordered_set;

typedef string geno;
typedef unsigned int pheno;

// don't touch these globals, which define how all functions work
unsigned int numOfColors = 3;
pheno CLUSTER_COLOR = 0;

size_t alphabet_size = 3;
string alphabet = "ABC"; // must be in lexicographical order

unsigned int length=0;

const double GRAY_UNUSED=0.0;
double gray_fraction=GRAY_UNUSED;
vector<double> cdf;

/** Initialize the alphabet to a different size */
void initialize_alphabet_size(unsigned int new_alpha_size) {
  const string max_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  if (new_alpha_size > 52) {
    std::cerr << "ERROR: maximum alphabet size is 52. defaulting to 52" << endl;
    new_alpha_size = 52;
  }
  ::alphabet = max_alphabet.substr(0,new_alpha_size);
  ::alphabet_size = ::alphabet.size();
}

/** Initialize the number of possible phenotypes */
void initialize_numOfColors(const unsigned int num) { 
  ::numOfColors = num; 
  // define: CLUSTER_COLOR=last in the list
  CLUSTER_COLOR = ::numOfColors-1;
}

/** creates string with all chars equal, implicitly defining length of
    all string in the genotype space */
void initialize_length(size_t length) { 
  ::length = length;
}

/** initialize gray_fraction
    
    @param g if g=0, then gray is not used as a color, and all colors
    are equally likely. If 0<g<1, then g is probability of a random
    point being gray, and all other colors are equally likely (i.e.,
    (1-g)/(numOfColors-1).

 */
void initialize_gray_fraction(const double g) {
  ::gray_fraction = g;

  if( ::gray_fraction == GRAY_UNUSED)
    return;

  // define: first color is gray
  // pdf[0] = gray_fraction
  // pdf[i] = (1-g)/(1-numOfCOlors) , when i!=0
  vector<double> pdf(::numOfColors, ((1-::gray_fraction) / (::numOfColors-1)) );
  pdf[0] = ::gray_fraction;

  vector<double> cdf;
  cdf.reserve(pdf.size());
  double running_sum = 0;
  for(vector<double>::iterator it = pdf.begin(); it != pdf.end(); ++it) {
    running_sum += *it;
    cdf.push_back(running_sum);
  }
  ::cdf = cdf;
}

/**
   Returns mutants
*/
vector<string> mut(const string g) { 
  vector<string> result;
  result.reserve((::alphabet_size -1) * ::length );
  for(int pos = ::length - 1; pos != -1; --pos) { 
    string::iterator alphabet_end;
    for(string::iterator alternative = alphabet.begin(), alphabet_end = alphabet.end();
	alternative != alphabet_end; ++alternative) {
      if (*alternative != g[pos]) {
	string mutant(g);
	mutant.reserve(::length); // may improve perf
	mutant[pos] = *alternative;
	result.push_back(mutant);
	TRACE(	cout << "generated mutant " << mutant << endl);
      }
    }
  }
  return result;
}


/** Generates random phenotype.

    If g == GRAY_UNUSED, then all phenotypes are uniformly distributed. 

    Otherwise, then the "gray" phenotype has probability g, and all
    the others are uniformly distributed.
*/
inline
pheno colorOf(const geno & g) {
  if( ::gray_fraction == GRAY_UNUSED )
    return rand() % ::numOfColors;
  else {
    const double real = (double)rand() / (RAND_MAX + 1);
    for(unsigned int i = 0; i < ::numOfColors; ++i) {
      if( real < ::cdf[i] )
	return i;
    }
  }
  // assert(ERROR)
  return ::numOfColors-1;
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

/* intersection of s1 and keys(m)
   
   @param[in] s1 a mathematical set of Ts
   @param[in] m an unordered_map with keys of T

*/
template <class T, class TVal>
vector<T> set_intersection(const vector<T> & s1,
			   const unordered_map<T,TVal> & m) {
  vector<T> result;
  for(typename vector<T>::const_iterator it_add = s1.begin(); it_add != s1.end(); ++it_add) {
    if (m.find(*it_add) != m.end() ) {
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

   Returns a std::unordered_map of the observed nodes and their colors.
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
    TRACE(cout << "Traversing from node " << cursor << endl);
    // compute the neighbors not previously observed
    vector<geno> all_neighbors(mut(cursor));
    vector<geno> new_neighbors(set_difference(all_neighbors, observed));
    TRACE(cout << "\tFound " << cursor << " had neighbors " << all_neighbors << " of which the previously unobserved were " << new_neighbors << endl);
    // if there are some new nodes to observe ...
    if(!new_neighbors.empty()) {
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

  const geno origin = string(::length, alphabet[0]);
  return search(origin);
}

extern "C"
void reseed(unsigned int seed) {
  std::srand(seed);
}

extern "C" 
struct cluster_measures {
  unsigned int cluster_size;	// s
  unsigned int perimeter_size;	// t
  unsigned int colors;		// E
  unsigned int exits_size;	// u
  double robustness;            // r
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

  // for every item ...
  for(unordered_map<geno,pheno>::const_iterator it = m.begin(); it != m.end(); ++it) {
    const geno g = it->first;
    const pheno p = it->second;
    TRACE(cout << "searching g=" << g << endl);
    // .. in the perimeter ...
    if(p != CLUSTER_COLOR) {
      TRACE(cout << "\twhich is in the perimeter" << endl);
      // ... tally it, track its pheno
      ++perimeter_size;
      perimeter_colors.insert(p);
      TRACE(cout << "\tadded its color" << endl);
      vector<geno> mutants( set_intersection( mut(g),m) );
      TRACE(cout << "\tcalculated its (observed) mutants: " << mutants << endl);
      for(vector<geno>::iterator it_mut = mutants.begin(); it_mut != mutants.end(); ++it_mut) {
	if(m.find(*it_mut)->second == CLUSTER_COLOR) {
	  ++exits_size;
	}
      }
    }
  }

  results.perimeter_size = perimeter_size;
  results.cluster_size = m.size() - results.perimeter_size;
  results.colors = perimeter_colors.size();
  results.exits_size = exits_size;
  const double sl  = (results.cluster_size * ::length);
  results.robustness = double(sl - exits_size) / sl;
  return results;
}

extern "C"
cluster_measures calculate_measures(const unsigned int length, const unsigned int alphabetsize, const unsigned int numOfColors) {
  return calculate_measures_from_run(doRun(length,alphabetsize,numOfColors));
}

extern "C"
struct mean_cluster_measures {
  double mean_cluster_size;	// s
  double mean_perimeter_size;	// t
  double mean_colors;	        // E
  double mean_exits_size;	// u
  double mean_robustness;       // r
};

// calculate means of the cluster measures
mean_cluster_measures calculate_statistics(const unsigned int length, 
					   const unsigned int alphabetsize, 
					   const unsigned int numOfColors,
					   const unsigned int samples) {
  using namespace boost::accumulators;

  accumulator_set<unsigned int, stats<tag::mean> > cluster_size_acc;
  accumulator_set<unsigned int, stats<tag::mean> > perimeter_size_acc;
  accumulator_set<unsigned int, stats<tag::mean> > colors_acc;
  accumulator_set<unsigned int, stats<tag::mean> > exits_size_acc;
  accumulator_set<double,       stats<tag::mean> > robustness_acc;
  
  cluster_measures r;
  for(unsigned int i = 0; i < samples; ++i ) {
    r = calculate_measures(length,alphabetsize,numOfColors);
    cluster_size_acc(r.cluster_size);
    perimeter_size_acc(r.perimeter_size);
    colors_acc(r.colors);
    exits_size_acc(r.exits_size);
    robustness_acc(r.robustness);
  }  

  mean_cluster_measures result =
    {
      mean(cluster_size_acc),
      mean(perimeter_size_acc),
      mean(colors_acc),
      mean(exits_size_acc),
      mean(robustness_acc)
    };
  return  result;
}


string usage() {
  return 
    "Usage: clusters LENGTH ALPHABETSIZE COLORS [[SAMPLES] SILENT]\n"
    "Prints measures from one search with given parameters, or else means over samples searches"
    "\n"
    "\n"
    "Arguments: \n"
    " length        length of genotype strings \n"
    " alphabetsize  possible symbols in position of a genotype string \n"
    " colors        number of possible 'phenotype' colors \n"
    " samples       number of searches to perform \n"
    " silent        only print out numbers \n";
}

int main(int argc, char *argv[])
{
  unsigned int seed;
  unsigned int length;
  unsigned int alphabetsize;
  unsigned int numOfColors;

  unsigned int samples =0;

  bool silent = false;
  seed = 0;
  
  if(argc==4 || argc==5 || argc==6) {
    if(argc ==6)
      silent = true;

    if(argc == 5 || argc == 6)
      samples           = std::atoi(argv[4]);
    
    length		= std::atoi(argv[1]);
    alphabetsize	= std::atoi(argv[2]);
    numOfColors		= std::atoi(argv[3]);
  }
  else {
    std::cerr << usage() << endl;
    exit(1);
  }

  //  srand(time(NULL)); // seed the random number generator
  std::srand(seed); // seed the random number generator

  const bool NORMAL_EXECUTION = true;
  if(NORMAL_EXECUTION) {

    if(!silent) {
      cout << "searching with:" << endl;
      cout << "\tlength = " << length << endl;
      cout << "\talphabetsize = " << alphabetsize << endl;
      cout << "\tnumOfColors = " << numOfColors << endl;
    }

    // display one search
    if (samples == 0) {

      unordered_map<geno,pheno> mm(doRun(length,alphabetsize,numOfColors));
      cout << mm << endl;
    
      cluster_measures results = calculate_measures_from_run(mm);
      cout << "cluster_size		= " << results.cluster_size << endl;
      cout << "results.perimeter_size	= " << results.perimeter_size << endl;
      cout << "results.colors		= " << results.colors << endl;
      cout << "results.exits_size	= " << results.exits_size << endl;
      cout << "results.robustness	= " << results.robustness << endl;
    }
    else {
      if(!silent)
	cout << "\tsamples = " << samples << endl;

      mean_cluster_measures results = calculate_statistics(length,alphabetsize,numOfColors,samples);
    
      if(!silent) {
	cout << "mean cluster_size		= " << results.mean_cluster_size << endl;
	cout << "mean results.perimeter_size	= " << results.mean_perimeter_size << endl;
	cout << "mean results.colors		= " << results.mean_colors << endl;
	cout << "mean results.exits_size		= " << results.mean_exits_size << endl;
	cout << "mean results.robustness		= " << results.mean_robustness << endl;
      }
      else {
	cout << results.mean_cluster_size << endl;
	cout << results.mean_perimeter_size << endl;
	cout << results.mean_colors << endl;
	cout << results.mean_exits_size << endl;
	cout << results.mean_robustness << endl;
      }
    }
  }
  
  // benchmark 10 random searches
  if (false) {
    for(int i =0; i < 10; ++i) {
      doRun(4,4,5);
    }
  }

  // check 2nd call of random
  if (false) {
    for(int i =0; i < 3; ++i) {
      std::srand(seed); // seed the random number generator
      (void) doRun(10,4,5);
      const unordered_map<geno,pheno> result1 = doRun(10,4,5);
      std::srand(seed); // seed the random number generator
      (void) doRun(10,4,5);
      const unordered_map<geno,pheno> result2 = doRun(10,4,5);
      if( result1 != result2) 
	cout << "doRun() identical on 1st call after re-seeding" << endl;
      else
	cout << "doRun() NOT identical on 1st call after re-seeding" << endl;
    }
  }

  return 0;
}
