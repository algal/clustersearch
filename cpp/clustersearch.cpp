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


#include "boost/program_options.hpp"


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

namespace configs {
  // don't touch these globals, which define how all functions work
  unsigned int numOfColors = 3;
  pheno CLUSTER_COLOR = 0;
  
  size_t alphabet_size = 3;
  string alphabet = "ABC"; // must be in lexicographical order
  
  unsigned int length=0;
  
  const double GRAY_UNUSED=0.0;
  double gray_fraction=GRAY_UNUSED;
  vector<double> cdf;
}


/** Initialize the alphabet to a different size */
void initialize_alphabet_size(unsigned int new_alpha_size) {
  const string max_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  if (new_alpha_size > 52) {
    std::cerr << "ERROR: maximum alphabet size is 52. defaulting to 52" << endl;
    new_alpha_size = 52;
  }
  configs::alphabet = max_alphabet.substr(0,new_alpha_size);
  configs::alphabet_size = configs::alphabet.size();
}

/** creates string with all chars equal, implicitly defining length of
    all string in the genotype space */
void initialize_length(size_t length) { 
  configs::length = length;
}

/** Initialize the number of possible phenotypes */
void initialize_numOfColors(const unsigned int num) { 
  configs::numOfColors = num; 
  // define: CLUSTER_COLOR=last in the list
  configs::CLUSTER_COLOR = configs::numOfColors-1;
}

// parse pdfstr into a vector<double>
vector<double> pdfstr_to_pdf(const string& pdfstr) {
  vector<string> tokens;

  {
    std::stringstream ss(pdfstr);
    string item;
    while(std::getline(ss, item, ',')) {
      tokens.push_back(item);
    }
  }

  vector<double> result;
  for(vector<string>::iterator it = tokens.begin(); it != tokens.end(); ++it) {
    std::istringstream i(*it);
    double val;
    if(!(i >> val)) {
      std::cerr << "ERROR: passed invalid pdf argument" << std::endl;
      exit(1);
    }
    
    result.push_back(val);
  }
  return result;
}

/** initialize gray_fraction
    
    @param g if g=0, then gray is not used as a color, and all colors
    are equally likely. If 0<g<1, then g is probability of a random
    point being gray, and all other colors are equally likely (i.e.,
    (1-g)/(numOfColors-1).
*/
void initialize_pdf(const double g, const string pdfstr = "") {
  vector<double> pdf;

  if(pdfstr == "") {
    configs::gray_fraction = g;

    if( configs::gray_fraction == configs::GRAY_UNUSED)
      return;

    // define: first color is gray
    // pdf[0] = gray_fraction
    // pdf[i] = (1-g)/(numOfCOlors-1) , when i!=0
    pdf.assign(configs::numOfColors, ((1-configs::gray_fraction) / (configs::numOfColors-1)) );
    pdf[0] = configs::gray_fraction;

  }
  else {
    TRACE(cout << "Initializing pdf from pdfstr=" << pdfstr << endl);
    configs::gray_fraction = configs::GRAY_UNUSED;
    pdf = pdfstr_to_pdf(pdfstr);
  }

  vector<double> cdf;
  cdf.reserve(pdf.size());
  double running_sum = 0;
  for(vector<double>::iterator it = pdf.begin(); it != pdf.end(); ++it) {
    running_sum += *it;
    cdf.push_back(running_sum);
  }
  configs::cdf = cdf;
  TRACE(cout << "Initialized configs::cdf to: " << configs::cdf << endl);
}

void initialize(const unsigned int alphabetsize, 
		const unsigned int length, 
		const unsigned int numOfColors,
		const double gray_fraction=0.0) {
  initialize_alphabet_size(alphabetsize);
  initialize_numOfColors(numOfColors);
  initialize_length(length);
  initialize_pdf(gray_fraction);
}

//   Returns mutants
vector<string> mut(const string g) { 
  TRACE(cout << " generating mutants of " << g << ":" << endl);
  vector<string> result;
  result.reserve((configs::alphabet_size -1) * configs::length );
  for(int pos = configs::length - 1; pos != -1; --pos) { 
    string::iterator alphabet_end;
    for(string::iterator alternative = configs::alphabet.begin(), 
	  alphabet_end = configs::alphabet.end();
	alternative != alphabet_end; ++alternative) {
      if (*alternative != g[pos]) {
	string mutant(g);
	mutant.reserve(configs::length); // may improve perf
	mutant[pos] = *alternative;
	result.push_back(mutant);
	TRACE(	cout << "  generated mutant " << mutant << endl);
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
  if( configs::gray_fraction == configs::GRAY_UNUSED )
    return rand() % configs::numOfColors;
  else {
    const double real = ((double)rand()) / ((double) RAND_MAX);
    TRACE(cout << "colorOf: generated real=" << real << endl);
    for(unsigned int i = 0; i < configs::numOfColors; ++i) {
      if( real < configs::cdf[i] )
	return i;
    }
  }
  // assert(ERROR)
  return configs::numOfColors-1;
}



// helper for set_difference
// (to be made function local with c++0x)
template <class T, class TVal>
struct is_contained_in {
  const unordered_map<T,TVal> & mm;
  is_contained_in(const unordered_map<T,TVal> & mmm) : mm(mmm) {}
  bool operator()(const T & item) { return (mm.find(item) != mm.end()); }  
};
/* removes any keys in m from s1.

   @param[inout] s1 a mathematical set of Ts
   @param[in] m an unordered_map with keys of T
*/
template <class T, class TVal>
vector<T> set_difference(vector<T> & s1, const unordered_map<T,TVal> & m) {
  s1.erase(std::remove_if(s1.begin(),s1.end(), is_contained_in<T,TVal>(m)),s1.end());
  return s1;
}


/* intersection of s1 and keys(m)
   
   @param[in] s1 a mathematical set of Ts
   @param[in] m an unordered_map with keys of T
*/
template <class T, class TVal>
vector<T> set_intersection(const vector<T> & s1, const unordered_map<T,TVal> & m) {
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
  observed[root]=configs::CLUSTER_COLOR;

  list<geno> to_traverse;
  to_traverse.push_back(root);

  geno cursor;
  while( !to_traverse.empty() ) {
    // process head of queue of nodes till to be traversed
    cursor = to_traverse.front();
    to_traverse.pop_front();
    TRACE(cout << "Traversing from node " << cursor << endl);
    // compute the neighbors not previously observed
    TRACE(cout << " Found " << cursor << endl);
    vector<geno> neighbors(mut(cursor));
    TRACE(cout << "  had neighbors: " << neighbors << endl); 
    set_difference(neighbors, observed);
    TRACE(cout << "  of which the previously unobserved are: " << neighbors << endl);
    // if there are some new nodes to observe ...
    if(!neighbors.empty()) {
      for(vector<geno>::iterator g = neighbors.begin(); 
	  g != neighbors.end(); ++g) {
	// "discover" the colors and record the visit
	if ((observed[*g] = colorOf(*g)) == configs::CLUSTER_COLOR) {
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
   network) and its mutational neighborhood, 'picked' from a random
   string graph with string LENGTH, built from an alphabet of
   ALPHABETSIZE, where every node is painted one color out of
   NUMOFCOLORS.
   
   'Picks' this cluster by doing a breadth-first search, assigning
   colors randomly as it progresses.
   
   NOTE:
   - modifies global variables alphabet and colors
   - alphabetsize and numOfColors must be < 26 

*/
unordered_map<geno,pheno> doRun(const unsigned int length, const unsigned int alphabetsize, const unsigned int numOfColors, const double gray=0.0) {
  initialize(alphabetsize,length,numOfColors,gray);
  const geno origin = string(configs::length, configs::alphabet[0]);
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

  TRACE(cout << "Processing map to calculate measures" << endl);

  // for every item ...
  for(unordered_map<geno,pheno>::const_iterator it = m.begin(); it != m.end(); ++it) {
    const geno g = it->first;
    const pheno p = it->second;
    TRACE(cout << "Searching g=" << g << endl);
    // .. in the perimeter ...
    if(p != configs::CLUSTER_COLOR) {
      TRACE(cout << "\twhich is in the perimeter" << endl);
      // ... tally it, track its pheno
      ++perimeter_size;
      perimeter_colors.insert(p);
      TRACE(cout << "\tadded its color" << endl);
      vector<geno> mutants( set_intersection( mut(g),m) );
      TRACE(cout << "\tcalculated its (observed) mutants: " << mutants << endl);
      for(vector<geno>::iterator it_mut = mutants.begin(); it_mut != mutants.end(); ++it_mut) {
	if(m.find(*it_mut)->second == configs::CLUSTER_COLOR) {
	  ++exits_size;
	}
      }
    }
  }

  results.perimeter_size = perimeter_size;
  results.cluster_size = m.size() - results.perimeter_size;
  results.colors = perimeter_colors.size();
  results.exits_size = exits_size;
  unsigned int mutations  = results.cluster_size * configs::length * (configs::alphabet_size - 1 );
  results.robustness = double(mutations - exits_size) / (double) mutations;
  return results;
}

extern "C"
cluster_measures calculate_measures(const unsigned int length, 
				    const unsigned int alphabetsize, 
				    const unsigned int numOfColors,
				    const double gray) {
  return calculate_measures_from_run(doRun(length,alphabetsize,numOfColors,gray));
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
					   const unsigned int samples,
					   const double gray) {
  using namespace boost::accumulators;

  accumulator_set<unsigned int, stats<tag::mean> > cluster_size_acc;
  accumulator_set<unsigned int, stats<tag::mean> > perimeter_size_acc;
  accumulator_set<unsigned int, stats<tag::mean> > colors_acc;
  accumulator_set<unsigned int, stats<tag::mean> > exits_size_acc;
  accumulator_set<double,       stats<tag::mean> > robustness_acc;
  
  cluster_measures r;
  for(unsigned int i = 0; i < samples; ++i ) {
    r = calculate_measures(length,alphabetsize,numOfColors,gray);
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


int main(int argc, char *argv[])
{
  namespace po = boost::program_options;

  unsigned int alphabetsize;
  unsigned int length;
  unsigned int numOfColors;
  double gray;
  string pdfstr;
  unsigned int samples;
  unsigned int seed;
  unsigned int verbosity;
  string mode;
  const unsigned int VERBOSITY_NONE = 0;
  const unsigned int VERBOSITY_LOW = 1;
  const unsigned int VERBOSITY_HIGH = 2;

  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help"										, "produce help message")
    ("alpha",  po::value<unsigned int>(&alphabetsize)	->default_value(2)		, "alphabet size")
    ("length", po::value<unsigned int>(&length)		->default_value(8)		, "string length")
    ("colors", po::value<unsigned int>(&numOfColors)	->default_value(3)		, "number of possible colors")
    ("gray",   po::value<double      >(&gray)		->default_value(0.0)		, 
     "probability a string is 'gray'\n"
     "Values:\n"
     "  gray=0: \tall colors have equal probability.\n"
     "  else  : \tone color has probability gray, and all other colors including the cluster color share an equal probability.")
    ("pdf", po::value<string>(&pdfstr)			->default_value("")	, "pdf")
    ("samples", po::value<unsigned int>(&samples)	->default_value(1)		, "number of searches to perform")
    ("seed",    po::value<unsigned int>(&seed)		->default_value(0)		, "initial pseudorandom seed (non-negative integer)")
    ("verbose", po::value<unsigned int>(&verbosity)	->default_value(1)		, 
     "verbosity\n"
     "  verbose=0: just results\n"
     "  verbose=1: results, inputs, labels\n"
     "  verbose=2: results, inputs, labels, observed nodes")
    ("mode", po::value<string>(&mode)			->default_value("stats")	, 
     "stats, data, or bench\n"
     "  stats: calculate means\n"
     "  data : dump raw search results\n"
     "  bench: silently perform 10000 searches")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    
  
  if (vm.count("help")) {
    cout << "Usage: clusters" << endl
	 << "Prints cluster measures from a search over a random string graph" << endl
	 << endl << desc << endl
	 << "Graph nodes are strings of length LENGTH, made from an alphabet of" << endl
	 << "ALPHA symbols, where each string is randomly assigned one of" << endl
	 << "COLORS possible colors." << endl
	 << endl
	 << "By default colors are dstributed uniformly; but if GRAY is non-zero," << endl
	 << "then GRAY defines the probability that a randomly chosen node is" << endl
	 << "GRAY, and the rest are distributed uniformly." << endl
	 << endl
	 << "The last color, color (COLORS-1), is the color of the cluster we searched." << endl
	 << "If GRAY is non-zero, the first color (color 0) is gray." << endl
      ;
    return 1;
  }


  //  srand(time(NULL)); // seed the random number generator
  std::srand(seed); // seed the random number generator

  if(verbosity > VERBOSITY_NONE) {
    cout << "searching with:" << endl;
    cout << "\talphabetsize = " << alphabetsize << endl;
    cout << "\tlength = " << length << endl;
    cout << "\tnumOfColors = " << numOfColors << endl;
    cout << "\tgray = " << gray << endl;
    cout << "\tseed = " << seed << endl;
    cout << "\tmode = " << mode << endl;
    cout << "\tpdf  = " << pdfstr << endl;
  }

  // display one search
  if (mode == "data") {
    if( verbosity > VERBOSITY_NONE) 
      cout << endl << "As mode=data, dumping results from " << samples << " searches" << endl;
      
    for(unsigned int i = 0; i < samples; ++i) {
      unordered_map<geno,pheno> mm(doRun(length,alphabetsize,numOfColors,gray));

      if(verbosity == VERBOSITY_HIGH) {
	cout << "Where cluster has color 0, and gray (if defined) has the maximum color, observed nodes as follows: " << endl;
	cout << mm << endl;
      }    

      cluster_measures results = calculate_measures_from_run(mm);
      if( verbosity > VERBOSITY_NONE ) {
	cout << endl;
	cout << "results.cluster_size   = s = " << results.cluster_size << endl;
	cout << "results.perimeter_size = t = " << results.perimeter_size << endl;
	cout << "results.colors_seen    = E = " << results.colors << endl;
	cout << "results.exits_size     = u = " << results.exits_size << endl;
	cout << "results.robustness     = r = " << results.robustness << endl;
      } 
      else if( verbosity  == VERBOSITY_NONE ) {
	cout << results.cluster_size	<< "\t";;
	cout << results.perimeter_size	<< "\t";;
	cout << results.colors		<< "\t";;
	cout << results.exits_size	<< "\t";;
	cout << results.robustness	<< endl;
      }
    }
  }
  else if (mode=="stats") {
    if(verbosity > VERBOSITY_NONE)
      cout << endl << "Mode=stats. Calculating statistics over " << samples << " searches." << endl;

    mean_cluster_measures results = calculate_statistics(length,alphabetsize,numOfColors,samples,gray);
    
    if( verbosity > VERBOSITY_NONE ) {
      cout << "mean results.cluster_size      = s = " << results.mean_cluster_size << endl;
      cout << "mean results.perimeter_size    = t = " << results.mean_perimeter_size << endl;
      cout << "mean results.colors_seen       = E = " << results.mean_colors << endl;
      cout << "mean results.exits_size        = u = " << results.mean_exits_size << endl;
      cout << "mean results.robustness        = r = " << results.mean_robustness << endl;
    }
    else if( verbosity  == VERBOSITY_NONE ) {
      cout << results.mean_cluster_size << endl;
      cout << results.mean_perimeter_size << endl;
      cout << results.mean_colors << endl;
      cout << results.mean_exits_size << endl;
      cout << results.mean_robustness << endl;
    }
  }
  else if (mode=="bench") {
    // benchmark 1000 random searches
    for(int i =0; i < 10000; ++i) {
      doRun(length,alphabetsize,numOfColors);
    }
  }
  else if(mode=="bench2") {
    // check 2nd call of random
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
