#include <cstring>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>
#include "boost/tr1/unordered_map.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"

#include "printable.hpp"

using std::tr1::unordered_map;
using std::string;
using std::vector;

typedef std::string geno;
typedef unsigned int pheno;

const  pheno cluster_color = 0;

// true if Hamming distance=1
inline
bool unit_distance(const char *  xstr,const char  *  ystr) {
    size_t distance = 0;
    size_t len = std::strlen(xstr);
    for(size_t pos =0; pos != len; ++pos) {
        if ( xstr[pos] != ystr[pos] ) {
            ++distance;
            if( distance > 1 ) {
                return false;
            }
        }
    }

    if (distance == 0)
        return false;
    return true;
}

bool unit_distance(const std::string & x,const std::string & y) {
    return unit_distance(x.c_str(),y.c_str());
}

typedef std::pair<geno,geno> Edge;

vector< Edge > 
get_edges(const unordered_map<geno,pheno> & m) {
    vector< Edge > edges;

    // brute force O(n^2) enumeration of non-self, undirected edges
    const unordered_map<geno,pheno>::const_iterator mend = m.end();
    for(unordered_map<geno,pheno>::const_iterator itx = m.begin(); itx != mend; ++itx) {
        unordered_map<geno,pheno>::const_iterator ity = itx;
        std::advance(ity,1);
        for(; ity != mend; ++ity) {
            if ( unit_distance(itx->first, ity->first ) ) {
                edges.push_back( Edge(itx->first, ity->first ) ); 
            }
        }
    }
    return edges;
}

std::string graphToDot(const unordered_map<geno,pheno> & m) {
  std::ostringstream result;

  result << "graph mygraph {" "\n";

  vector< Edge > edges( get_edges( m ) );
  for(vector<Edge>::const_iterator it = edges.begin(); it != edges.end(); ++it) {
    result << " " << it->first << " -- " << it->second << ";" "\n";
  }
  result << "}" "\n";
  return result.str();
}

int main() {
    unordered_map<geno,pheno> m;
    m["AA"]=0;
    m["AB"]=0;
    m["BA"]=1;
    m["BB"]=1;

    std::cout << graphToDot(m) << std::endl;
}
