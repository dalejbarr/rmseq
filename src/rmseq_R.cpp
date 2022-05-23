#include <Rcpp.h>
#include "rmseq.hpp"
#include "Permuter_R.hpp"

//#include <iostream>
//using std::cout;
//using std::endl;

using namespace Rcpp;
using std::vector;
using std::string;

// [[Rcpp::export(.seq_n)]]
CharacterVector seq_n(CharacterVector labels,
		      int n_buckets,
		      int n_reps,
		      bool edge_streaks,
		      bool repeats) {

  Permuter_R perm;

  vector<string> vvs = seq_n(as<vector<string>>(labels),
			     n_buckets,
			     n_reps,
			     edge_streaks,
			     repeats,
			     &perm);
  
  return wrap(vvs);
}

// [[Rcpp::export(.com_n)]]
List com_n(CharacterVector labels,
	   int n_buckets,
	   int n_reps,
	   bool edge_streaks,
	   bool repeats) {

  Permuter_R perm;

  vector<vector<string>> vvs = com_n(as<vector<string>>(labels),
				     n_buckets,
				     n_reps,
				     edge_streaks,
				     repeats,
				     &perm);

  return wrap(vvs);
}

// [[Rcpp::export(.seq_n_ix)]]
IntegerVector seq_n_ix_R(CharacterVector labels,
			 int n_buckets,
			 int n_reps,
			 bool edge_streaks = false,
			 bool repeats = false) {

  Permuter_R perm;
  
  vector<size_t> vv_ix = seq_n_ix(as<vector<string>>(labels),
				  n_buckets, n_reps,
				  edge_streaks, repeats, &perm);

  IntegerVector result(vv_ix.size());
  for (size_t i = 0; i < vv_ix.size(); i++) {
    result[i] = static_cast<int>(vv_ix[i]) + 1;
  }

  return wrap(result);
}

// [[Rcpp::export(.com_n_ix)]]
List com_n_ix_R(CharacterVector labels,
		int n_buckets,
		int n_reps,
		bool edge_streaks = false,
		bool repeats = false) {

  Permuter_R perm;
  
  vector<vector<size_t>> vv_ix = com_n_ix(as<vector<string>>(labels),
					  n_buckets, n_reps,
					  edge_streaks, repeats, &perm);

  List result(vv_ix.size());
  for (size_t i = 0; i < vv_ix.size(); i++) {
    IntegerVector iv(vv_ix[i].size());
    for (size_t j = 0; j < vv_ix[i].size(); j++) {
      iv[j] = static_cast<int>(vv_ix[i][j]) + 1;
    }
    result[i] = iv;
  }

  return wrap(result);
}

