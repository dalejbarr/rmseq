#include "rmseq.hpp"

//#include <iostream>
//using std::cout;
//using std::endl;

using std::string;
using std::vector;

vector<string> render_string(const vector<size_t> &vix,
			     const vector<string> &labels) {
  // render a vector<string> from a vector<size_t> of indices
  vector<string> vs(vix.size());
  for (size_t j = 0; j < vix.size(); j++) {
    vs[j] = labels[vix[j]];
  }

  return vs;
}

vector<vector<string>> render_string_set(const vector<vector<size_t>> &vvix,
					 const vector<string> &labels) {
  // render a vector<vector<string>> from vector<vector<size_t> of indices
  vector<vector<string>> vvs(labels.size());

  for (size_t i = 0; i < labels.size(); i++) {
    vvs[i] = render_string(vvix[i], labels);
  }

  return vvs;
}

// create a single complement set from a nominal variable
vector<vector<string>> com_n(const vector<string> &labels,
			     const int &n_buckets,
			     const int &n_reps,
			     const bool &edge_streaks,
			     const bool &repeat,
			     Permuter * pperm /*= nullptr*/) {

  vector<vector<size_t>> vvix = com_n_ix(labels, n_buckets, n_reps,
					 edge_streaks, repeat, pperm);
  
  vector<vector<string>> vvs = render_string_set(vvix, labels);
  
  return vvs;
}

// create a single complement set from a nominal variable
vector<vector<size_t>> com_n_ix(const vector<string> &labels,
				const int &n_buckets,
				const int &n_reps,
				const bool &edge_streaks,
				const bool &repeat,
				Permuter * pperm /*= nullptr*/) {

  vector<size_t> vs_ix = seq_n_ix(labels, n_buckets, n_reps,
				  edge_streaks, repeat, pperm);
  
  vector<vector<size_t>> vvs_ix = complement_set(vs_ix, labels.size());
    
  return vvs_ix;
}

// single sequence
vector<string> seq_n(const vector<string> &labels,
		     const int &n_buckets,
		     const int &n_reps,
		     const bool &edge_streaks,
		     const bool &repeat,
		     Permuter * pperm /*= nullptr*/) {
  
  vector<size_t> seq = seq_n_ix(labels, n_buckets, n_reps,
				edge_streaks, repeat, pperm);

  return render_string(seq, labels);  
}

// single sequence
vector<size_t> seq_n_ix(const vector<string> &labels,
			const int &n_buckets,
			const int &n_reps,
			const bool &edge_streaks,
			const bool &repeat,
			Permuter * pperm /*= nullptr*/) {

  if (pperm == nullptr) {
    // allocate randomizer
  }

  vector<size_t> vs_ix(labels.size() * n_buckets * n_reps);  

  // basis vector for a single bucket
  vector<size_t> basis(labels.size() * n_reps);
  for (size_t i = 0; i < labels.size(); i++) {
    for (size_t j = 0; j < (size_t) n_reps; j++) {
      basis[i * n_reps + j] = i;
    }
  }

  if (repeat) {
    vector<size_t> scramble = pperm->permute(basis.size());
    for (size_t i = 0; i < (size_t) n_buckets; i++) {
      for (size_t j = 0; j < scramble.size(); j++) {
	vs_ix[i * scramble.size() + j] = basis[scramble[j]];
      }
    }
  } else {
    if (edge_streaks) { // "free" method
      for (size_t i = 0; i < (size_t) n_buckets; i++) {
	vector<size_t> scramble = pperm->permute(basis.size());
	for (size_t j = 0; j < scramble.size(); j++) {
	  vs_ix[i * scramble.size() + j] = basis[scramble[j]];
	}
      }
    } else { // avoid streaks at edges
      for (size_t i = 0; i < (size_t) n_buckets; i++) {
	vector<size_t> scramble = pperm->permute(basis.size());
	if (i > 0) { // correct for overlap with end of previous bucket
	  size_t last_value = vs_ix[i * scramble.size() - 1];
	  if (scramble[0] == last_value) { // first matches last; correct
	    // get an index to be swapped (1 to n-1)
	    size_t swaptarg_ix = pperm->permute(basis.size() - 1)[0] + 1;
	    // swap the values
	    size_t tmp = scramble[0];
	    scramble[0] = scramble[swaptarg_ix];
	    scramble[swaptarg_ix] = tmp;
	  }
	} // if (i > 0)
	for (size_t j = 0; j < scramble.size(); j++) {
	  vs_ix[i * scramble.size() + j] = basis[scramble[j]];
	}
      }
    }
  }
  
  return vs_ix;
}

vector<vector<size_t>> complement_set(const vector<size_t> &src,
					  const size_t &n_levels) {
    // create complement vectors
  vector<vector<size_t>> vvs_ix(n_levels);

  for (size_t i = 0; i < n_levels; i++) {
    if (i == 0) {
      vvs_ix[0] = src;
    } else {
      vvs_ix[i] = vector<size_t>(src.size());
    }
  }
  
  for (size_t i = 1; i < n_levels; i++) {
    for (size_t j = 0; j < src.size(); j++) {
      vvs_ix[i][j] = (vvs_ix[i - 1][j] + 1) % n_levels;
    }
  }

  return vvs_ix;  
}
