#ifndef GUARD_rmseq_hpp
#define GUARD_rmseq_hpp

#include <vector>
#include <string>
#include "Permuter.hpp"

std::vector<std::string>
  render_string(const std::vector<std::size_t> &vix,
		const std::vector<std::string> &labels);

std::vector<std::vector<std::string>>
  render_string_set(const std::vector<std::vector<std::size_t>> &vvix,
		    const std::vector<std::string> &labels);

std::vector<std::vector<std::string>>
  com_n(const std::vector<std::string> &labels,
	const int &n_buckets,
	const int &n_reps,
	const bool &edge_streaks,
	const bool &repeat,
	Permuter * pperm = nullptr);

std::vector<std::vector<std::size_t>>
  com_n_ix(const std::vector<std::string> &labels,
	   const int &n_buckets,
	   const int &n_reps,
	   const bool &edge_streaks,
	   const bool &repeat,
	   Permuter * pperm = nullptr);

std::vector<std::string> seq_n(const std::vector<std::string> &labels,
			       const int &n_buckets,
			       const int &n_reps,
			       const bool &edge_streaks,
			       const bool &repeat,
			       Permuter * pperm = nullptr);

std::vector<std::size_t> seq_n_ix(const std::vector<std::string> &labels,
				  const int &n_buckets,
				  const int &n_reps,
				  const bool &edge_streaks,
				  const bool &repeat,
				  Permuter * pperm = nullptr);

std::vector<std::vector<std::size_t>>
  complement_set(const std::vector<std::size_t> &src,
		 const std::size_t &n_levels);

#endif
