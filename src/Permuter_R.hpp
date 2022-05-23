#ifndef Permuter_R_hpp
#define Permuter_R_hpp

#include "Permuter.hpp"

class Permuter_R : public Permuter {
public:
  Permuter_R() {}
  virtual ~Permuter_R();
  virtual std::vector<std::size_t> permute(std::size_t n);
  virtual std::vector<int> permute(int n);
};

#endif
