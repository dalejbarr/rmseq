#ifndef GUARD_Permuter_hpp
#define GUARD_Permuter_hpp

#include <vector>

// abstract class for randomly permuting indices
class Permuter {
public:
  Permuter() {}
  virtual ~Permuter();
  virtual std::vector<std::size_t> permute(std::size_t n) = 0;
  virtual std::vector<int> permute(int n) = 0;
};

#endif
