#include "solution.hpp"

static int getSumOfDigits(int n) {
  int sum = 0;
  while (n != 0) {
    sum = sum + n % 10;
    n = n / 10;
  }
  return sum;
}

int solution(const hash_map_t *hash_map, const std::vector<int> &lookups) {
  int result = 0;
  constexpr int LOOKAHEAD = 16;

  {
    int i = 0;
    for (; i < lookups.size() - LOOKAHEAD; ++i) {
      const int val = lookups[i];
      if (hash_map->find(val))
        result += getSumOfDigits(val);
      
      #ifdef SOLUTION
      hash_map->prefetch(lookups[i + LOOKAHEAD]);
      #endif
    }

    for (; i < lookups.size(); ++i) {
      const int val = lookups[i];
      if (hash_map->find(val))
        result += getSumOfDigits(val);
    }
  }

  return result;
}
