# Eigen Iterators
Defines an iterator pointing to the first/last element in a Eigen::Array object

Since the Eigen::Array objects have no member functions begin/end defined nor have iterator member types, this snippet of code implements a random-access iterator of an unspecified type whose value_type is T and its reference type is T& for the mutable version and const T& for the constant iterator declaration.

This is an adaptation of begin/end helper functions defined for std::valarray, which in turn is an overload of std::begin/end, defined in <iterator>.

## Parameters:
x -> A Eigen::Array object.

## Demo
```c++
#include <algorithm>
#include "Eigen/Core"

#define NA 1000000.0

using namespace std;
using namespace Eigen;

int main() {
  ArrayXd side(side_par);
  side << 2.0, 2.0, NA, NA, 2.0;
  
  replace_if(begin(side), end(side), [](double i) { return (i == NA); }, 1.0);
  
  if( any_of(begin(side), end(side), [](double i) { return (fabs(i) != 1.0); }) ) {
    return 1;
  }
  return 0;
 }
 ```

# Eigen median, mean, variance
post-condition: After returning, the elements in v may be reordered and the resulting order is implementation defined.

https://stackoverflow.com/questions/1719070/what-is-the-right-approach-when-using-stl-container-for-median-calculation/1719155#1719155         
This algorithm handles both even and odd sized inputs efficiently using the STL nth_element (amortized O(N)) algorithm
and the max_element algorithm (O(n)). Note that nth_element has another guaranteed side effect, namely that all of the elements
before n are all guaranteed to be less than v[n], just not necessarily sorted.


# Eigen fast linear regression solver

# Eigen distance

# Eigen contingency tables