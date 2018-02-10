# EigenIterators
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
