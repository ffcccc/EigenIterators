# EigenIterators
Defines an iterator pointing to the first/last element in a Eigen::Array object

Since the Eigen::Array objects have no member functions begin/end defined nor have iterator member types, this snippet of code implements a random-access iterator of an unspecified type whose value_type is T and its reference type is T& for the mutable version and const T& for the constant iterator declaration.

This is an adaptation of begin/end helper functions defined for std::valarray, which in turn is an overload of std::begin/end, defined in <iterator>.

Parameters:
x -> A Eigen::Array object.

