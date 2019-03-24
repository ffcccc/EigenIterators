# Eigen::Utils

**EigenUtis** represents the effort of organizing a growing set of snippets into a coherent collection of algorithms that support the amazing ([Eigen Matrix](http://eigen.tuxfamily.org/index.php)) library.

**EigenUtils** makes it easier to create your computational algorithms, by eliminating a lot of the setup burden and integrating Eigen code with the rest of the code:

- You don't need to install/configure
- You don't need to install runtime dependencies
- You don't need no more to copy your tensors from the Eigen representatio to another (e.g. STL) to apply standard library's operators
- It's easy to try out, you can just delete your forked repository if you don't like it

In a few minutes you'll be set up with a minimal, responsive header-only set of libs giving you more time to spend on writing code.

![Jekyll Now Theme Screenshot](/images/jekyll-now-theme-screenshot.jpg "Jekyll Now Theme Screenshot")


## EigenUtils Features

Among the subjects covered by the lib you can find:

## iter -> STL Iterators for Eigen containers 
Defines an iterator pointing to the first/last element in a Eigen::Array object

Since the Eigen::Array objects have no member functions begin/end defined nor have iterator member types, this snippet  implements a random-access iterator of an unspecified type whose value_type is T and its reference type is T& for the mutable version and const T& for the constant iterator declaration.

This is an adaptation of begin/end helper functions defined for std::valarray, which in turn is an overload of std::begin/end, defined in <iterator>.

### Parameters:
x -> A Eigen::Array object.

### Usage
```c++
...
#define NA 1000000.0

  ArrayXd side(side_par);
  side << 2.0, 2.0, NA, NA, 2.0;
  
  replace_if(begin(side), end(side), [](double i) { return (i == NA); }, 1.0);
  if( any_of(begin(side), end(side), [](double i) { return (fabs(i) != 1.0); }) ) {
    return 1;
  }
  return 0;
...
 
```

## stats -> median, mean, variance
This section started facing the problem of optimize computation of median in a vectorized environment. An interesting hint came from the [stackoverflow community](https://stackoverflow.com/questions/1719070/what-is-the-right-approach-when-using-stl-container-for-median-calculation/1719155#1719155):  
         
>This algorithm handles both even and odd sized inputs efficiently using the STL nth_element (amortized O(N)) algorithm
and the max_element algorithm (O(n)). Note that nth_element has another guaranteed side effect, namely that all of the elements
before n are all guaranteed to be less than v[n], just not necessarily sorted.

### Usage
```c++
...
 
```

## lm
Fast linear regression solver. Thanks to Rcpp/REigen (link)[].
### Usage
```c++
...
 
```

## lr
Logistic regression solver, by means of gradient descent of the error function. IRLS to come...
### Usage
```c++
...
 
```

## dist
Distance measures
### Usage
```c++
...
 
```

## corr
Correlation measures.
### Usage
```c++
...
 
```

## norm
Normal distribution calculations. Mutauated from Plink and R
### Usage
```c++
...
 
```

## numDeriv
First and second order numerical derivatives
### Usage
```c++
...
 
```

## assoc
Statistical association tests, contingency tables, binary similarity measures
### Usage
```c++
...
 
```


## Quick Start
I've created a main.cpp file in the root directory, containing a walkthrough and test application. It's a good starting point to check for use-cases and It's going to get more detailed and curated in the coverage of tests.

### Step 1) Clone EigenUtils to your working path

### Step 2) Customize the 'include path' of your project

![_config.yml](/images/config.png "_config.yml")

### Step 3) Use It

![First Post](/images/first-post.png "First Post")

## Requirements

1. Of course you will need a local copy of the Eigen repository
2. A restricted number of algoritms (e.g. the T-test code) still depend on the Boost library for some math


## Credits

- [Jekyll](https://github.com/jekyll/jekyll) - Thanks to its creators, contributors and maintainers.

## Contributing

Issues and Pull Requests are greatly appreciated. If you've never contributed to an open source project before I'm more than happy to walk you through how to create a pull request.

You can start by [opening an issue](https://github.com/barryclark/jekyll-now/issues/new) describing the problem that you're looking to resolve and we'll go from there.

I want to keep Jekyll Now as minimal as possible. Every line of code should be one that's useful to 90% of the people using it. Please bear that in mind when submitting feature requests. If it's not something that most people will use, it probably won't get merged. :guardsman:

