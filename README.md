# Template CXX
Template CXX project configuration with CMake for a cross-platform application.

## Dependencies
Relies on the [Boost C++](https://www.boost.org/) library.

To build using cmake, it needs the following entries set:
- BOOST_DIR CACHE PATH "Boost installation directory path".
  Example:
  ```
  -DBOOST_DIR="/path/to/boost/installation-folder/"
  ```
- BOOST_PRG_OPTS_LIB CACHE FILEPATH "Boost::program_options static lib filepath".
  Example:
  ```
  -DBOOST_PRG_OPTS_LIB="/path/to/libboost_program_options.a"
  ```

The particular Boost library used in this template is [Program Options](https://www.boost.org/doc/libs/1_84_0/doc/html/program_options.html).
