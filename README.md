# CBNE C++ Implementation
## Building

### Boost libary

Unfortunately, the C++ implementation makes use of the Boost library which is a bit of a beast. Building is thus not as straightforward as could be hoped for.

You can follow the instructions in the [Boost documentation](https://www.boost.org/doc/libs/1_68_0/more/getting_started/unix-variants.html) to build and install Boost libraries. We repeat them here for clarity. After having downloaded and unzipped the library navigate, in a terminal window, into the newly unzipped folder. Then run

```
./bootstrap.sh
```

```
./b2 isntall
```

Note the location the library has been installed to. Then change the following line in [CMakeLists.txt](https://github.com/CQCL-DEV/CTDA-Algs/blob/main/CBNE-Serial-C%2B%2B/CMakeLists.txt) so that the path points to the location you have installed the Boost libraries to.

```
set(BOOST_ROOT "/home/USER/boost_libraries/")
```

### CBNE

Now you can run `cmake` and build as per any CMake project. This can be done either from the command line or using the CMake extension in Visual Code. If building from the command line, navigate to `CTDA-Algs/Apres-C++/` then do 

```
mkdir build
cd build
cmake ..
make
```
If everything has gone right you should now have a working executable!

## Running

To run the code:

```
USAGE: 

   apres  [-h] [--version] [-e <double>] [-p <string>]


Where: 

   -p <string>,  --path <string>
     Path to .graphml input file

   -e <double>,  --epsilon <double>
     Epsilpn value

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

Please note that the implementation is under active development, so we should be seeing an increase in performance and further options.

## Running Unit Tests

If you want to build and run the unit tests, you will first need to install Catch2 v3 test framework. To do so, you can use the following commands.

```
$ git clone https://github.com/catchorg/Catch2.git
$ cd Catch2
$ cmake -Bbuild -H. -DBUILD_TESTING=OFF
$ sudo cmake --build build/ --target install
```

Then, in [CMakeLists.txt](https://github.com/CQCL-DEV/QTDA-Algs/blob/ahmed-apres/ctda/Apres-C%2B%2B/CMakeLists.txt) change `if(FALSE)` to `if(TRUE)` and the build as normal. The tests can be run using `ctest` or by running them wuth Catch2.

# Docker

You can also build and run using the supplied dockerfile. In this case, you can ignore all the instructions above. The instructions below assume that you have Docker installed on your machine.

In the root of this repo, run

```
sudo docker build -t <your_image_name:tag> -f Dockerfile .
```

This builds the image from the dockerfile. To actually run the application use:

```
sudo docker run -v <path to .graphml files on your machine>:/benchmarks --rm <your_image_name:tag> -p /benchmarks/<path to benchmark> <other options>
```
