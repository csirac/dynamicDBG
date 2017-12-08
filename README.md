# Dynamic de Bruijn Graph implementation

Copyright Alan Kuhnle, Victoria Crawford

Compact data structures for the representation of de Bruijn graphs are important for genome sequence assembly. A fully dynamic data structure supporting efficient and exact membership queries for de Bruijn graphs was recently proposed in [1].
This software is an implementation of the data structure in [1] with dynamic
edges (and not dynamic nodes).

[1]: Belazzougui, D., Gagie, T., Maekinen, V., and Previtali, M. (2016). Fully Dynamic de Bruijn Graphs. In String Processing and Information Retrieval, pages 145–152.

## Dependencies

This implementation depends on the Boost C++ library:
http://www.boost.org/doc/libs/1_59_0/more/getting_started/unix-variants.html

Installing Boost on a Debian-based system should require only the following command:
   ```
   sudo apt-get install libboost-all-dev
   ```
We use the BBhash implementation of a minimal perfect hash function, which is
included in our source code.
The implementation is available here: https://github.com/rizkg/BBHash/

## Dynamic de Bruijn Graph
The dynamic De Bruijn graph is implemented in the class FDBG in cpp-src/FDBG.cpp.

We provide a program cpp-src/BuildDataStructure.cpp that takes in
a fasta file and a k-mer length. A Makefile to compile this file is
provided in cpp-src, which may need to be edited to include the Boost
location on your system.

### Example dataset
In the cpp-src directory, a file yeast.fasta is present. This is the reference genome for yeast strain S288C in fasta format; it can be used to test our data structure with the BuildDataStructure.cpp example program.

### Instructions to run BuildDataStructure.cpp
- Install the Boost dependency (see dependencies above).
- Go to ./cpp-src
- Compile using the included Makefile, editing it if necessary.
    ```
    make
    ```
- After compilation, the executable bds should have been produced. This program requires two arguments.
The first argument is any fasta file, for example the included yeast.fasta file. The second argument
is the desired value of k for de Bruijn graph, must be between 1 and 31.
Example:
    ```
    ./bds yeast.fasta 20
    ```

### Description of output for BuildDataStructure
During execution, bds will print detailed logs to the screen. These logs indicate when each step of the construction of the data structure begins. In order, these steps are the following
- Parses fasta file
- Writes binary files for the k-mers found so that on future runs, these files are loaded instead
- Begins construction of De Bruijn graph
  + Construction of Karp-Rabin hash function
  + Injectivity test
  + Construction of minimal perfect hash function
  + Computation of r inverse modulo P
  + Construction of IN, OUT matrices
  + Forest construction
- Once construction is complete, the data structure is saved in binary format to a file.
- Tree height tests are conducted on the forest
- Edges are dynamically removed from the data structure, and then added back in;
  the forest is updated after each change.
