## Dynamic de Bruijn Graph implementation

Copyright Alan Kuhnle, Victoria Crawford 

Compact data structures for the representation of de Bruijn graphs are important for genome sequence assembly. A fully dynamic data structure supporting efficient and exact membership queries for de Bruijn graphs was recently proposed in [1]. For our project, we have partially implemented this data structure and evaluated its practical performance as compared to popular implementations based upon bloom filters when used for de Bruijn construction from sequence data. For a list of features that remain unimplemented, see the implementation to-do list below.

A detailed description of our implementation is contained in the paper.

[1]: Belazzougui, D., Gagie, T., Maekinen, V., and Previtali, M. (2016). Fully Dynamic de Bruijn Graphs. In String Processing and Information Retrieval, pages 145–152.

### Dependencies
The fully dynamic De Bruijn data structure is implemented in class FDBG in cpp-src/FDBG.cpp.
This implementation depends only on

Boost: http://www.boost.org/doc/libs/1_59_0/more/getting_started/unix-variants.html

Installing Boost on a Debian-based system should require only the following command:
   ```
   sudo apt-get install libboost-all-dev
   ```

An example program using the FDBG class is provided in cpp-src/BuildDataStructure.cpp. A Makefile to compile this file is provided in cpp-src, which may need to be edited to include the Boost location on your system.

Comparison of FDBG to various bloom filters is provided in file main.cpp, which has the additional
dependency

libbf: https://github.com/mavam/libbf/ 

For the minimal perfect hash function, we use the BBhash implementation; this consists of a single header file included in our repository. The implementation is available here: https://github.com/rizkg/BBHash/

### Provided example dataset
In the cpp-src directory, a file yeast.fasta is present. This is the reference genome for yeast strain S288C in fasta format; it can be used to test our data structure with the BuildDataStructure.cpp example program.

### Instruction to run BuildDataStructure example program
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


### Instructions to run main (comparison with bloom filters)
- Go to ./cpp-src
- Compile
    ```
    sh compile.sh main [random_mode = "shift"] [# of queries = 1000000] 
    ```
- Test
    ```
    ./main data/yeast.fasta 10 | tee log/terminal_log.txt
    ```

log.txt is written in log/.

terminal_log.txt save all system output for python parsing.

### Demo
Written in flask and bootstrap. Go to web and run 
```
python main.py
```
Open a web browser with an ip the program provided.

### Implementation to-do
The following is a list of items
we were unable to get to during this project.
1. Implement dynamic perfect hash function
2. Update injectivity test in hash construction to use the O(1) update of hash function
3. Use O(1) update to construct IN,OUT
4. Update $k$-mer representation to allow $k > 32$
5. Change underlying BitArray type to fixed-width integer type (\texttt{uint8\_t} would be best)
6. The forest construction, if run on a small connected component, does not guarantee it has sampled the best root. That is, a component may have tree below minimum height when a better root choice could ensure the minimum height. This is a bug that needs to be fixed.
The merge_tree function may not merge trees in all possible cases, as discussed in our report. Also, dynamic adding and deleting functions do not always take advantage of the Karp-Rabin update whenever possible.
