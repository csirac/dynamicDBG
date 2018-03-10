/**

Copyright 2017 Alan Kuhnle, Victoria Crawford

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/




#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#define BOOST_LOG_DYN_LINK 1
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <chrono>
#include "FDBG.cpp"
#include "TestUtil.cpp"
#include "formatutil.cpp"

using namespace std;


/**
 * Run like "./a.out kmers.fasta 8"
 */
int main(int argc, char* argv[]) {
  
   // fasta filename
   string filename = argv[1];
   int k = stoi(argv[2]);
   FDBG Graph;

   //Have we constructed this dataset before on these parameters?
   string dsfile = filename.substr( 0, filename.find_last_of( '.' ) ) + "fdbg" + to_string( k ) + ".bin";
   // get k-mers and edgemers from file
   unordered_set<kmer_t> kmers;
   unordered_set<kmer_t> edgemers;
   handle_mers( filename, k, kmers, edgemers );

   if (file_exists( dsfile )) {
     //Yes, so avoid reconstructing
     BOOST_LOG_TRIVIAL(info) << "Loading from " << dsfile;
     ifstream ifile_ds( dsfile.c_str(), ios::in | ios::binary );
     Graph.load( ifile_ds );
     ifile_ds.close();
   } else {
     BOOST_LOG_TRIVIAL(info) << "Building De Bruijn Graph ...";
     Graph.build( kmers, edgemers, kmers.size(), k);
     ofstream ofile( dsfile.c_str(), ios::out | ios::binary );
     Graph.save( ofile );
     ofile.close();
   }

	// CGAC
	kmer_t new_node_1 = 97;
	Graph.addNode(new_node_1);
	// GACG
	kmer_t new_node_2 = 134;
	Graph.addNode(new_node_2);
	// CGAC
	kmer_t new_node_3 = 26;
	Graph.addNode(new_node_3);

 	Graph.printHashFunction(kmers);
	Graph.printAddedHashFunction();

	Graph.printINandOUT(kmers);		
}

