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
#include "formatutil.cpp"
#include "TestUtil.cpp"

using namespace std;


/**
 * Run like "./a.out kmers.fasta 8"
 */
int main(int argc, char* argv[]) {
   // Set debug level
   boost::log::core::get()->set_filter(boost::log::trivial::severity
      >= boost::log::trivial::info);

   BOOST_LOG_TRIVIAL(info) << "Beginning to build data structure ...";

   // Check if the user put in the correct command line arguments
   if (argc < 3) {
         BOOST_LOG_TRIVIAL(fatal) << "Missing required arguments. Usage:"
				  << argv[0] << " <reads.fasta> <k>";
         exit(1);
   }
  
   // fasta filename
   string filename = argv[1];
   int k = stoi(argv[2]);
   FDBG Graph;

   //Have we constructed this dataset before on these parameters?
   string dsfile = filename.substr( 0, filename.find_last_of( '.' ) ) + "fdbg" + to_string( k ) + ".bin";
   // get k-mers and edgemers from file
   unordered_set<kmer_t> kmers;
   unordered_set<kmer_t> edgemers;
   auto start = std::chrono::system_clock::now();
   handle_mers( filename, k, kmers, edgemers );
   auto end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end-start;
   BOOST_LOG_TRIVIAL(info) << "Getting " << kmers.size() + edgemers.size() << " mers took " << elapsed_seconds.count() << " s";

   if (file_exists( dsfile )) {
     //Yes, so avoid reconstructing
     BOOST_LOG_TRIVIAL(info) << "Loading from " << dsfile;
     ifstream ifile_ds( dsfile.c_str(), ios::in | ios::binary );
     Graph.load( ifile_ds );
     ifile_ds.close();
     BOOST_LOG_TRIVIAL(info) << "Data structure built in " << Graph.construction_time << " s";
   } else {
     cout << "PAY ATTENTION!! IT IS SLOW!!" << endl;
     BOOST_LOG_TRIVIAL(info) << "Building De Bruijn Graph ...";
     Graph.build( kmers, edgemers, kmers.size(), k);

     BOOST_LOG_TRIVIAL(info) << "Data structure built in " << Graph.construction_time << " s";
     BOOST_LOG_TRIVIAL(info) << "Writing data structure to file " << dsfile;
     ofstream ofile( dsfile.c_str(), ios::out | ios::binary );
     Graph.save( ofile );
     ofile.close();
   }

   //   BOOST_LOG_TRIVIAL(info) << "Estimated size(bits) (Mb):" << Graph.estimateBitSize()/(8.0 * 1024 * 1024);
   //   BOOST_LOG_TRIVIAL(info) << "Size(Mb):" << Graph.bitSize() / (8.0 * 1024 * 1024);
   //   BOOST_LOG_TRIVIAL(info) << "Bits per element:" << Graph.bitSize() / static_cast<double>( Graph.n );


   /**
    * First batch of tree height tests
    */
   BOOST_LOG_TRIVIAL(info) << "Tree height tests...";
   
   // Compute data about trees in the forest
   //BOOST_LOG_TRIVIAL(info) << "Fixing trees...";
   //Graph.repairTrees();
   vector< size_t > treesBelowSize;
   vector< size_t > nCCs;
   
   unsigned num_trees;
   double avg_height; // average height of trees
   unsigned num_above; // number above the supposed max
   unsigned num_below; // number below the supposed min
   map< unsigned, unsigned > height_dist;
   unsigned minHeight;
   unsigned maxHeight;
   Graph.getTreeData(num_trees, avg_height, num_above, num_below ); //, height_dist, minHeight, maxHeight); 
   
   BOOST_LOG_TRIVIAL(info) << "There are " << num_trees << " trees";
   BOOST_LOG_TRIVIAL(info) << "The average height of a tree is " << avg_height;
   BOOST_LOG_TRIVIAL(info) << "The number of trees above the max height is  "
       << num_above;
   BOOST_LOG_TRIVIAL(info) << "The number of trees below the min size is  "
       << num_below;

   treesBelowSize.push_back( num_below );

   vector< uint64_t > nSCC;
   uint64_t NSCC;
   
   nCCs.push_back( Graph.numberConnectedComponents( kmers, NSCC ) );
   nSCC.push_back( NSCC );
   
   BOOST_LOG_TRIVIAL(info) << "Remove edges test ...";

   // How many random edges will be randomly removed
  
   unsigned remove_edge_count = min( 0.1*edgemers.size(), 1000000.0 ); 
   unordered_set<kmer_t> removed_edges;

   BOOST_LOG_TRIVIAL(info) << "Removing " << remove_edge_count << " random edges ...";

   double t_elapsed = removeRandomEdges(remove_edge_count, Graph, edgemers, removed_edges);

   double avgRemovalTime = t_elapsed / remove_edge_count;

   BOOST_LOG_TRIVIAL(info) << "Avg. removal time: " << avgRemovalTime << " s";

   /**
    * tree height tests after edge removal
    */
   BOOST_LOG_TRIVIAL(info) << "Tree height tests after edge removals ...";

 
   Graph.getTreeData(num_trees, avg_height, num_above, num_below); 

   BOOST_LOG_TRIVIAL(info) << "There are " << num_trees << " trees";
   BOOST_LOG_TRIVIAL(info) << "The average height of a tree is " << avg_height;
   BOOST_LOG_TRIVIAL(info) << "The number of trees above the max height is  "
      << num_above;
   BOOST_LOG_TRIVIAL(info) << "The number of trees below the min size is  "
      << num_below;
   treesBelowSize.push_back( num_below );
   nCCs.push_back( Graph.numberConnectedComponents( kmers, NSCC ) );
   nSCC.push_back( NSCC );
   
   BOOST_LOG_TRIVIAL(info) << "Add edges test (adding those edges back in) ...";

   unordered_set<kmer_t>::iterator removed_edge_it;
   kmer_t prefix;
   kmer_t suffix;
   t_elapsed = 0.0;
   for (removed_edge_it = removed_edges.begin(); removed_edge_it != removed_edges.end();
      removed_edge_it++) {

      Graph.split_edge(*removed_edge_it, prefix, suffix);
      clock_t t_start = clock();
      Graph.newDynamicAddEdge(prefix, suffix);
      t_elapsed += double (clock() - t_start) / CLOCKS_PER_SEC;
   }

   double avgAddTime = t_elapsed / removed_edges.size();
   
   BOOST_LOG_TRIVIAL(info) << "Avg. Addition time: " << avgAddTime << " s";

   /**
    * tree height tests after edge addition
    */
   BOOST_LOG_TRIVIAL(info) << "Tree height tests after edge additions ...";

  
   Graph.getTreeData(num_trees, avg_height, num_above, num_below ); 

   BOOST_LOG_TRIVIAL(info) << "There are " << num_trees << " trees";
   BOOST_LOG_TRIVIAL(info) << "The average height of a tree is " << avg_height;
   BOOST_LOG_TRIVIAL(info) << "The number of trees above the max height is  "
      << num_above;
   BOOST_LOG_TRIVIAL(info) << "The number of trees below the min size is  "
      << num_below;

   treesBelowSize.push_back( num_below );
   nCCs.push_back( Graph.numberConnectedComponents( kmers, NSCC ) );
   nSCC.push_back( NSCC );
   
   BOOST_LOG_TRIVIAL(info) << "Query testing...";

   std::random_device rd;
   std::mt19937 gen(rd());
   std::uniform_int_distribution<u_int64_t> sample_dis(0, pow(2, 2*k) - 1);
   size_t nQuery = 1000000;
   kmer_t kk;
   t_elapsed = 0.0;
   clock_t t_start;
   for (unsigned i = 0; i < nQuery; ++i) {
     //generate a random k-mer
     kk = sample_dis( gen ) ;
     //     BOOST_LOG_TRIVIAL( info ) << "Testing k-mer: " << get_kmer_str( kk, k ) << ' ' <<
     t_start = clock();
     Graph.detect_membership(kk);
	  t_elapsed += double(clock() - t_start);
   }
   t_elapsed = t_elapsed/ CLOCKS_PER_SEC;

   double avgQueryTime = t_elapsed / nQuery;
   BOOST_LOG_TRIVIAL(info) << "Average query time: " << avgQueryTime << " s";

   // Now we do the average query time but with only kmers that are in the graph

   nQuery = min(0.1*kmers.size(), 1000000.0);

   BOOST_LOG_TRIVIAL(info) << "Query testing for " << nQuery
			   << " kmers in the graph...";
   t_elapsed = randomMembership(nQuery, Graph, kmers);
   avgQueryTime = t_elapsed/nQuery;

   BOOST_LOG_TRIVIAL(info) << "Average query time for kmers that are in the graph: "
			   << avgQueryTime << " s";

//   for (i = kmers.begin(); i != kmers.end(); ++i) {
//		BOOST_LOG_TRIVIAL(debug) << 
//		BOOST_LOG_TRIVIAL(debug) << Graph.inefficient_detect_membership( *i ) << ' ' << Graph.detect_membership( *i );
//		if (!Graph.detect_membership( *i )) {
// 			BOOST_LOG_TRIVIAL(fatal) << "Membership test failed.";
// 			exit(1);
//		}
//	}

   //   BOOST_LOG_TRIVIAL(info) << "Membership tests passed!";

   // ////Test time to compute hash function
   // nQuery = 1000;
   // t_elapsed = 0.0;
   // t_start = clock();
   // for (unsigned i = 0; i < nQuery; ++i) {
   //   //generate a random k-mer
   //   kk = sample_dis( gen ) ;
   //   Graph.getHash(kk);
   // }
   // t_elapsed += double (clock() - t_start) / CLOCKS_PER_SEC;

   // double avgHashTime = t_elapsed / nQuery;
   // BOOST_LOG_TRIVIAL(info) << "Average hash time: " << avgHashTime << " s";

   return 0;
}

