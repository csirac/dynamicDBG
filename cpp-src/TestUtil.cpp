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




#define NDEBUG
#include <random>
#include <math.h>
#include <set>
#include <unordered_set>

/**
 * Functions to help test the FDBG data structure
 */

/**
 * Remove random edges from the graph
 * count is the number of edges removed
 *
 * Returns elapsed time in seconds (of only the edge removal update)
 */
double removeRandomEdges(const unsigned& count, FDBG& Graph, unordered_set<kmer_t> edgemers,
   unordered_set<kmer_t>& removed) {

   assert(edgemers.size() >= count);

   removed.clear();

   kmer_t prefix;
   kmer_t suffix;
   double time_elapsed = 0;
 
   // get random number generator
   std::random_device rd;
   std::mt19937 rand (rd());

   vector < kmer_t > vEdges( edgemers.begin(), edgemers.end() );
   std::uniform_int_distribution<u_int64_t> unif_dist (0, edgemers.size() - 1);
   while (removed.size() < count) {
      int randnum = unif_dist(rand);

      // get the prefix and suffix
      Graph.split_edge( vEdges[ randnum ], prefix, suffix);

      clock_t t_start = clock();
      if ( Graph.dynamicRemoveEdge(prefix, suffix) ) {
	 time_elapsed += double (clock() - t_start) / CLOCKS_PER_SEC;
	 removed.insert( vEdges[ randnum ] );
      }
   }

   return time_elapsed;
}

/**
 * Test the membership of random nodes in the graph
 *
 * Returns elapsed time in seconds
 */
double randomMembership(const unsigned& count, FDBG& Graph, const unordered_set<kmer_t>& kmers) {

   double time_elapsed = 0;

   // get random number generator
   std::random_device rd;
   std::mt19937 rand (rd());
   std::uniform_int_distribution<u_int64_t> unif_dist (0, kmers.size() - 1);

   vector<kmer_t> vKmers( kmers.begin(), kmers.end() );
   
   for (int i = 0; i < count; i++) {
      int randnum = unif_dist(rand);
      // auto it = kmers.begin();

      // for (int i = 0; i < randnum; i++) {
      // 	 it++;			
      // 	 assert(it != kmers.end());
      // }

      // Test membership
      clock_t t_start = clock();
      bool member = Graph.detect_membership( vKmers[ randnum ] );
      time_elapsed += double (clock() - t_start);
      assert (member);
   }

   return time_elapsed/CLOCKS_PER_SEC;
}

/**
 * Remove and then add back in random kmers along with their edges
 *
 * Returns elapsed time in seconds for removals and for additions
 */
pair<double, double> randomDynamicNodes(const unsigned& count, FDBG& Graph, const unordered_set<kmer_t>& kmers) {

   double time_elapsed_remove = 0;
   double time_elapsed_add = 0;

   // get random number generator
   std::random_device rd;
   std::mt19937 rand (rd());
   std::uniform_int_distribution<u_int64_t> unif_dist (0, kmers.size() - 1);

   vector<kmer_t> vKmers( kmers.begin(), kmers.end() );

   for (int i = 0; i < count; i++) {

      int randnum = unif_dist(rand);
		kmer_t kmer = vKmers[randnum];

		//BOOST_LOG_TRIVIAL(info) << "Removing and then adding back in " << kmer;

		// get the node's neighbors so the edges can be added back
		vector<kmer_t> neighbors;
		vector<bool> in_or_out;
		Graph.get_neighbors(kmer, neighbors, in_or_out);

      // Remove the node
		//size_t depth;
		//kmer_t root;
		//unsigned size_before = Graph.getTreeSize(kmer, depth, root);
      clock_t t_start = clock();
      Graph.removeNode(kmer);
		//Graph.isolateNode(kmer);
		//for (int i = 0; i < neighbors.size(); i++) {
		//	if (in_or_out[i]) {
		//		// was in IN
		//		Graph.dynamicRemoveEdge(neighbors[i], kmer);
		//	}
		//	else {
		//		// was in OUT
		//		Graph.dynamicRemoveEdge(kmer, neighbors[i]);
		//	}
		//}
      time_elapsed_remove += double (clock() - t_start);
      //assert (!Graph.detect_membership(kmer));

	   // Add the node along with its edges back in
      t_start = clock();
      Graph.addNode(kmer);
		//unsigned size = Graph.getTreeSize(kmer, depth, root);
      assert (Graph.detect_membership(kmer));
		//assert (size == 1);
		//assert (root == kmer);
		for (int i = 0; i < neighbors.size(); i++) {
			if (in_or_out[i]) {
				// was in IN
				Graph.newDynamicAddEdge(neighbors[i], kmer);
			}
			else {
				// was in OUT
				Graph.newDynamicAddEdge(kmer, neighbors[i]);
			}
		}
		//for (int i = 0; i < neighbors.size(); i++) {
		//	assert(Graph.getTreeSize(neighbors[i], depth, root) >= min(Graph.alpha, size_before));
		//}
      time_elapsed_add += double (clock() - t_start);
      assert (Graph.detect_membership(kmer));
   }

   return make_pair(time_elapsed_remove/CLOCKS_PER_SEC, time_elapsed_add/CLOCKS_PER_SEC);
}























//
