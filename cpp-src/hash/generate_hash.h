#ifndef GEN_HASH
#define GEN_HASH
#define NDEBUG

#include "BooPHF.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <unordered_set>
#include <boost/multiprecision/cpp_int.hpp>

using namespace boost::multiprecision;
#include "HashUtil.cpp"

using namespace std;

typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t  > boophf_t;

//For now, this is our big int class. If we need more than 1024 bits, we can increase
//LARGE_BITS should be a power of 2
#define LARGE_BITS 1024
typedef number<cpp_int_backend< LARGE_BITS, LARGE_BITS, unsigned_magnitude, checked, void> > largeUnsigned;

/**
 * Take in a set of k-mers, generate a hash function
 */
class generate_hash {

 public:
   u_int64_t n_kmer_orig; //the number of k-mers not considering deletions/insertions
	u_int64_t max_hash; // the max value this hash function has ever been able to take
							  // used to produce and check for invalid hash values
   unsigned k_kmer; //the lengths of the k-mers (max 32)

   /* 
    * Will store the precomputed powers of 
    * Needs to be a vector will have k distinct powers
    * Using 128 bit type to prevent overflow
    * Stores powers in order r^1, r^2, ..., r^k
    */
   vector< largeUnsigned > powersOfR;
   vector< u_int64_t > powersOfRModP;

   //image of k-mers through Karp-Rabin hash function
   //discarded after use
   std::unordered_set<u_int64_t> KRHash; 

	// A map of kmers to their hash values for newly added
	// nodes. A hack until we have a dynamic hash function
	// implementation
	std::map<kmer_t, u_int64_t> new_nodes; 

   boophf_t* bphf; //MPHF we will generate

   u_int64_t r; // the base for our Karp-Rabin Hash function
   u_int64_t rinv; //inverse of r modulo Prime
   u_int64_t Prime; // the prime for our Karp-Rabin Hash function

   const static short sigma = 4; // alphabet size


	// TODO	
   void save( ostream& of ) {
     of.write ( (char*) &n_kmer_orig, sizeof( u_int64_t ) );
     of.write ( (char*) &max_hash, sizeof( u_int64_t ) );
     of.write ( (char*) &k_kmer, sizeof( unsigned ) );
     of.write ( (char*) &r, sizeof( u_int64_t ) );
     of.write ( (char*) &rinv, sizeof( u_int64_t ) );
     of.write ( (char*) &Prime, sizeof( u_int64_t ) );

     cerr << "writing BBHash..." << endl;
     bphf->save( of );
   }

	// TODO
   void load( istream& of ) {
     of.read ( (char*) &n_kmer_orig, sizeof( u_int64_t ) );
     of.read ( (char*) &max_hash, sizeof( u_int64_t ) );
     of.read ( (char*) &k_kmer, sizeof( unsigned ) );
     of.read ( (char*) &r, sizeof( u_int64_t ) );
     of.read ( (char*) &rinv, sizeof( u_int64_t ) );
     of.read ( (char*) &Prime, sizeof( u_int64_t ) );
     //     cerr << "n k r Prime " << n_kmer << ' ' << k_kmer << ' ' << r << ' ' << Prime << '\n';
     
     //     this->bphf = new boomphf::mphf<u_int64_t, hasher_t>(n_kmer, KRHash_vec, 4, 2.0, true, false);
     bphf = new boomphf::mphf< u_int64_t, hasher_t >();
     //     cerr << "loading BBHash..." << endl;
     
     bphf->load( of );
     //     cerr << "done\n";
     precomputePowers_mod();
   }
   
   /**
    * Create hash function out of n k-mers of length k
    */
   generate_hash( unordered_set< kmer_t >& kmers , u_int64_t n , unsigned k ) {
      //  std::srand(std::time(NULL));
      std::srand( 0 );

      construct_hash_function( kmers,  n,  k );
   }

   generate_hash() {
      //default constructor
      std::srand( 0 );
   }

   generate_hash( istream& of ) {
     std::srand( 0 );
     load( of );
   }

   void construct_hash_function( unordered_set< kmer_t >& kmers , u_int64_t n , unsigned k ) {
      this->n_kmer_orig = n; // number of k-mers
		this->max_hash = n-1; // max possible hash value
      this->k_kmer = k; // length of each k-mer

      BOOST_LOG_TRIVIAL(info) << "Constructing the hash function ...";
      build_KRHash(kmers); // build KR hash function
      build_minimalPerfectHash(); // build minimal perfect hash function
      
      //compute inverse of r modulo prime
      int256_t kr = findInverse( static_cast< int256_t > (r), static_cast< int256_t > (Prime) );

      while ( kr < 0 )
	 kr = kr + Prime;

      if (kr > Prime)
	 kr = kr % Prime;

      if ( (kr * r) % Prime == 1 ) {
	 BOOST_LOG_TRIVIAL(info) << "r-inverse correctly computed modulo prime.";
      } else {
	 BOOST_LOG_TRIVIAL(fatal) << "r-inverse incorrectly computed modulo prime!";
	 exit(1);
      }
	  
      rinv = static_cast< u_int64_t >( kr );
   }

   /*
    * Once we know k (k_kmer) and r
    * we can precompute the powers of r
    */
   void precomputePowers() {
      largeUnsigned ri;
      powersOfR.clear();
      //NEED 1 to k. Not 0 to (k - 1)
      for (unsigned i = 1; i <= k_kmer; ++i) {
	 ri = mypower( r, i );
	 powersOfR.push_back( ri );
      }
   }

   void precomputePowers_mod() {
      u_int64_t ri;
      powersOfRModP.clear();
      //NEED 1 to k. Not 0 to (k - 1)
      for (unsigned i = 1; i <= k_kmer; ++i) {
	 ri = mypower_mod( r, i );
	 powersOfRModP.push_back( ri );
      }
   }

   bool verifyPowers() {
      for ( unsigned i = 0;
	    i < powersOfR.size();
	    ++i ) {
	 u_int64_t res = static_cast< u_int64_t >( powersOfR[ i ] % Prime );
	 if (res != powersOfRModP[ i ]) {
	    BOOST_LOG_TRIVIAL( fatal ) << "Hash computations overflowing.";
	    exit(1);
	 }
      }

      BOOST_LOG_TRIVIAL( info ) << "Powers of r modulo p are computed correctly.";
   }

	/**
	 * Add a node to the hash function
	 * HACK: For now, we keep the new node in a map from kmer to its hash val
	 * Don't actually have a dynamic hash function
	 * Returns new hash value
	 */
	u_int64_t add_node(const kmer_t& new_kmer)
	{

		assert(this->new_nodes.find(new_kmer) == this->new_nodes.end());

		// check if this node was previously deleted

		//auto it = this->deleted_nodes.find(this->f(new_kmer));
		
		//if (it != this->deleted_nodes.end()) {

			// this was a new node
		//	this->deleted_nodes.erase(it);
		//}

		// Add it to the new nodes map with the next available hash value
		u_int64_t new_hash = this->max_hash + 1;
		this->max_hash++;

		this->new_nodes.insert(make_pair(new_kmer, new_hash));
		
		assert(this->new_nodes.find(new_kmer) != this->new_nodes.end());

		return new_hash;
	}

	//void remove_node(const kmer_t& node, const u_int64_t& hash) {
	void remove_node(const kmer_t& node) {
		

		auto it = this->new_nodes.find(node);
		
		if (it != this->new_nodes.end()) {

			// this was a new node
			this->new_nodes.erase(it);
		}
		//else {

			// one of the original nodes
		//	this->deleted_nodes.insert(hash);	
		//}
		
		assert(this->new_nodes.find(node) == this->new_nodes.end());
		
	}

	/**
	 * Checks whether this kmer has been stored as one of the new nodes
	 * And returns the hash if it has
	 * o/w returns this->max_hash + 1
	 * which is never a valid hash
	 */
	u_int64_t new_node_hash(const kmer_t& kmer) {

		if (!this->new_nodes.empty()) {

			map<kmer_t, u_int64_t>::iterator found_kmer = this->new_nodes.find(kmer);

			if (found_kmer != this->new_nodes.end()) {
				return found_kmer->second;	
			}
		}

		
		return this->max_hash + 1;

	}

	/**
	 * Returns whether the hash value is in the range of possible hash values
	 * Between 0 and this->max_hash
	 */
	bool hash_in_range(const u_int64_t& hash) {

		if ((hash < 0) || (hash > this->max_hash)) {

			return false;
		}
		else {

			return true;
		}

	}

   /**
	 * Find the hash value of a k-mer
	 * Checks added and original nodes
	 * Will still return hash value of deleted nodes, that must be
	 * checked separately.
	 */
	u_int64_t get_hash_value(const kmer_t& seq)
	{
		u_int64_t res = this->new_node_hash(seq);

		if (res == this->max_hash + 1) {
			// not an added node
			u_int64_t krv = generate_KRHash_val_mod(seq, this->k_kmer);
			res = this->bphf->lookup(krv); // still need only 64 bits for kmer_t
	
			//if (!this->isValidOriginalHash(res)) {
				// can't be valid
			//	res = -1;
			//}

		}
		
		return res;
	}


   /**
    * Find the hash value of a k-mer
    * Allows f( v ) notation
    */
   u_int64_t operator()(const kmer_t& seq) {
      return get_hash_value( seq );
   }


   // Task4: generate_KRHash_val
   // data is a k-mer
   // k is the length of the k-mer
   // r is the base
   // P is the prime
   void build_KRHash( unordered_set< kmer_t >& kmers ){

      BOOST_LOG_TRIVIAL(info) << "Constructing Karp-Rabin hash function ...";

      // prime we will mod out by
      const u_int64_t tau = 1;

      BOOST_LOG_TRIVIAL(info) << "Theoretical prime lower bound: "
			<< tau*this->k_kmer*this->n_kmer_orig*this->n_kmer_orig;

      //Having problem with overflows because of large primes
      //Even though the theoretical bound is above, let's try smaller ones.
      double smallerPrime = this->n_kmer_orig*this->n_kmer_orig;
      Prime = getPrime((u_int64_t) smallerPrime );

      BOOST_LOG_TRIVIAL(info) << "Trying prime: " << Prime;
      
      unsigned n_failures = 0;
      // keep generating new base until we find one that is injective over our k-mers
      bool f_injective;
      u_int64_t v1, v2; // holder for KRH value

      do
	 {
	    ++n_failures;
	    if ( n_failures == 5 ) {
	       BOOST_LOG_TRIVIAL(info) << "Trying a larger prime... ";
	       Prime = getPrime( Prime * 2 );
	       n_failures = 0;
	       BOOST_LOG_TRIVIAL(info) << "Trying prime: " << Prime;
	    }

	    
	    f_injective = true; //assume f is injective until evidence otherwise
	    this->r = randomNumber((u_int64_t) 1, Prime-1);
	    //Once we have a candidate base r
	    //we should avoid recomputing its powers all the time
	    
	    //	    precomputePowers();
	    precomputePowers_mod();
	    //	    verifyPowers();
	    
	    for ( unordered_set< kmer_t >::iterator
		     it1 = kmers.begin(); it1 != kmers.end();
		  ++it1 ) {
	       v1 = generate_KRHash_val_mod( *it1, k_kmer );
	       //	       v2 = generate_KRHash_val( *it1, k_kmer, Prime);
	       //	       if (v1 != v2) {
	       //		  BOOST_LOG_TRIVIAL(fatal) << "Error in Karp-rabin computations.";
	       //		  exit(1);
	       //	       }
	       
	       //BOOST_LOG_TRIVIAL(trace) << "hash of kmer: " << v;
	       if (KRHash.find(v1) == KRHash.end())
		  {
		     // this is a new value
		     KRHash.insert(v1);
		  }
	       else // not injective
		  {
			 BOOST_LOG_TRIVIAL(trace) << "Base " <<this->r << " with prime "
						  << Prime << " failed injectivity.";
			 KRHash.clear(); // clear it out and start over
			 f_injective = false;
			 break;
		      }

		}
	     } while (!f_injective);


	  BOOST_LOG_TRIVIAL(info) << "Base " << this->r << " with prime " << Prime << " is injective.";

	  return;
   }


   /*
    * Computes large powers with large ints
    */
   largeUnsigned mypower( const u_int64_t& base, unsigned exponent ) {
      largeUnsigned rvalue( 1 );
      while (exponent > 0) {
	 rvalue *= static_cast< largeUnsigned >(base);
	 --exponent;
      }

      return rvalue;
   }

   /*
    * Computes powers mod p with u_int64_t and integer exponents
    * Uses larger arithmetic type to prevent overflow
    * Returns a u_int64_t because it mods out before returning
    */
   u_int64_t mypower_mod( const u_int64_t& base, unsigned exponent ) {
      uint128_t rvalue( 1 );
      while (exponent > 0) {
	 rvalue = rvalue * base;
	 rvalue = rvalue % Prime;
	 --exponent;
      }

      return static_cast< u_int64_t >( rvalue );
   }
   
   /**
    * Given a kmer, find out its KRH
    * MODULO Prime
    * Need to use 128 bits since 4 * x might overflow
	 * Note that this shouldn't be used for newly added nodes
	 * that are not part of the main hash function
    */
   u_int64_t generate_KRHash_val_mod(const kmer_t& kmer,
				     const unsigned& k ) {

		uint128_t val = 0; // what will be the KRH value

		// go through each bp and add value
		for (unsigned i = 0; i < k; ++i) {
	 		// val += baseNum(kmer.at(i)) * pow(r, i);
	 		val = val + ((access_kmer( kmer, k, i) *
				 static_cast< uint128_t >( powersOfRModP[i] )));
	 		val = val % Prime;
		}

		return static_cast<u_int64_t>(val);
   }

   
   /**
    * Given a kmer, find out its KRH using base r and prime P
    */
   u_int64_t generate_KRHash_val(const kmer_t& kmer,
				 const unsigned& k,
				 const u_int64_t& P){
  
		//  BOOST_LOG_TRIVIAL(trace) << "Generating KRHash val";
		//use 128 bits to prevent overflow
		largeUnsigned val = 0; // what will be the KRH value

		// go through each bp and add value
		for (unsigned i = 0;
		i < k;
		++i) {
	 // val += baseNum(kmer.at(i)) * pow(r, i);
			 val +=
		  static_cast< largeUnsigned > ( access_kmer( kmer, k, static_cast<unsigned>(i)) ) *
		  powersOfR[i]; //powersOfR[i] = r^{i + 1}
		}

		val = val % P;

		return static_cast<u_int64_t>(val);

   }

   /**
    * Given a kmer, find out its KRH using base r and prime P
    * HOWEVER: does not mod out by P. So returns large unsigned
    */
   largeUnsigned generate_KRHash_raw(const kmer_t& kmer,
				     const unsigned& k ) {
      //  BOOST_LOG_TRIVIAL(trace) << "Generating KRHash val";
      //use 128 bits to prevent overflow
      largeUnsigned val = 0; // what will be the KRH value

      // go through each bp and add value
      for (unsigned i = 0;
	   i < k;
	   ++i) {
	 // val += baseNum(kmer.at(i)) * pow(r, i);
	       val +=
		  static_cast< largeUnsigned > ( access_kmer( kmer, k, static_cast<unsigned>(i)) ) *
		  powersOfR[i]; //powersOfR[i] = r^{i + 1}
      }

      return val;
   }

   /*
    * This function takes as input a Karp-Rabin value (KR_val)
    * Then updates it by subtracting the value from 'first' character source kmer
    * Then dividing by r (at this point, it has shifted last k-1 characters up)
    * Finally adding the last term corresponding to the 'last' character
    *
    * target k-mer is OUT neighbor of source k-mer
    *
 	 * Note that this shouldn't be used for newly added nodes
	 * that are not part of the main hash function
      */
   void update_KRHash_val_OUT
      ( largeUnsigned& KR_val,       //KR hash of source kmer
	const unsigned& first,   //character at front of source k-mer
	const unsigned& last ) { //last character in target k-mer
      //   BOOST_LOG_TRIVIAL(debug) << "Updating a KR value by OUT...";
      //   BOOST_LOG_TRIVIAL(debug) << "First of source: " << first;
      //   BOOST_LOG_TRIVIAL(debug) << "Last of target: " << last;

      //   largeUnsigned before_div = KR_val;
      //   BOOST_LOG_TRIVIAL(debug) << "Value before division: " << before_div;

      //   BOOST_LOG_TRIVIAL(debug) << "Division check: " << (KR_val * static_cast< largeUnsigned >( r ) == before_div );

      //   KR_val = before_div;
      //largeUnsigned q;
      //   largeUnsigned rem;
      //   divide_qr( KR_val, static_cast< largeUnsigned >( r ),  q, rem );

      //   BOOST_LOG_TRIVIAL(debug) << "Division check 2: " << (q * static_cast< largeUnsigned >( r ) == before_div );

      //   BOOST_LOG_TRIVIAL(debug) << "Remainder: " << rem;
      KR_val = KR_val / static_cast< largeUnsigned >( r );
      KR_val = KR_val - static_cast< largeUnsigned >( first );
      KR_val = KR_val + static_cast< largeUnsigned >(last) * powersOfR[ k_kmer - 1 ]; // last * r^k
   }

   /*
    * This function takes as input a Karp-Rabin value (KR_val)
    *
    * target k-mer is IN neighbor of source k-mer
 	 * Note that this shouldn't be used for newly added nodes
	 * that are not part of the main hash function
      */
   void update_KRHash_val_IN
      ( largeUnsigned& KR_val,       //KR hash of source kmer
	const unsigned& first,   //character at front of target k-mer
	const unsigned& last ) { //last character in source k-mer
      //   BOOST_LOG_TRIVIAL(debug) << "Updating a KR value by IN...";
      //   BOOST_LOG_TRIVIAL(debug) << "First of target: " << first;
      //   BOOST_LOG_TRIVIAL(debug) << "Last of source: " << last;

      KR_val = KR_val - last * powersOfR[ k_kmer - 1 ]; // last * r^k
      KR_val = KR_val * r;
      KR_val = KR_val + first * r;
   }

   typedef number<cpp_int_backend<128, 128, unsigned_magnitude, checked, void> > moduloInt;
   u_int64_t update_KRHash_val_OUT_mod
      ( const u_int64_t& KR_val,       //KR hash of source kmer (mod P)
	const unsigned& first,   //character at front of source k-mer
	const unsigned& last ) { //last character in target k-mer
      //	   BOOST_LOG_TRIVIAL(debug) << "Updating a KR value by OUT(mod)...";
      //	   BOOST_LOG_TRIVIAL(debug) << "First of source: " << first;
      //	   BOOST_LOG_TRIVIAL(debug) << "Last of target: " << last;
      moduloInt rinv = this->rinv;
      moduloInt Prime = this->Prime;
      moduloInt kr = KR_val;
      moduloInt llast = last;
      moduloInt ffirst = first;
      moduloInt rk = powersOfRModP[ k_kmer - 1 ];
      moduloInt four = 4;
	   
      moduloInt sub_val = four*Prime - ffirst * r;
      kr = (kr + sub_val );
      kr = (kr * rinv);
      kr = (kr + llast * rk);

	   //	   moduloInt q, rem;
	   //	   divide_qr( kr, Prime, q, rem );
	   //	   HashInt kr2 = integer_modulus( kr, this->Prime );
      kr = kr % Prime;

	   //	   while (kr > Prime)
	   //	      kr = kr - Prime;

	   //	   kr = kr - (kr / Prime)*Prime;
	   
      return static_cast< u_int64_t >( kr ); 
   }

   /*
    * This function takes as input a Karp-Rabin value (KR_val) modulo Prime
    *
    * target k-mer is IN neighbor of source k-mer
      */
   u_int64_t update_KRHash_val_IN_mod
      ( const u_int64_t& KR_val,       //KR hash of source kmer
	const unsigned& first,   //character at front of target k-mer
	const unsigned& last ) { //last character in source k-mer
	   //	   BOOST_LOG_TRIVIAL(debug) << "Updating a KR value by IN(mod)...";
	   //	   BOOST_LOG_TRIVIAL(debug) << "First of target: " << first;
	   //	   BOOST_LOG_TRIVIAL(debug) << "Last of source: " << last;

      moduloInt r = this->r;
      moduloInt Prime = this->Prime;
      moduloInt kr = KR_val;
      moduloInt llast = last;
      moduloInt ffirst = first;
      moduloInt rk = powersOfRModP[ k_kmer - 1 ];
      moduloInt four = 4;
	   
      moduloInt sub_val = four*Prime - llast*rk;

      kr = (kr + sub_val) ; // last * r^k
      kr = (kr * r) ;
      kr = (kr + ffirst * r) ;

      //	   moduloInt q, rem;
      //	   divide_qr( kr, Prime, q, rem );
      //	   HashInt kr2 = integer_modulus( kr, this->Prime );
      kr = kr % Prime;

      //	   while (kr > Prime)
      //	      kr = kr - Prime;
      
      //	   kr = kr - (kr / Prime)*Prime;
	   
      return static_cast< u_int64_t > (kr) ;
	   
   }

   
   /*
    * Looks up the minimal perfect hash value, given the Karp-Rabin (raw value)
    * KR raw value means not modded out by the prime yet.
	 * If in added nodes, just returns the hash for that kmer. Doesn't use KR.
        */
   u_int64_t perfect_from_KR( const kmer_t& kmer, const largeUnsigned& KR_val ) {

		u_int64_t hash = this->new_node_hash(kmer);

		if (hash == this->max_hash + 1) {
			// not a new node
      	largeUnsigned KR2 = KR_val % Prime;
      	hash = this->bphf->lookup( static_cast< u_int64_t >(KR2) );
		}

		return hash;

   }

   /*
    * Looks up the minimal perfect hash value, given the Karp-Rabin (value modulo Prime)
	 * If in added nodes, just returns the hash for that kmer. Doesn't use KR.
    */
   u_int64_t perfect_from_KR_mod( const kmer_t& kmer, const u_int64_t& KR_val ) {

		u_int64_t hash = this->new_node_hash(kmer);
      
		if (hash == this->max_hash + 1) {
			// not a new node
			hash = this->bphf->lookup( KR_val );
		}

		return hash;
   }

   /**
    * Build a minimal perfect hash function on the set of integers that our kmers are
    * mapped to via KRH
    */
   void build_minimalPerfectHash(){

      //std::sort(KR_hash_val, KR_hash_val+n_kmer);
      //u_int64_t jj = 0;
      //for (int ii = 1; ii < n_kmer; ii++) {
      //    if (KR_hash_val[ii] != KR_hash_val[jj])
      //        KR_hash_val[++jj] = KR_hash_val[ii];
      //}
      //printf("Found %lli duplicated items from KR_hash_val.  \n", n_kmer-(jj + 1) );

      //auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(KR_hash_val), static_cast<const u_int64_t*>(KR_hash_val+n_kmer));

      //            bphf = new boomphf::mphf<u_int64_t, hasher_t>(n_kmer, data_iterator, nthreads, gammaFactor);

      BOOST_LOG_TRIVIAL(info) << "Building minimal perfect hash function ...";

      std::vector<u_int64_t> KRHash_vec (KRHash.begin(),
					 KRHash.end());

      // MPHF for our KRHash function values
      this->bphf = new boomphf::mphf<u_int64_t, hasher_t>(this->n_kmer_orig, KRHash_vec, 4, 2.0, true, false);

      BOOST_LOG_TRIVIAL(info) << "Minimal perfect hash function created with "
			      << (float) (bphf->totalBitSize())/this->n_kmer_orig << " bits per element.";

      //      KRHash.clear();
   }


};

#endif
