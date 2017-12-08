#ifndef HASH_UTIL
#define HASH_UTIL

#include <ctime>
#include <iostream>
#include "../formatutil.cpp"

using namespace std;
typedef uint64_t kmer_t;

template <typename T>
T getPrime(T);
template <typename T>
bool isPrime(T);
template <typename T>
T randomNumber(const T, const T);
short baseNum(char);
unsigned access_kmer( kmer_t mer,  unsigned i );
void print_kmer( kmer_t mer, unsigned k,ostream& os);
string get_kmer_str( kmer_t mer, unsigned k);
/**Utility functions for generating our hash function*/

/*                                                                                                      
 * Finds inverse of a modulo b
 * Since we know a and b it uses Euclidean Algorithm to find the inverse of a                           
 * Note ai + bj = 1 Considering i is inverse of a and j is inverse of b                                 
 */
template< typename T >
T findInverse(T a, T b) {
   T x[3];
   T y[3];
   T quotient  = a / b;
   T remainder = a % b;

   x[0] = 0;
   y[0] = 1;
   x[1] = 1;
   y[1] = quotient * -1;

   int i = 2;
   while( (b % (a%b)) != 0 ) {
      a = b;
      b = remainder;
      quotient = a / b;
      remainder = a % b;
      x[i % 3] = (quotient * -1 * x[(i - 1) % 3]) + x[(i - 2) % 3];
      y[i % 3] = (quotient * -1 * y[(i - 1) % 3]) + y[(i - 2) % 3];

      ++i;
   }

   //x[i - 1 % 3] is inverse of a
   //y[i - 1 % 3] is inverse of b
   return x[(i - 1) % 3];
}

/**
 * Return the next prime greater than n 
 */
template <typename T>
T getPrime(T n) {

    T inc = 1;

    while (1){
        if (isPrime(n+inc)) { 
            return n+inc;
        }
        inc += 1;
    }

    return 0;
}

/**
 * Check whether n is a prime. n must be greater than 1.
 */
template <typename T>
bool isPrime(T n) {

    T root, i;

    if (n%2 == 0 || n%3 == 0)
        return false;

    root = sqrt(n);

    for (i=5; i<=root; i+=6)
    {
        if (n%i == 0)
           return false;
    }

    for (i=7; i<=root; i+=6)
    {
        if (n%i == 0)
           return false;
    }

    return true;
}

/**
 * Generate a random number between [start, end]. 
 * Note that start and end is included.
 */
template <typename T>
T randomNumber(const T start, const T end){
    double myRand = std::rand()/(1.0 + RAND_MAX); 
    T range = end - start + 1;
    T random_num = (myRand * range) + start;
    return random_num;
}

/**
 * get number in {0, 1, 2, 3} corresponding to base letter
 */
short baseNum(char base) {

   switch (base) {
      case 'A':
         return 0;
      case 'C':
         return 1;
      case 'G':
         return 2;
      case 'T':
         return 3;
   }

   return -1;
}

// Only Alan knows why this is here
#endif
