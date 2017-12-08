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




#ifndef FORMATUTIL_CPP
#define FORMATUTIL_CPP

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#define BOOST_LOG_DYN_LINK 1
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <chrono>
#include <unordered_set>

using namespace std;

typedef uint64_t kmer_t;

unsigned access_kmer( kmer_t mer, unsigned k, unsigned i ) {
  mer = mer >> 2*(k - i - 1);
  kmer_t mask = static_cast<kmer_t>(3);
  mer = mask & mer;
  return static_cast<unsigned>(mer); 
}

string get_kmer_str( kmer_t mer, unsigned k) {

  string kmer_str = "";

  for (unsigned i = 0; i < k; ++i) {

    unsigned kk = access_kmer( mer, k, i );

    switch( kk ) {
      case 0:
        kmer_str += 'A';
        break;
      case 1:
        kmer_str += 'C';
        break;
      case 2:
        kmer_str += 'G';
        break;
      case 3:
        kmer_str += 'T';
        break;
    }
  }

  return kmer_str;
}

void print_kmer( kmer_t mer, unsigned k,ostream& os) {

  os << get_kmer_str(mer, k);

}

/*
 * Sets the i'th position of a mer of length k as indicated by character c
 * c \in {A,C,G,T}
 */
void set_kmer( kmer_t& mer, unsigned k, unsigned i, char c ) {
  //clear i-th position
  kmer_t op = 3; //0...011
  op = op << 2*(k - i - 1);  //11 in i'th spot, zeros elsewhere
  op = ~op;      //00 in i'th spot, ones elsewhere
  mer = mer & op; //i'th position of mer cleared.

  //set i'th position
  kmer_t val;
  switch( c ) {
    case 'A':
      val = 0;
      break;
    case 'C':
      val = 1;
      break;
    case 'G':
      val = 2;
      break;
    case 'T':
      val = 3;
      break;
  }

  val = val << 2*(k - i - 1);  //correct bits in i'th spot, zeros elsewhere
  mer = mer | val;
}

kmer_t mer_string_to_binary( string& r, unsigned& i, unsigned& K ) {
   kmer_t mer = 0;

   kmer_t val;
   for (unsigned j = 0; j < K; ++j) {
      //set j'th position
      char c = r[i + j];
      switch( c ) {
      case 'A':
	 val = 0;
	 break;
      case 'C':
	 val = 1;
	 break;
      case 'G':
	 val = 2;
	 break;
      case 'T':
	 val = 3;
	 break;
      default:
	 cerr << "unexpected character: " << c 
	      << ", replacing with A" << endl;
	 val = 0;
	 break;
      }

      val = val << 2*(K - j - 1);  //correct bits in j'th spot, zeros elsewhere
      mer = mer | val;
   }

   return mer;
}


void getKmers(vector< string >& reads, unsigned K, unordered_set<kmer_t>& kmers) {
   BOOST_LOG_TRIVIAL(info) << "Erasing all 'N' characters...";
   kmers.clear();
   for (auto r : reads) {
      r.erase (std::remove(r.begin(), r.end(), 'N'), r.end());
      if (r.size() < K) continue;
      for (unsigned i = 0; i < r.size() - K + 1; i++) {
	 kmer_t kmer_bin = mer_string_to_binary(r, i, K);
	 kmers.insert(kmer_bin);
      }
   }
}

void parseFasta( vector< string >& result, string filename ) {
   result.clear();
   ifstream ifile( filename.c_str() );
   string line;
   string read;
   while ( getline( ifile, line ) ) {
      if (line[0] == ' ') {
	 //blank line, skip
	 continue;
      }
      if (line[0] == '>') {
	 //if a read was in progress, finish it
	 if (read.size() > 0)
	    result.push_back( read );
	 //Start a new read
	 read.clear();
      } else {
	 //add this line to the read
	 read.append( line );
      }
   }

   //don't forget the last read
   if (read.size() > 0)
      result.push_back( read );
   
   ifile.close();
}

/* Test if a file exists */
bool file_exists(const std::string &name)
{
  if (FILE *file = fopen(name.c_str(), "r"))
  {
    fclose(file);
    return true;
  }
  else
  {
    return false;
  }
}

void write_to_bin(string outfile, unordered_set<kmer_t> &kmers)
{
  BOOST_LOG_TRIVIAL(info) << "Writing binary k-mers to file: " << outfile;
  ofstream binfile(outfile.c_str(), ios::out | ios::binary);

  kmer_t *kdata = new kmer_t[kmers.size()];
  auto it = kmers.begin();
  unsigned i = 0;
  while (it != kmers.end())
  {
    kdata[i] = *it;
    ++it;
    ++i;
  }

  binfile.write((char *)kdata, sizeof(kmer_t) * kmers.size());

  delete[] kdata;

  binfile.close();

  BOOST_LOG_TRIVIAL(info) << "Verifying file " << outfile << "...";
  vector<kmer_t> kmer_vec(kmers.size(), 0);
  ifstream ifile(outfile.c_str(), ios::in | ios::binary);
  ifile.read(reinterpret_cast<char *>(kmer_vec.data()), sizeof(kmer_t) * kmers.size());
  ifile.close();
  auto it2 = kmer_vec.begin();
  it = kmers.begin();
  while (it2 != kmer_vec.end())
  {
    if (*it2 != *it)
    {
      BOOST_LOG_TRIVIAL(fatal) << "File corrupted...";

      exit(1);
    }

    ++it2;
    ++it;
  }

  BOOST_LOG_TRIVIAL(info) << "File has been successfully verified.";
}

void read_from_bin(string filename, unordered_set<kmer_t> &out)
{
  BOOST_LOG_TRIVIAL(info) << "Reading mers from binary file " << filename << " ...";

  ifstream ifile(filename.c_str(), ios::in | ios::binary);
  ifile.seekg(0, ifile.end);
  size_t nbytes = ifile.tellg();
  ifile.seekg(0, ifile.beg);

  BOOST_LOG_TRIVIAL(info) << "File size (Gb): " << nbytes / (1024.0*1024.0*1024.0);
  
  vector<kmer_t> kmer_vec(nbytes / sizeof(kmer_t), 0);

  ifile.read(reinterpret_cast<char *>(kmer_vec.data()), nbytes);
  ifile.close();

  unordered_set<kmer_t> kmers(kmer_vec.begin(), kmer_vec.end());

  out.swap(kmers);

  BOOST_LOG_TRIVIAL(info) << "Input finished, " << kmer_vec.size() << " mers read.";
}

/*
 * Creates a (fake) Fasta file from a set of reads
 *
 */
void writeFasta( string filename, vector< string >& reads ) {
   ofstream ofile( filename.c_str() );

   for (size_t i = 0; i < reads.size(); ++i) {
      ofile << ">Somethingsomething Location, length=" << reads[i].size() << endl;
      ofile << reads[i] << endl;
   }

   ofile.close();
}

void get_kmers_fasta_or_bin(string fasta_file, unsigned k, unordered_set<kmer_t> &out)
{
  // fasta filename
  string filename = fasta_file;

  // Check if the k-mers from this file have already been stored in binary format
  string kfile = filename.substr(0, filename.find_last_of('.')) + to_string(k) + ".bin";

  if (!file_exists(kfile))
  {
    BOOST_LOG_TRIVIAL(info) << "Binary file for " << k << "-mers does not exist ...";
    // get reads from file
    BOOST_LOG_TRIVIAL(info) << "Parsing fasta ...";
    vector<string> reads;
    parseFasta(reads, filename);
    getKmers(reads, k, out);

    //save to binary file for future use
    write_to_bin(kfile, out);
  }
  else
  {
    //k-mers are already stored in the binary format!
    read_from_bin(kfile, out);
  }
}

/*
 * This function checks to see if binary information from k, k+1 mers
 * is already present. If so, it simply loads the binary format.
 * If not, it parses the fasta file and writes the k, k+1-mer binary
 * files for future use.
 */

void handle_mers(string fasta_file, unsigned k,
                 unordered_set<kmer_t> &kmers,   //k mers from fasta
                 unordered_set<kmer_t> &edgemers //k + 1 mers from fasta
                 )
{

  get_kmers_fasta_or_bin(fasta_file, k, kmers);
  get_kmers_fasta_or_bin(fasta_file, k + 1, edgemers);
}

#endif
