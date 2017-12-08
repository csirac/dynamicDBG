#include <iostream>
#include <cmath>

/**
 * Implements an array of bits by storing integers
 */
class BitArray {

   private:

      // pointer to all our ints, which hold the bits
      unsigned* ints;

      // the number of bits in an int
      unsigned int_size;

      // the number of ints that ints points to
      unsigned num_ints;

   public:

      // Create a bit array of n bits
      BitArray(size_t num_bits)  {

	allocate( num_bits );
      }

  BitArray() {
    //default constructor
  }

  void allocate( size_t num_bits ) {
    // size of an unsigned is given in bytes
    this->int_size = 8*sizeof(unsigned);

    this->num_ints = ceil((num_bits)/((double) this->int_size));

    // create an array of that number of ints
    this->ints = new unsigned [num_ints];

    // Make sure all ints are 0
    this->clearInts();
  }

      void clearInts() {

         for (unsigned i = 0; i < this->num_ints; i++) {
            this->ints[i] = 0;
         }

      }

      // Get which int a certain bit is in
      unsigned int_num(size_t bit_num) {

         return bit_num/this->int_size;

      }

      // Get what index of its int a certain bit is in
      unsigned int_index(size_t i) {

         return i % this->int_size;

      }

      // get the ith bit
      bool get(size_t i) {

         // What int this bit_num is in
         unsigned bit_int = this->ints[this->int_num(i)];
         // What spot in that int this bit_num is
         unsigned bit_ind = this->int_index(i);

         // get this bit in the last spot of an int
         unsigned num = bit_int >> (this->int_size - bit_ind - 1);
         // Then get higher order bits off
         num &= 1;

         return (bool) num; 
      }

      // set the ith bit to v
      void set(size_t i, bool v) {
 
         // which int in our int array
         unsigned int_num = this->int_num(i);

         // what index from the right
         unsigned int_ind_right = this->int_size - this->int_index(i) - 1;

         // Set that spot
         unsigned op = 1;
         op = op << int_ind_right;
         
         if (v) {
            // set to 1
            this->ints[int_num] |=  op;
         }
         else {
            // set to 0
            op = ~op;
            this->ints[int_num] &= op;
         }

      }

      // Get the amount of bits this is taking up
      size_t total_bit_size() {

	return static_cast<size_t>(this->num_ints)*this->int_size;

      }

      // Print a string of bits to screen
      void print() {

         unsigned index = 0;

         for (int i = 0; i < this->num_ints; i++) {

             for (int j = 0; j < this->int_size; j++) {
                std::cout << this->get(index);
                index++;
             }

             std::cout << std::endl;
         }

      }

  /* writes BitArray to a binary file stream */
  void save( ostream& of ) {
    of.write( (char *) (&num_ints), sizeof( unsigned ) );
    of.write( (char *) (&int_size), sizeof( unsigned ) );

    unsigned* ptr = ints;
    for (unsigned i = 0; i < num_ints; ++i) {
      of.write( (char *)(ptr++), sizeof(unsigned) );
    }
  }

  void load( istream& of ) {
    of.read( (char *) (&num_ints), sizeof( unsigned ) );
    of.read( (char *) (&int_size), sizeof( unsigned ) );
    //    cerr << "Bitarray int_size " << int_size << endl;
    //    cerr << "Bitarray num_ints " << num_ints << endl;
    ints = new unsigned [ num_ints ];
    unsigned* ptr = ints;
    for (unsigned i = 0; i < num_ints; ++i) {
      of.read( ((char*) (ptr++)),
	       sizeof (unsigned) );
    }
    
  }
  

  ~BitArray() {

    delete[] this->ints;
  }
};



