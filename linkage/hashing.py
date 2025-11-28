import hashlib  # A standard Python library
import random  # For random hashing

import bitarray  # Efficient bit-arrays, available from:


# https://pypi.org/project/bitarray/


# =============================================================================

class RandomHashing():
    """Random-hashing for Bloom filters was proposed and used by:
       - R. Schnell and C. Borgs, Randomized response and balanced Bloom
         filters for privacy preserving record linkage, Workshop on Data
         Integration and Applications, held at ICDM, Barcelona, 2016.
  """

    # ---------------------------------------------------------------------------

    def __init__(self, hash_funct, bf_len, num_hash_funct,
                 get_q_gram_pos=False):
        """Initialise the random-hashing class by providing the required
       parameters.

       Input arguments:
         - hash_funct      The hash function to be used to encode q-grams.
         - bf_len          The length in bits of the Bloom filters to be
                           generated.
         - num_hash_funct  The number of hash functions to be used.
         - get_q_gram_pos  A flag, if set to True then the bit positions of
                           where q-grams are hash into are returned in a
                           dictionary.

       Output:
         - This method does not return anything.
    """

        # Initialise the class variables
        #
        self.hash_funct = hash_funct

        assert bf_len > 1, bf_len
        self.bf_len = bf_len

        assert num_hash_funct > 0
        self.num_hash_funct = num_hash_funct

        assert get_q_gram_pos in [True, False]
        self.get_q_gram_pos = get_q_gram_pos

    # ---------------------------------------------------------------------------

    def hash_q_gram_set(self, q_gram_set, salt_str=None):
        """Hash the given q-gram set according to the parameter using when
       initialising the class.

       Input arguments:
         - q_gram_set  The set of q-grams (strings) to be hashed into a Bloom
                       filter.
         - salt_str    An optional string used for salting, if provided this
                       string will be concatenated with every q-gram in the
                       given q-gram set. If set to None no salting will be
                       done.

       Output:
         - bf               A Bloom filter with bits set according to the input
                            q-gram set and random hashing parameters.
         - q_gram_pos_dict  [Only returned if 'get_q_gram_pos' is set to True]
                            A dictionary which has q-grams as keys and where
                            values are sets with the positions these q-grams
                            are hashed to.
    """

        bf_len = self.bf_len  # Short-cuts
        k = self.num_hash_funct

        get_q_gram_pos = self.get_q_gram_pos
        if get_q_gram_pos:
            q_gram_pos_dict = {}

        # Initialise the Bloom filter to have only 0-bits
        #
        bf = bitarray.bitarray(bf_len)
        bf.setall(0)

        # Bloom filter length minus 1 using for position index in Bloom filter
        #
        bf_len_m1 = bf_len - 1

        # Hash the q-grams into the Bloom filter
        #
        for q_gram in q_gram_set:

            if get_q_gram_pos:
                q_gram_pos_set = set()

            if salt_str is not None:  # If a salt is given concatenate with q-gram
                q_gram = q_gram + salt_str

            # Use q-gram itself to see random number generator
            #
            hex_str = self.hash_funct(q_gram.encode("utf-8")).hexdigest()
            random.seed(int(hex_str, 16))

            for i in range(k):
                pos_i = random.randint(0, bf_len_m1)
                bf[pos_i] = 1

                if get_q_gram_pos:
                    q_gram_pos_set.add(pos_i)

            if get_q_gram_pos:
                q_gram_pos_dict[q_gram] = q_gram_pos_set

        if get_q_gram_pos:
            return bf, q_gram_pos_dict
        else:
            return bf


# =============================================================================
# Some testing code if called from the command line

if (__name__ == '__main__'):

    print('Running some tests:')
    print()

    # Define three hash functions
    #
    bf_hash_funct1 = hashlib.sha1
    bf_hash_funct2 = hashlib.md5
    bf_hash_funct3 = hashlib.sha256

    # Define Bloom filter hashing parameters
    #
    bf_len = 1000
    k = 10

    test_q_gram_set = {'he', 'el', 'll', 'lo', 'o ', ' w', 'wo', 'or', 'rl', 'ld'}

    print('  Test q-gram set:', test_q_gram_set)
    print()

    print('  Testing random hashing...'),  # - - - - - - - - - - - - - - - - - - -

    # Initialise the random hashing class
    #
    RH = RandomHashing(bf_hash_funct1, bf_len, k)

    rh_bf1 = RH.hash_q_gram_set(test_q_gram_set)
    assert len(rh_bf1) == bf_len
    assert rh_bf1.count(1) > 0

    rh_bf2 = RH.hash_q_gram_set(test_q_gram_set)
    assert len(rh_bf2) == bf_len
    assert rh_bf2.count(1) > 0

    assert rh_bf1 == rh_bf2

    rh_bf3 = RH.hash_q_gram_set(test_q_gram_set, 'test salt str')
    assert len(rh_bf3) == bf_len
    assert rh_bf3.count(1) > 0

    assert rh_bf1 != rh_bf3

    # assert dh_bf1 != rh_bf1  # Highly unlikely the two hashing approaches will
    # assert dh_bf2 != rh_bf2  # generate the same Bloom filters
    # assert dh_bf3 != rh_bf3

    # Check g-gram position dictionary
    #
    RH2 = RandomHashing(bf_hash_funct1, bf_len, k, True)

    rh_bf, rh_q_gram_pos_dict = RH2.hash_q_gram_set(test_q_gram_set)

    assert len(rh_q_gram_pos_dict) == len(test_q_gram_set)
    for (q_gram, pos_set) in rh_q_gram_pos_dict.iteritems():
        assert q_gram in test_q_gram_set
        assert len(pos_set) > 0 and len(pos_set) <= k
        for pos in pos_set:
            assert pos >= 0 and pos < bf_len

    print('OK')
    print()

# =============================================================================
# End.
