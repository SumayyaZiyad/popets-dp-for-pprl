import random   # For random hashing

import bitarray  # Efficient bit-arrays, available from:
                 # https://pypi.org/project/bitarray/

PAD_CHAR = chr(1)   # Used for q-gram padding

# =============================================================================

class BLIP():
  """BLoom-and-flIP (BLIP) hardening was proposed and used in PPRL by:
       - R. Schnell and C. Borgs, Randomized response and balanced Bloom
         filters for privacy preserving record linkage, Workshop on Data
         Integration and Applications, held at ICDM, Barcelona, 2016.

     BLIP was originally proposed as a non-interactive differentially
     private approach to randomize BFs in the context of privacy-preserving
     comparisons ofuser profiles in social networks by:
       - M. Alaggan, S. Gambs, and A.M. Kermarrec, BLIP: non-interactive
         differentially-private similarity computation on bloom filters,
         Symposium on Self-Stabilizing Systems, 2012.

     Note that this class implements both bit flipping methods proposed
     by Alaggan et al. and Schnell and Borgs.
  """

  # ---------------------------------------------------------------------------
  def __init__(self, sel_method='sch', blip_prob=0.5, random_seed=42):
    """Initialise the BLIP hardening class by providing the required
       parameters.

       Input arguments:
         - sel_method      The method of which bit flipping method to be used
                           in the hardening technique, which can either be
                           'ala' (based on the method proposed by
                           Alaggan et al., 2012) or
                           'sch' (Schnell and Borgs, 2016).

         - blip_prob       The probability value that used to flip the bit
                           values in certain bit positions in Bloom filters
                           based on differential privacy characteristics.

         - random_seed     The value used to seed the random generator used to
                           generate a random value. If no random
                           shuffling should be done set the value of this
                           argument to None. Default value is set to 42.

       Output:
         - This method does not return anything.
    """

    self.type = 'BLP'  # To identify the hardening method

    # Initialise the class variables
    assert sel_method in ['ala','sch'], sel_method
    assert 0 <= blip_prob <= 1, blip_prob

    self.sel_method  = sel_method
    self.random_seed = random_seed
    self.blip_prob   = blip_prob

  # ---------------------------------------------------------------------------

  def harden_bf(self, bf):
    """Harden the provided Bloom filter by flipping bits in certain
       positions.

       Input arguments:
         - bf  A Bloom filter assumed to have its bits set from an encoded
               q-gram set.

       Output:
         - blip_bf  The new Bloom filter after bit flipping has been applied.
    """

    bf_len = len(bf)

    # Initialise bitarray for a new Bloom filter
    #
    blip_bf = bitarray.bitarray(bf_len)

    for pos in range(bf_len):

      rand = random.random()
      new_bit = 0

      if rand <= self.blip_prob:
        # check if the rand value is at least flip probability

        if self.sel_method == 'ala':  # flip bits according to Alaggan et al.
          if bf[pos] == 0:
            new_bit = 1
          else:
            new_bit = 0

        elif self.sel_method == 'sch': # flip bits according to Schnell et al.
          bit_val = random.choice([0,1])
          new_bit = bit_val

      else: # no flipping is required
        new_bit = bf[pos]

      # Set the current bit position according to bit flipping
      #
      blip_bf[pos] = new_bit

    assert len(bf) == len(blip_bf)

    return blip_bf

# =============================================================================
# End.
