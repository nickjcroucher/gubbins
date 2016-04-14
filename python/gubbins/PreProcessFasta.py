import os
import hashlib
from Bio import AlignIO
from collections import defaultdict

class PreProcessFasta(object):
  def __init__(self, input_filename, verbose = False):
      self.input_filename = input_filename
      self.verbose = verbose


  def _hash_sequences(self):
      sequence_hash_to_taxa = defaultdict(list)
      with open(self.input_filename) as input_handle:
          alignments = AlignIO.parse(input_handle, "fasta")
          for alignment in alignments:
              for record in alignment:
                  sequence_hash = hashlib.md5()
                  sequence_hash.update(str(record.seq).encode('utf-8') )
                  hash_of_sequence  = sequence_hash.digest()
                  sequence_hash_to_taxa[hash_of_sequence].append(record.id)
                   
                  if self.verbose:
                      print("Sample "+ str(record.id) + " has a hash of " +str(hash_of_sequence) )
                  
      return sequence_hash_to_taxa
          
