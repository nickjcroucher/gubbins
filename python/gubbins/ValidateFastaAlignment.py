import os
import re
import sys
from Bio import AlignIO
from collections import Counter
from Bio.Align import MultipleSeqAlignment
# Gubbins imports
from gubbins.utils import process_sequence_names

class ValidateFastaAlignment(object):

    def __init__(self, input_filename):
      self.input_filename = input_filename

    def is_input_fasta_file_valid(self):
      try:
          if not self.does_each_sequence_have_the_same_length():
              print("Each sequence must be the same length")
              return False
          if not self.are_sequence_names_unique():
              print("All sequence names in the fasta file must be unique")
              return False
          if not self.does_each_sequence_have_a_name_and_genomic_data():
              print("Each sequence must have a name and some genomic data")
              return False
      except:
          return False
      return True

    def does_each_sequence_have_a_name_and_genomic_data(self):
      with  open(self.input_filename, "r") as input_handle:
        alignments = AlignIO.parse(input_handle, "fasta")
        number_of_sequences = 0
        for alignment in alignments:
            for record in alignment:
                number_of_sequences +=1
                if record.name is None or record.name == "":
                  sys.stderr.write("Error with the input FASTA file: " + record.name + " is blank\n")
                  return False
                if record.seq is None or record.seq == "":
                  sys.stderr.write("Error with the input FASTA file: " + record.name + " is empty\n")
                  return False
                if re.search('[^ACGTNacgtn-]', str(record.seq))  != None:
                  sys.stderr.write("Error with the input FASTA file: " + record.name + " contains disallowed characters, only ACGTNacgtn- are permitted\n")
                  return False
      return True

    def does_each_sequence_have_the_same_length(self):
      try:
        with open(self.input_filename) as input_handle:
          alignments = AlignIO.parse(input_handle, "fasta")
          sequence_length = -1
          for alignment in alignments:
              for record in alignment:
                 if sequence_length == -1:
                   sequence_length = len(record.seq)
                 elif sequence_length != len(record.seq):
                   print("Error with the input FASTA file: The sequences are not of the same length, this is not an alignment: "+record.name)
                   return False
          input_handle.close()
      except:
        print("Unexpected error:", sys.exc_info()[0])
        print("Error with the input FASTA file: It is in the wrong format, check it is an alignment")
        return False
      return True

    def are_sequence_names_unique(self):
        any_modified_names = False
        with open(self.input_filename) as input_handle:
            alignment = AlignIO.read(input_handle, "fasta")
            sequence_names = []
            for record in alignment:
                # Remove disallowed characters
                if '#' in record.name or ':' in record.name:
                    record.name = process_sequence_names(record.name)
                    record.id = process_sequence_names(record.id)
                    record.description = process_sequence_names(record.namdescriptione)
                    any_modified_names = True
                # Store modified names
                sequence_names.append(record.name)
        if [k for k,v in list(Counter(sequence_names).items()) if v>1] != []:
          return False
        # Update alignment if names changed
        if any_modified_names:
            with open(self.input_filename, "w") as output_handle:
                AlignIO.write(alignment,output_handle, "fasta")
        return True
      
