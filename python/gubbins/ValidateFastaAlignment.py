import os
import re
import sys
from Bio import AlignIO
from collections import Counter
from Bio.Align import MultipleSeqAlignment

class ValidateFastaAlignment(object):
    def __init__(self, input_filename):
      self.input_filename = input_filename

    def is_input_fasta_file_valid(self):
      self.check_special_chars_in_names()

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
                  print("Error with the input FASTA file: One of the sequence names is blank")
                  return False
                if record.seq is None or record.seq == "":
                  print("Error with the input FASTA file: One of the sequences is empty")
                  return False
                if re.search('[^ACGTNacgtn-]', str(record.seq))  != None:
                  print("Error with the input FASTA file: One of the sequences contains odd characters, only ACGTNacgtn- are permitted")
                  print(record.id)
                  return False
        input_handle.close()
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
                   print("Error with the input FASTA file: The sequences dont have the same lengths this isnt an alignment: "+record.name)
                   return False
          input_handle.close()
      except:
        print("Unexpected error:", sys.exc_info()[0])
        print("Error with the input FASTA file: It is in the wrong format so check its an alignment")
        return False
      return True

    def are_sequence_names_unique(self):
      with open(self.input_filename) as input_handle:
        alignments = AlignIO.parse(input_handle, "fasta")
        sequence_names = []
        for alignment in alignments:
            for record in alignment:
                sequence_names.append(record.name)

        if [k for k,v in list(Counter(sequence_names).items()) if v>1] != []:
          return False
        input_handle.close()
      return True

    def check_special_chars_in_names(self):
        # Remove any #s and :s from the isolate names
        with open(self.input_filename) as input_handle:
            alignments = AlignIO.parse(input_handle, "fasta")
            sequence_names = []
            for alignment in alignments:
                for record in alignment:
                    sequence_names.append(record.name)
            input_handle.close()

        if any(["#" in name for name in sequence_names]) | any([":" in name for name in sequence_names]):
            print("Removing #s and/or :s from sequence names")
            with open(self.input_filename, 'r') as input_handle:
                alignments_seq = AlignIO.read(input_handle, "fasta")

                for record in alignments_seq:
                    record.name = record.name.replace("#","_").replace(":","_")
                    record.id = record.id.replace("#", "_").replace(":", "_")
                    record.description = record.description.replace("#", "_").replace(":", "_")

            with open(self.input_filename, "w") as output_handle:
                AlignIO.write(alignments_seq,output_handle, "fasta")

        return True


      
