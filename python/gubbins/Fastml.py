import os
import re
import subprocess

class Fastml(object):
  def __init__(self, fastml_exec = None):
      self.fastml_exec = fastml_exec
      self.fastml_version = None
      self.fastml_model = None
      self.fastml_parameters = self.__calculate_parameters__()
      
  def __calculate_parameters__(self):
      if(self.which(self.fastml_exec) == None):
        return None
      
      if re.search('nucgtr', str(self.__run_without_options__())):
          self.fastml_version = 3
          self.fastml_model = 'g'
      else:
          self.fastml_version = 2
          self.fastml_model = 'n'
          
      return self.fastml_exec + " -qf -b -a 0.00001 -m"+self.fastml_model+" "

      
  def __run_without_options__(self):
      return subprocess.Popen(self.fastml_exec, stdout = subprocess.PIPE, shell=True).communicate()[0]
      
  def which(self,program):
      executable = program.split(" ")
      program = executable[0]
      def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
      fpath, fname = os.path.split(program)
      if fpath:
        if is_exe(program):
          return program
      else:
        for path in os.environ["PATH"].split(os.pathsep):
          exe_file = os.path.join(path, program)
          if is_exe(exe_file):
            return exe_file
    
      return None
        