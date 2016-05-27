# encoding: utf-8
# Wellcome Trust Sanger Institute
# Copyright (C) 2013  Wellcome Trust Sanger Institute
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

import os
import sys
import subprocess
import re

class RAxMLExecutable(object):
	def __init__(self, threads,  verbose = False ):
		self.verbose = verbose
		self.threads = threads
		self.single_threaded_executables = ['raxmlHPC-AVX','raxmlHPC-SSE3','raxmlHPC']
		self.multi_threaded_executables = ['raxmlHPC-PTHREADS-AVX','raxmlHPC-PTHREADS-SSE3','raxmlHPC-PTHREADS']
		
		self.raxml_executable = self.select_executable_based_on_threads()
		self.tree_building_parameters = ' -f d -p 1 -m GTRGAMMA '
		self.internal_sequence_parameters = ' -f A -p 1 -m GTRGAMMA '
		
	def tree_building_command(self):
		command = self.raxml_executable + self.threads_parameter() + self.tree_building_parameters
		if self.verbose:
			print("Tree building command: "+command)
		return command
		
	def internal_sequence_reconstruction_command(self):
		command = self.raxml_executable + self.threads_parameter() + self.internal_sequence_parameters
		if self.verbose:
			print("Internal sequence reconstruction command: "+command)
		return command
	
	def choose_executable_from_list(self,list_of_executables):
		flags = []
		if os.path.exists('/proc/cpuinfo'):
			output = subprocess.Popen('grep flags /proc/cpuinfo', stdout = subprocess.PIPE, shell=True).communicate()[0].decode("utf-8")
			flags =  output.split()
			
		for executable in list_of_executables:
			if os.path.exists('/proc/cpuinfo'):
				if re.search('AVX', executable) and 'avx' not in flags:
					continue
				elif re.search('SSE3', executable) and 'ssse3'  not in flags:
					continue
			
			if self.which(executable) != None:
				return executable
		  
		return None
		
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
	
	def threads_parameter(self):
		if self.threads > 1:
			return " -T " + str(self.threads) + " "
		else:
			return ""
	
	def select_executable_based_on_threads(self):
		single_threaded_exec = self.choose_executable_from_list(self.single_threaded_executables)
		multi_threaded_exec =  self.choose_executable_from_list(self.multi_threaded_executables)
		
		if self.threads == 1:
			if single_threaded_exec != None:
				return single_threaded_exec
			else:
				print("Trying multithreaded version of RAxML because no single threaded version of RAxML could be found. Just to warn you, this requires 2 threads.\n")
				self.threads = 2
		
		if self.threads > 1:
			if multi_threaded_exec != None:
				return multi_threaded_exec
			else:
				sys.exit("No usable version of RAxML could be found, please ensure one of these executables is in your PATH:\nraxmlHPC-PTHREADS-AVX\nraxmlHPC-PTHREADS-SSE3\nraxmlHPC-PTHREADS\n raxmlHPC-AVX\nraxmlHPC-SSE3\nraxmlHPC")
			
