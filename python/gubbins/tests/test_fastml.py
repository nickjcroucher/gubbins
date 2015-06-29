#! /usr/bin/env python3
# encoding: utf-8

"""
Tests if we can detect which version of fastml is running so we can choose the correct model
"""

import unittest
import re
import os
import subprocess
from gubbins import Fastml

os.environ["PATH"] += os.pathsep + 'gubbins/tests/bin'

class TestFastml(unittest.TestCase):
    
  def test_no_fastml_installed(self):
      fastml_check = Fastml('exec_doesnt_exist')
      self.assertEqual(fastml_check.fastml_version(), None)
      self.assertEqual(fastml_check.fastml_model, None)
      self.assertEqual(fastml_check.fastml_parameters(), None)
      
  def test_fastml_3_installed(self):
      fastml_check = Fastml('dummy_fastml3')
      self.assertEqual(fastml_check.fastml_version(), 3)
      self.assertEqual(fastml_check.fastml_model,'g')
      self.assertEqual(fastml_check.fastml_parameters(), 'dummy_fastml3 -qf -b -a 0.00001 -mg')
    
  def test_fastml_2_installed(self):
      fastml_check = Fastml('dummy_fastml2')
      self.assertEqual(fastml_check.fastml_version(), 2)
      self.assertEqual(fastml_check.fastml_model, 'n')
      self.assertEqual(fastml_check.fastml_parameters(), 'dummy_fastml2 -qf -b -a 0.00001 -mn')
    
  def test_custom_fastml_2_installed(self):
      fastml_check = Fastml('dummy_custom_fastml2')
      self.assertEqual(fastml_check.fastml_version(), 2)
      self.assertEqual(fastml_check.fastml_model,'g')
      self.assertEqual(fastml_check.fastml_parameters(), 'dummy_custom_fastml2 -qf -b -a 0.00001 -mg')

if __name__ == "__main__":
  unittest.main()