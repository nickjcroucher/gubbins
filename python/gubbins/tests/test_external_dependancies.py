#! /usr/bin/env python
# encoding: utf-8

"""
Integration testing of external dependancies. Likely to be the most brittle tests, but the most important.
"""

import unittest
import re
from gubbins import common

class TestExternalDependancies(unittest.TestCase):
  
  def test_which(self):
    # the location of ls varies depending on OS so just check end
    assert re.match('.*/ls$', common.GubbinsCommon.which('ls')) != None
    # Strip parameters
    assert re.match('.*/ls$', common.GubbinsCommon.which('ls -alrt')) != None
    assert common.GubbinsCommon.which('non_existant_program') == None
    
  #def test_parse_and_run(self):
  #  assert 1 == 0
  #
  #def test_pairwise_comparison(self):
  #  assert 1 == 0
  #
  #def test_delete_files_based_on_list_of_regexes(self):
  #  assert 1 == 0
  #
  #def test_use_bundled_exec(self):
  # this wont work properly if the python script is installed.
  #  assert 1 == 0
        

if __name__ == "__main__":
  unittest.main()

