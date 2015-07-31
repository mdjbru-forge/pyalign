### * Description

### ** Requirements

# sudo pip install nose
# sudo pip install coverage

### ** Usage

# nosetests ./
# nosetests ./ --with-coverage --cover-html

### * Setup

### ** Import

import unittest
import sys
import os
sys.path.insert(0, os.path.abspath(".."))

import pyalign as mod

### * Tests

class TestRunMafft(unittest.TestCase) :

    def setUp(self) :
        self.oldOsSystem = mod.os.system
        def f(command) :
            self.command = command
        mod.os.system = f

    def tearDown(self) :
        mod.os.system = self.oldOsSystem
    
    def test_run_000(self) :
        mod.runMafft("inputToto", "outputToto", 4)
        result = self.command
        expected = "mafft --thread 4 inputToto | awk '{if (substr($0,1,1)==\">\"){if (p){printf \"\\n\";} print $0} else printf(\"%s\",$0);p++;} END {print \"\\n\"}' > outputToto"
        self.assertEqual(result, expected)
