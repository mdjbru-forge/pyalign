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

### ** Test files

testDir = "tests/test_ungapFastaFile_files"
def tf(f) :
    # Convert a file name to a path to this file in the test directory
    return os.path.join(testDir, f)
def rmTf() :
    # Remove the output files produced during test (start with "_")
    [os.remove(tf(x)) for x in os.listdir(testDir) if x.startswith("_")]
    

### * Tests

class TestUngapFastaFile(unittest.TestCase) :

    def assertFileEqual(self, f1, f2) :
        with open(f1, "r") as fi1 :
            with open(f2, "r") as fi2 :
                self.assertEqual(fi1.read(), fi2.read())

    def tearDown(self) :
        rmTf()
                
    def test_noGap_000(self) :
        mod.ungapFastaFile(tf("inputFile1.fa"), tf("_testFile1.fa"))
        self.assertFileEqual(tf("_testFile1.fa"), tf("outputFile1.fa"))

    def test_gapInNames_000(self) :
        mod.ungapFastaFile(tf("inputFile2.fa"), tf("_testFile2.fa"))
        self.assertFileEqual(tf("_testFile2.fa"), tf("outputFile2.fa"))

    def test_gapInSeqs_000(self) :
        mod.ungapFastaFile(tf("inputFile3.fa"), tf("_testFile3.fa"))
        self.assertFileEqual(tf("_testFile3.fa"), tf("outputFile3.fa"))
        
