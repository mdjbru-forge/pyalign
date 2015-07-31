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

class TestRandomTag(unittest.TestCase) :

    def setUp(self) :
        self.allowedChar = "0123456789abcdef"
    
    def test_tagLength_000(self) :
        result = mod.randomTag(6)
        self.assertEqual(len(result), 6)

    def test_tagLength_001(self) :
        result = mod.randomTag(10)
        self.assertEqual(len(result), 10)

    def test_tagLength_zero(self) :
        result = mod.randomTag(0)
        self.assertEqual(result, "")

    def test_tagComposition_000(self) :
        result = mod.randomTag(10)
        check = [x in self.allowedChar for x in result]
        self.assertTrue(all(check))

