#!/usr/bin/env python

import unittest
import os
import shutil
import sys
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

class Tests(unittest.TestCase):
    def run_modeller_script(self, script_name, model_name, resrng):
        """Run a Modeller script and test the output model"""
        os.chdir(os.path.join(TOPDIR, 'comparative_modeling'))
        # Run script
        p = subprocess.check_call(["python", script_name, "--test"])
        # Make sure PDB was produced with the requested residue range
        with open('%s.B99990001.pdb' % model_name) as fh:
            pdb_lines = [x for x in fh.readlines() if x.startswith('ATOM')]
        rng = (int(pdb_lines[0][22:26]), int(pdb_lines[-1][22:26]))
        self.assertEqual(rng, resrng)

    def test_model_352_504(self):
        """Test generation of comparative model for region 352-504"""
        self.run_modeller_script('model_ecm29_352.py',
                                 'ecm29_352_504', (352, 504))

    def test_model_686_1738(self):
        """Test generation of comparative model for region 686-1738"""
        self.run_modeller_script('model_ecm29_686.py',
                                 'ecm29_684_1749', (686, 1738))

    def test_simple(self):
        """Test model building"""
        os.chdir(os.path.join(TOPDIR, 'production_scripts'))
        p = subprocess.check_call(["python", "smodeling.py", "--test"])
        # todo: assert outputs, run analysis

if __name__ == '__main__':
    unittest.main()
