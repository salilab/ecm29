#!/usr/bin/env python

import unittest
import os
import shutil
import sys
import subprocess
import ihm.reader

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

    def test_mmcif(self):
        """Test generation of mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'production_scripts'))
        if os.path.exists("ecm29.cif"):
            os.unlink("ecm29.cif")
        # Potentially override methods that need network access
        env = os.environ.copy()
        env['PYTHONPATH'] = os.path.join(TOPDIR, 'test', 'mock') \
                            + ':' + env.get('PYTHONPATH', '')
        p = subprocess.check_call(
                ["python", "smodeling.py", "--mmcif", "--dry-run"], env=env)
        # Check output file
        self._check_mmcif_file('ecm29.cif')

    def _check_mmcif_file(self, fname):
        with open(fname) as fh:
            s, = ihm.reader.read(fh)
        self.assertEqual(len(s.citations), 1)
        self.assertEqual(s.citations[0].doi, '10.1074/jbc.M117.803619')
        self.assertEqual(len(s.software), 3)
        self.assertEqual(len(s.orphan_starting_models), 20)
        # Should be 1 state
        self.assertEqual(len(s.state_groups), 1)
        state1, = s.state_groups[0]
        # Should be 2 models
        self.assertEqual(sum(len(x) for x in state1), 2)
        # Check # of spheres and atoms in each model
        m1 = state1[0][0]
        m2 = state1[1][0]
        for m in m1, m2:
            self.assertEqual(len(m._spheres), 8142)
            self.assertEqual(len(m._atoms), 0)
        # Should be 2 ensembles
        self.assertEqual([e.num_models for e in s.ensembles], [11980, 6261])
        # Just one restraint - crosslinks
        xl, = s.restraints
        self.assertEqual(len(xl.experimental_cross_links), 63)
        self.assertEqual(len(xl.cross_links), 63)


if __name__ == '__main__':
    unittest.main()
