#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_cdenrichrgenestoterm
----------------------------------

Tests for `cdenrichrgenestoterm` module.
"""

import os
import sys
import unittest
import tempfile
import shutil
from unittest.mock import MagicMock
import pandas as pd


from cdenrichrgenestoterm import cdenrichrgenestoterm


class TestCdenrichrgenestoterm(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def get_default_genesets(self):
        return 'GO_Biological_Process_2018,' \
               'GO_Cellular_Component_2018,' \
               'GO_Molecular_Function_2018,' \
               'KEGG_2019_Human,Reactome_2016,' \
               'WikiPathways_2019_Human,' \
               'Human_Phenotype_Ontology,Jensen_DISEASES'

    def test_read_inputfile(self):
        temp_dir = tempfile.mkdtemp()
        try:
            tfile = os.path.join(temp_dir, 'foo')
            with open(tfile, 'w') as f:
                f.write('hellothere')
            res = cdenrichrgenestoterm.read_inputfile(tfile)
            self.assertEqual('hellothere', res)
        finally:
            shutil.rmtree(temp_dir)

    def test_parse_args(self):
        myargs = ['inputarg']
        res = cdenrichrgenestoterm._parse_arguments('desc',
                                                      myargs)
        self.assertEqual('inputarg', res.input)
        self.assertEqual(0.05, res.cutoff)
        self.assertEqual('/tmp', res.tmpdir)
        self.assertEqual('GO_Biological_Process_2018,' +
                         'GO_Cellular_Component_2018,' +
                         'GO_Molecular_Function_2018,' +
                         'KEGG_2019_Human,Reactome_2016,' +
                         'WikiPathways_2019_Human,' +
                         'Human_Phenotype_Ontology,' +
                         'Jensen_DISEASES', res.genesets)

    def test_run_gprofiler_no_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            tfile = os.path.join(temp_dir, 'foo')
            myargs = [tfile]
            theargs = cdenrichrgenestoterm._parse_arguments('desc',
                                                              myargs)
            try:
                cdenrichrgenestoterm.run_enrichr(tfile,
                                                 theargs)
                self.fail('Expected FileNotFoundError')
            except FileNotFoundError:
                pass
        finally:
            shutil.rmtree(temp_dir)

    def test_run_gprofiler_empty_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            tfile = os.path.join(temp_dir, 'foo')
            open(tfile, 'a').close()
            myargs = [tfile]
            theargs = cdenrichrgenestoterm._parse_arguments('desc',
                                                              myargs)
            res = cdenrichrgenestoterm.run_enrichr(tfile,
                                                   theargs)
            self.assertEqual(None, res)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_gprofiler_with_empty_result(self):
        temp_dir = tempfile.mkdtemp()
        try:
            enrichr = MagicMock()
            mockres = MagicMock()
            mockres.res2d = pd.DataFrame()
            enrichr.enrichr = MagicMock(return_value=mockres)
            tfile = os.path.join(temp_dir, 'foo')
            with open(tfile, 'w') as f:
                f.write('a,b,c')
            myargs = [tfile, '--tmpdir', temp_dir]
            theargs = cdenrichrgenestoterm._parse_arguments('desc',
                                                            myargs)
            res = cdenrichrgenestoterm.run_enrichr(tfile,
                                                   theargs,
                                                   enrichr=enrichr)
            gs = self.get_default_genesets()
            self.assertEqual(None, res)
            enrichr.enrichr.assert_called_once_with(gene_list=['A', 'B', 'C'],
                                                    cutoff=0.05,
                                                    gene_sets=gs,
                                                    no_plot=True,
                                                    outdir=temp_dir)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_gprofiler_with_valid_result(self):
        temp_dir = tempfile.mkdtemp()
        try:
            enrichr = MagicMock()
            mockres = MagicMock()
            df = pd.DataFrame(columns=['Term',
                                       'Gene_set',
                                       'Adjusted P-value',
                                       'Genes', 'Overlap'],
                              data=[['term1', 'set1',
                                     0.5,
                                     'B;C', '2/25'],
                                    ['term2', 'set2',
                                     0.3, 'A;C', '7/9']])
            mockres.res2d = df
            enrichr.enrichr = MagicMock(return_value=mockres)
            tfile = os.path.join(temp_dir, 'foo')
            with open(tfile, 'w') as f:
                f.write('a,b,c')
            myargs = [tfile, '--tmpdir', temp_dir]
            theargs = cdenrichrgenestoterm._parse_arguments('desc',
                                                            myargs)
            res = cdenrichrgenestoterm.run_enrichr(tfile,
                                                   theargs,
                                                   enrichr=enrichr)
            gs = self.get_default_genesets()
            self.assertEqual('term2', res['name'])
            self.assertEqual('set2', res['source'])
            self.assertEqual(0.3, res['p_value'])
            self.assertEqual('', res['description'])
            self.assertEqual(['A', 'C'], res['intersections'])
            self.assertEqual(9, res['term_size'])
            enrichr.enrichr.assert_called_once_with(gene_list=['A', 'B', 'C'],
                                                    cutoff=0.05,
                                                    gene_sets=gs,
                                                    no_plot=True,
                                                    outdir=temp_dir)
        finally:
            shutil.rmtree(temp_dir)

    def test_main_invalid_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            tfile = os.path.join(temp_dir, 'foo')
            myargs = ['prog', tfile]
            res = cdenrichrgenestoterm.main(myargs)
            self.assertEqual(2, res)
        finally:
            shutil.rmtree(temp_dir)

    def test_main_empty_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            tfile = os.path.join(temp_dir, 'foo')
            open(tfile, 'a').close()
            myargs = ['prog', tfile]
            res = cdenrichrgenestoterm.main(myargs)
            self.assertEqual(0, res)
        finally:
            shutil.rmtree(temp_dir)


if __name__ == '__main__':
    sys.exit(unittest.main())
