#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_integration_of_enrichr
----------------------------------

Integration tests for `enrichr` module.
"""

import os
import sys
import unittest
import tempfile
import shutil
import json
import subprocess
import stat
from unittest.mock import MagicMock
from cdenrichrgenestoterm import cdenrichrgenestoterm


SKIP_REASON = 'CDENRICHR_DOCKER_IMAGE, CDENRICHR_DOCKER, CDENRICHR_TMPDIR environment ' \
              'variable(s) not set to a docker image' \
              ' cannot run integration tests of Enrichr with Docker'


@unittest.skipUnless(os.getenv('CDENRICHR_DOCKER_IMAGE') is not None and
                     os.getenv('CDENRICHR_DOCKER') is not None and
                     os.getenv('CDENRICHR_TMPDIR') is not None, SKIP_REASON)
class TestCdhidefInDocker(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def run_enrichr_docker(self, cmdargs, temp_dir=None):
        """
        Runs hidef command as a command line process
        :param cmd_to_run: command to run as list
        :type cmd_to_run: list
        :return: (return code, standard out, standard error)
        :rtype: tuple
        """
        cmd = [os.getenv('CDENRICHR_DOCKER'), 'run', '--rm', '-v',
               temp_dir+':'+temp_dir,
               os.getenv('CDENRICHR_DOCKER_IMAGE')]
        cmd.extend(cmdargs)
        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)

        out, err = p.communicate()

        return p.returncode, out, err

    def test_run_enrichr_small_gene_list_may_take_one_minute(self):
        temp_dir = tempfile.mkdtemp(dir=os.getenv('CDENRICHR_TMPDIR'))
        try:
            gene_list = 'NR1H3,MNDA,APOE,NCF2,CRIP3,NCF1,NCF4,CD48,LILRA6,' \
                        'LILRA4,CD44,LILRA5,CD47,LILRA2,LILRA1,ADA2,CD40,' \
                        'LINGO3,RAB20,CD14,MPP1,CD19,PTPN7,PTPN6,FTHL17,' \
                        'LHFPL2,UNC13D,COL23A1'
            input_file = os.path.join(temp_dir, 'input.txt')
            with open(input_file, 'w') as f:
                f.write(gene_list + '\n')
            ecode, out, err = self.run_enrichr_docker([input_file],
                                                      temp_dir=temp_dir)
            self.assertEqual(0, ecode)

            res = json.loads(out)
            self.assertEqual(6, len(res.keys()))
            self.assertEqual('Chronic granulomatous disease', res['name'])
            self.assertEqual('Jensen_DISEASES', res['source'])
            self.assertEqual('', res['description'])
            self.assertEqual(5, res['term_size'])
            self.assertEqual(3, len(res['intersections']))

        finally:
            shutil.rmtree(temp_dir)

    def test_run_enrichr_two_gene_items_may_take_one_minute(self):
        temp_dir = tempfile.mkdtemp(dir=os.getenv('CDENRICHR_TMPDIR'))
        try:
            gene_list = 'MTOR,TP53'
            input_file = os.path.join(temp_dir, 'input.txt')
            with open(input_file, 'w') as f:
                f.write(gene_list + '\n')
            ecode, out, err = self.run_enrichr_docker([input_file,
                                                       '--genesets',
                                                       'Jensen_DISEASES'],
                                                      temp_dir=temp_dir)
            self.assertEqual(0, ecode)

            res = json.loads(out)
            print(out)
            self.assertEqual(6, len(res.keys()))
            self.assertEqual('Neurofibromatosis', res['name'])
            self.assertEqual('Jensen_DISEASES', res['source'])
            self.assertEqual('', res['description'])
            self.assertEqual(38, res['term_size'])
            self.assertEqual(2, len(res['intersections']))

        finally:
            shutil.rmtree(temp_dir)


if __name__ == '__main__':
    sys.exit(unittest.main())
