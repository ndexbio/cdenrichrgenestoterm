#!/usr/bin/env python

import os
import sys
import argparse
import json
from contextlib import redirect_stdout

with redirect_stdout(sys.stderr):
    import gseapy


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


def _parse_arguments(desc, args):
    """
    Parses command line arguments
    :param desc:
    :param args:
    :return:
    """
    help_fm = Formatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_fm)
    parser.add_argument('input',
                        help='comma delimited list of genes in file')
    parser.add_argument('--cutoff', type=float, default=0.05,
                        help='Cutoff value')
    parser.add_argument('--tmpdir', default='/tmp',
                        help='Temp directory to hold output from task')
    parser.add_argument('--genesets', default='GO_Biological_Process_2018,'
                                              'GO_Cellular_Component_2018,'
                                              'GO_Molecular_Function_2018,'
                                              'KEGG_2019_Human,Reactome_2016,'
                                              'WikiPathways_2019_Human,'
                                              'Human_Phenotype_Ontology,'
                                              'Jensen_DISEASES',
                        help='Gene sets to enrich against. '
                             'Should be comma delimited')
    return parser.parse_args(args)


def read_inputfile(inputfile):
    """

    :param inputfile:
    :return:
    """
    with open(inputfile, 'r') as f:
        return f.read()


def run_enrichr(inputfile, theargs,
                enrichr=gseapy):
    """
    todo
    :param inputfile:
    :return:
    """
    genes = read_inputfile(inputfile)
    genes = genes.strip(',').strip('\n').upper().split(',')
    if genes is None or (len(genes) == 1 and len(genes[0].strip()) == 0):
        sys.stderr.write('No genes found in input')
        return None
    with redirect_stdout(sys.stderr):
        res = enrichr.enrichr(gene_list=genes, gene_sets=theargs.genesets,
                              cutoff=theargs.cutoff,
                              no_plot=True, outdir=theargs.tmpdir)
    df_result = res.res2d
    if df_result.shape[0] == 0:
        return None
    df_result.sort_values('Adjusted P-value',
                          ascending=True, inplace=True)
    df_result.reset_index(drop=True, inplace=True)
    theres = {'name': df_result['Term'][0],
              'source': df_result['Gene_set'][0],
              'p_value': df_result['Adjusted P-value'][0],
              'description': '',
              'term_size': int(df_result['Overlap'][0][df_result['Overlap'][0].index('/')+1:]),
              'intersections': df_result['Genes'][0].split(';')}

    return theres


def main(args):
    """
    Main entry point for program

    :param args: command line arguments usually :py:const:`sys.argv`
    :return: 0 for success otherwise failure
    :rtype: int
    """
    desc = """
        Running Enrichr via gseapy .

        Takes file with comma delimited list of genes as input and
        outputs best matching term (as determined by Adjusted P value)
        if any in json format:
        
        {
         "name": "TERM",
         "source": "SOURCE OF TERM",
         "p_value": Adjusted P-value from Enrichr,
         "description": "EMPTY STRING",
         "intersections": "List of Genes that intersect"
        }
        
    """

    theargs = _parse_arguments(desc, args[1:])

    try:
        inputfile = os.path.abspath(theargs.input)
        theres = run_enrichr(inputfile, theargs)
        if theres is None:
            sys.stderr.write('No terms found\n')
        else:
            json.dump(theres, sys.stdout)
        sys.stdout.flush()
        return 0
    except Exception as e:
        sys.stderr.write('Caught exception: ' + str(e))
        return 2


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
