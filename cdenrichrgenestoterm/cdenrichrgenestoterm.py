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
    parser.add_argument('--maxpval', type=float, default=0.05,
                        help='Max p value')
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
                enrichr=gseapy,
                retry_count=2):
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
    cur_try = 1
    with redirect_stdout(sys.stderr):
        while cur_try <= retry_count:
            try:
                res = enrichr.enrichr(gene_list=genes, gene_sets=theargs.genesets,
                                      cutoff=theargs.maxpval,
                                      no_plot=True, outdir=theargs.tmpdir)
                break
            except Exception as e:
                sys.stderr.write('Try # ' + str(cur_try) + ' caught exception: ' + str(e))
                cur_try += 1
                if cur_try > retry_count:
                    sys.stderr.write('Retries exceeded')
                    return None
        sys.stderr.write('res object: ' + str(res) + '\n')
        sys.stderr.write('res.res2d object: ' + str(res.res2d) + '\n')

    df_result = res.res2d
    if df_result.shape[0] == 0:
        sys.stderr.write('Empty data frame\n')
        return None
    """
    Example output:
    >>> df_result
          Gene_set                      Term  Overlap   P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  Combined Score Genes
    0  Jensen_DISEASES                 Keratosis     1/27  0.460766               1.0            0                     0    1.638807    1.269854e+00   TAT
    1  Jensen_DISEASES  Hepatocellular carcinoma     1/29  0.484897               1.0            0                     0    1.525786    1.104393e+00   TAT
    2  Jensen_DISEASES    Palmoplantar keratosis     1/31  0.507950               1.0            0                     0    1.427348    9.668460e-01   TAT
    3  Jensen_DISEASES                 Keratitis     1/39  0.590316               1.0            0                     0    1.134559    5.980234e-01   TAT
    4  Jensen_DISEASES   Intellectual disability    1/296  0.998903               1.0            0                     0    0.149486    1.641426e-04   TAT
    5  Jensen_DISEASES                 Carcinoma  1/11318  0.999993               1.0            0                     0    0.003910    2.732093e-08   TAT
    """
    df_result.rename(columns={'Adjusted P-value': 'apv'}, inplace=True)
    # filter out any rows where min overlap is not met
    df_result.drop(df_result[df_result.apv > theargs.maxpval].index,
                   inplace=True)
    if df_result.shape[0] == 0:
        sys.stderr.write('Empty data frame after p value filter\n')
        return None
    df_result.sort_values('apv',
                          ascending=True, inplace=True)
    df_result.reset_index(drop=True, inplace=True)
    theres = {'name': df_result['Term'][0],
              'source': df_result['Gene_set'][0],
              'p_value': df_result['apv'][0],
              'description': '',
              'term_size': int(df_result['Overlap'][0][df_result['Overlap'][0].index('/')+1:]),
              'intersections': df_result['Genes'][0].split(';')}
    sys.stderr.write('About to return this fragment: ' + str(theres) + '\n')
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
         "sourceTermId": "IS THE ID FOR THE ENRICHED TERM/FUNCTIONAL CATEGORY IN ITS NATIVE NAMESPACE"
         "p_value": Adjusted P-value from Enrichr,
         "description": "EMPTY STRING",
         "intersections": "List of Genes that intersect"
        }
        
    """

    theargs = _parse_arguments(desc, args[1:])

    try:
        inputfile = os.path.abspath(theargs.input)
        theres = run_enrichr(inputfile, theargs)
        sys.stderr.flush()
        if theres is None:
            sys.stderr.write('No terms found\n')
        else:
            json.dump(theres, sys.stdout)
        sys.stdout.flush()
        return 0
    except Exception as e:
        sys.stderr.write('Caught exception: ' + str(e))
        return 2
    finally:
        sys.stderr.flush()


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
