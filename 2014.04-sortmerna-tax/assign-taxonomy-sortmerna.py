#!/usr/bin/env python
# File created on 11 Apr 2014
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from tempfile import mkdtemp
from collections import defaultdict, Counter
from os.path import join
from shutil import rmtree

from skbio.parse.sequences import parse_fasta
from qiime.util import (parse_command_line_parameters, make_option,
                        load_qiime_config, get_qiime_temp_dir,
                        qiime_system_call, create_dir, remove_files)
from qiime.assign_taxonomy import TaxonAssigner

qiime_config = load_qiime_config()

default_reference_seqs_fp = qiime_config['assign_taxonomy_reference_seqs_fp']
default_id_to_taxonomy_fp = qiime_config['assign_taxonomy_id_to_taxonomy_fp']

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
# Members of the tuple in script_usage are (title, description, example call)
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [
    make_option('-i', '--input_query_fp', type='existing_filepath',
                help='the query sequences'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='directory to store output and log files')
]
script_info['optional_options'] = [
     make_option('-t', '--id_to_taxonomy_fp', type="existing_filepath",
                 help='Path to tab-delimited file mapping sequences to their '
                 'taxonomy. Each assigned taxonomy is provided as a '
                 'semicolon-separated string. [default: %s]'
                 % default_id_to_taxonomy_fp,
                 default=default_id_to_taxonomy_fp),
     make_option('-r', '--reference_seqs_fp', type="existing_filepath",
                 help='Path to reference sequences.  These '
                 'are indexed to create the reference database. '
                 '[default: %s; REQUIRED if -b is not provided]'
                 % default_reference_seqs_fp,
                 default=default_reference_seqs_fp),
     make_option('-b', '--reference_seqs_idx_fp',
                 help='Path to pre-indexed reference sequences.  These '
                 'are indexed to create the reference database. '
                 '[default: computed on-the-fly]', default=None),
     make_option('--sortmerna_params', default=None,
                 help='string of parameters to pass to sortmerna'
                 ' [default: no additional parameters as passed]'),
     make_option('--indexdb_params', default=None,
                 help='string of parameters to pass to indexdb_rna'
                 ' [default: no additional parameters as passed]'),
     make_option('--min_consensus_fraction', default=0.51, type=float,
                 help='Minimum fraction of database hits that must have a '
                 'specific taxonomic assignment to assign that taxonomy '
                 'to a query [default: %default]'),
     make_option('--min_percent_id', default=0.0, type=float,
                 help='Minimum percent identity to consider an alignment '
                 'to be a hit [default: all alignments are considered hits]'),
     make_option('--best', default=None, type=int,
                 help='sortmerna\'s --best parameter [default: %default]'),
     make_option('--num_alignments', default=None, type=int,
                 help='sortmerna\'s --num_alignments parameter [default: %default]'),
     make_option('--min_lis', default=2, type=int,
                 help='sortmerna\'s min_lis parameter [default: %default]'),
]
script_info['version'] = __version__

def call_cmd(cmd, HALT_EXEC):
    if HALT_EXEC:
        print cmd
        exit(0)
    else:
        stdout, stderr, exit_status = qiime_system_call(cmd)
        if exit_status != 0:
            print "indexdb_rna failed!\nSTDOUT\n%s\nSTDERR\n%s\n" \
                   % (stdout, stderr)
            exit(1)
    return cmd

def sortmerna_indexdb(input_fp, output_fp, params="", HALT_EXEC=False):
    """
    """

    cmd = "indexdb_rna --ref %s,%s -v %s" % (input_fp, output_fp, params)
    return call_cmd(cmd, HALT_EXEC)


def sortmerna_map(query_fp, refseqs_fasta_fp, refseqs_index_fp, blast_output_fp,
                  min_lis=2, best=None, num_alignments=None, params="",
                  HALT_EXEC=False):
    """
    """
    if best is not None:
        params = " ".join([params, "--best", str(best), "--min_lis",
                           str(min_lis)])
    elif num_alignments is not None:
        params = " ".join([params, "--num_alignments", str(num_alignments)])

    cmd = "sortmerna --ref %s,%s --reads %s --blast 1 --aligned %s -v %s" \
          % (refseqs_fasta_fp, refseqs_index_fp, query_fp, blast_output_fp,
             params)
    return call_cmd(cmd, HALT_EXEC)

def blast_to_tax_assignments(blast_output_f, id_to_taxonomy_map, min_percent_id=0.0):
    """
    """

    result = defaultdict(list)
    for line in blast_output_f:
        fields = line.split('\t')
        subject_id = fields[1]
        percent_id = float(fields[2])
        if percent_id > min_percent_id:
            subject_tax = id_to_taxonomy_map[subject_id]
            result[fields[0]].append([e.strip() for e in subject_tax.split(';')])
    return result

def tax_assignments_to_consensus_assignments(query_to_assignments,
                                             min_consensus_fraction=0.51):
    """
    """
    for query_id, assignments in query_to_assignments.iteritems():
        # this call will get cleaned up
        consensus_assignment = \
            get_consensus_assignment(assignments,
            min_consensus_fraction=min_consensus_fraction)
        query_to_assignments[query_id] = consensus_assignment

    return query_to_assignments

def get_consensus_assignment(assignments, unassignable_label='unassigned',
                             min_consensus_fraction=0.51):
    """ compute the consensus assignment from a list of assignments

        This code was pulled almost exactly from QIIME's
         UclustConsensusTaxonAssigner._get_consensus_assignment method.
    """
    num_input_assignments = len(assignments)
    consensus_assignment = []

    # if the assignments don't all have the same number
    # of levels, the resulting assignment will have a max number
    # of levels equal to the number of levels in the assignment
    # with the fewest number of levels. this is to avoid
    # a case where, for example, there are n assignments, one of
    # which has 7 levels, and the other n-1 assignments have 6 levels.
    # A 7th level in the result would be misleading because it
    # would appear to the user as though it was the consensus
    # across all n assignments.
    num_levels = min([len(a) for a in assignments])

    # iterate over the assignment levels
    for level in range(num_levels):
        # count the different taxonomic assignments at the current level.
        # the counts are computed based on the current level and all higher
        # levels to reflect that, for example, 'p__A; c__B; o__C' and
        # 'p__X; c__Y; o__C' represent different taxa at the o__ level (since
        # they are different at the p__ and c__ levels).
        current_level_assignments = \
            Counter([tuple(e[:level + 1]) for e in assignments])
        # identify the most common taxonomic assignment, and compute the
        # fraction of assignments that contained it. it's safe to compute the
        # fraction using num_assignments because the deepest level we'll
        # ever look at here is num_levels (see above comment on how that
        # is decided).
        tax, max_count = current_level_assignments.most_common(1)[0]
        max_consensus_fraction = max_count / num_input_assignments
        # check whether the most common taxonomic assignment is observed
        # in at least min_consensus_fraction of the sequences
        if max_consensus_fraction >= min_consensus_fraction:
            # if so, append the current level only (e.g., 'o__C' if tax is
            # 'p__A; c__B; o__C', and continue on to the next level
            consensus_assignment.append((tax[-1], max_consensus_fraction))
        else:
            # if not, there is no assignment at this level, and we're
            # done iterating over levels
            break

    # construct the results
    # determine the number of levels in the consensus assignment
    consensus_assignment_depth = len(consensus_assignment)
    if consensus_assignment_depth > 0:
        # if it's greater than 0, generate a list of the
        # taxa assignments at each level
        assignment_result = [a[0] for a in consensus_assignment]
        # and assign the consensus_fraction_result as the
        # consensus fraction at the deepest level
        consensus_fraction_result = \
            consensus_assignment[consensus_assignment_depth - 1][1]
    else:
        # if there are zero assignments, indicate that the taxa is
        # unknown
        assignment_result = [unassignable_label]
        # and assign the consensus_fraction_result to 1.0 (this is
        # somewhat arbitrary, but could be interpreted as all of the
        # assignments suggest an unknown taxonomy)
        consensus_fraction_result = 1.0

    return (
        assignment_result, consensus_fraction_result, num_input_assignments
    )

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    input_query_fp = opts.input_query_fp
    id_to_taxonomy_fp = opts.id_to_taxonomy_fp
    reference_seqs_fp = opts.reference_seqs_fp
    reference_seqs_idx_fp = opts.reference_seqs_idx_fp
    output_dir = opts.output_dir
    indexdb_params = opts.indexdb_params or ""
    sortmerna_params = opts.sortmerna_params or ""
    min_lis = opts.min_lis
    best = opts.best
    num_alignments = opts.num_alignments
    min_consensus_fraction = opts.min_consensus_fraction
    min_percent_id = opts.min_percent_id

    if (best is None) and (num_alignments is None):
        option_parser.error("Either --best or --num_alignments must be passed.")
    elif (best is not None) and (num_alignments is not None):
        option_parser.error("--best and --num_alignments cannot both be passed.")
    else:
        pass

    create_dir(output_dir)
    dirs_to_remove = []

    qiime_temp_dir = get_qiime_temp_dir()
    output_fp = join(output_dir, 'assignments.tsv')
    command_log_fp = join(output_dir, 'cmds.log')
    output_f = open(output_fp, 'w')
    command_log_f = open(command_log_fp, 'w')

    id_to_taxonomy_f = open(id_to_taxonomy_fp, 'U')
    id_to_taxonomy_map = \
        TaxonAssigner._parse_id_to_taxonomy_file(id_to_taxonomy_f)
    id_to_taxonomy_f.close()

    if reference_seqs_idx_fp is None:
        index_dir = mkdtemp(dir=qiime_temp_dir)
        create_dir(index_dir)
        dirs_to_remove.append(index_dir)
        reference_seqs_idx_fp = join(index_dir, 'ref.idx')
        cmd = sortmerna_indexdb(reference_seqs_fp, reference_seqs_idx_fp,
                                params=indexdb_params)
        command_log_f.write(cmd)
        command_log_f.write('\n')

    blast_output_basename = join(output_dir, 'log')
    blast_output_fp = '.'.join([blast_output_basename, 'blast'])
    cmd = sortmerna_map(input_query_fp, reference_seqs_fp,
                        reference_seqs_idx_fp, blast_output_basename,
                        min_lis=min_lis, best=best, num_alignments=num_alignments,
                        params=sortmerna_params)
    command_log_f.write(cmd)
    command_log_f.write('\n')

    query_to_assignments = blast_to_tax_assignments(
        open(blast_output_fp), id_to_taxonomy_map, min_percent_id=min_percent_id)
    results = tax_assignments_to_consensus_assignments(
        query_to_assignments, min_consensus_fraction=min_consensus_fraction)
    output_f.write('#OTU ID\ttaxonomy\tconfidence\tnum hits\n')
    # this is ugly... we need results for all input sequences, but if there
    # are no hits from sortmerna for a sequence we won't have it in results.
    # we'll want to find a nicer solution than having to iterate over the input
    # file. maybe sortmerna can create a failures file...
    for query_id, _ in parse_fasta(open(input_query_fp)):
        query_id = query_id.split()[0]
        assignment = results[query_id]
        if len(assignment) == 0:
            assignment = [['Unassigned'], 0.0, 0]
        output_f.write('\t'.join([query_id, '; '.join(assignment[0]),
                                  str(assignment[1]), str(assignment[2])]))
        output_f.write('\n')

    # clean up time...
    output_f.close()
    command_log_f.close()
    map(rmtree, dirs_to_remove)

if __name__ == "__main__":
    main()
