#!/usr/bin/env python
# File created on 02 Jul 2014
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from glob import glob
from os.path import join

from qiime.util import parse_command_line_parameters, make_option, create_dir
from cogent.parse.fasta import MinimalFastaParser

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
# Members of the tuple in script_usage are (title, description, example call)
script_info['script_usage'] = [("","Join all files whose names start with 'in' and end with 'fasta' in a new fasta file.",'%prog -i "in*fasta" -o ./out/')]
script_info['output_description']= ""
script_info['required_options'] = [
    make_option('-i', '--input_glob',
                help='pattern to match input filenames'),
    make_option('-o', '--output_dir', type="new_dirpath",
                help='the output directory')
]
script_info['optional_options'] = []
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    create_dir(opts.output_dir)
    input_fps = glob(opts.input_glob)
    output_fasta_f = open(join(opts.output_dir,'reference.fasta'), 'w')
    output_map_f = open(join(opts.output_dir,'reference_map.tsv'), 'w')
    i = 1
    for fp in input_fps:
        for seq_id, seq in MinimalFastaParser(open(fp,'U')):
            seq_id_fields = seq_id.strip().split()
            output_map_f.write('\t'.join([str(i), seq_id_fields[0]]))
            output_map_f.write('\n')
            output_fasta_f.write('>%s %s\n%s\n' % (str(i), seq_id, seq))
            i += 1
    output_fasta_f.close()
    output_map_f.close()




if __name__ == "__main__":
    main()
