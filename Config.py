#!/usr/bin/env python

import os
import sys
import re

'''
   In the event that a user does not have a .cafarc file in his main
   directory, this script will go ahead and create one for the user
   and then proceed with running the benchmark program.

   The .cafarc file basically details path and directory that will be used 
   throughout the program.

   If there is a change in the evidence codes being used or if uniprot-goa
   changes its ftp path, changes need to be made in the .cafarc file and
   it should get reflected throughout the program.
'''

def create():
    outfile_handle = open('.cafarc', 'w')
    outfile_handle.write('[WORKDIR]\n')
    
    work_dir_response = raw_input('Provide a path to your working directory (If left blank, defaults  to current directory) : ')
    if work_dir_response == '':
        outfile_handle.write('DEFAULT_PATH : .\n') 
    elif work_dir_response.startswith('.') or work_dir_response.startswith('/'):
        outfile_handle.write('DEFAULT_PATH : ' + work_dir_response + '\n')
    else:
        outfile_handle.write('DEFAULT_PATH : ' + './' + work_dir_response + '\n')
        
    outfile_handle.write('\n')

    outfile_handle.write('[FTP]\n')
    outfile_handle.write('HOSTNAME : ftp.ebi.ac.uk\n')
    outfile_handle.write('CURRENT_FILE_PATH : /pub/databases/GO/goa/UNIPROT\n')
    outfile_handle.write('OLD_FILE_PATH : /pub/databases/GO/goa/old/UNIPROT\n')
    outfile_handle.write('\n')
    
    outfile_handle.write('[DEFAULTS]\n')
    outfile_handle.write('EXP_EVIDENCE_CODES : ' + str(set(['EXP','IDA','IPI','IMP','IGI','IEP'])) + '\n')
    outfile_handle.write('ONTOLOGIES : ' + str(set(['F','P','C'])) + '\n')
    outfile_handle.write('TAXONOMY_FILENAME : names.dmp\n')
    
    outfile_handle.write('\n')

    outfile_handle.write('[SEQUENCE]\n')
    outfile_handle.write('BASE_URL : www.uniprot.org/uniprot/\n')

    outfile_handle.write('\n')

    outfile_handle.write('[REGEX]\n')
    outfile_handle.write('FTP_DATE : [a-zA-Z]+\_\d+\n')
    outfile_handle.write('FTP_FILE_START : gene_association\n')

if __name__ == '__main__':
    create()