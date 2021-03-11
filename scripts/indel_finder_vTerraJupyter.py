#! /usr/bin/env python

# import modules


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO

import os
import shutil
import glob
import sys

import pandas as pd
import re

# from datetime import date
# import time

######################
##### FUNCTIONS ######
######################

###############################################################################################
def get_ref_seq_record_id(ref_genome_path):
    ref = SeqIO.read(ref_genome_path, 'fasta')
    return {'ref_id': ref.id, 'ref_seq': ref.seq}

def create_temp_directories(wd):
    '''
    creates two temp directories for the alignment of each sequence to the reference genome 
    these temp directories will be deleted at the end
    '''
    # get date
    today_date = str(date.today())
    
    # create tempory directory to store fastas and alignments
    fasta_temp_dir = os.path.join(wd, 'temp-fastas_%s' % today_date)
    alignment_temp_dir = os.path.join(wd, 'temp-alignments_%s' % today_date)
    
    if os.path.exists(fasta_temp_dir):
        shutil.rmtree(fasta_temp_dir)
        
    if os.path.exists(alignment_temp_dir):
        shutil.rmtree(alignment_temp_dir)
        
    os.makedirs(fasta_temp_dir)
    os.makedirs(alignment_temp_dir)
    
    return {'alignment_temp_dir' : alignment_temp_dir, 'fasta_temp_dir': fasta_temp_dir}

def delete_temp_directory(directory_name):
    if os.path.exists(directory_name):
        print('removing the temp directory: %s' % directory_name)
        shutil.rmtree(directory_name)

def add_ref_genome(fasta_files_dir_path, wd, ref_genome_path, fasta_temp_dir) :
    '''
    generates a fasta file for each sequence that contains the reference sequence and sample sequence
    these are saved to the temp fasta files directory
    '''
    # read in the reference genome
    ref = SeqIO.read(ref_genome_path, 'fasta')

    os.chdir(fasta_files_dir_path)

    print('adding reference genome to each sequence and saving to temp fasta files directory')
    n = 0
    for file in glob.glob('*.fa*') :
        
        # reset_records list
        records = [ref]
        n = n +1
   
        # get record from fasta file
        record = SeqIO.read(file, 'fasta')
        print(n, record.id)

        #append record to reference genome
        records.append(record)

        #write records list to a fasta file
        path = os.path.join(fasta_temp_dir, '%s.fasta' % record.id)

        with open(path, 'w') as handle:
            SeqIO.write(records, handle, 'fasta')
            
def align_sequences(fasta_temp_dir, alignment_temp_dir, wd):      
    os.chdir(fasta_temp_dir)
    
    print('aligning each sample sequence to reference genome')
    n=0
    for file in glob.glob('*.fasta'):
        n = n + 1
        print(n)
        sample_seq_name = file.split('.fasta')[0]
#         for record in SeqIO.parse(file, 'fasta'):
#             if record.id != ref_id:
#                 sample_seq_name = record.id
                
      
        # create outpath file name for alignment
        alignment_file_name = os.path.join(alignment_temp_dir, '%s.alignment.fasta' % sample_seq_name)

        if not os.path.isfile(alignment_file_name):

            # do alignment
            mafft_cline = MafftCommandline (input=file)
            print(mafft_cline)
            stdout, stderr = mafft_cline()
            with open(alignment_file_name, 'w') as handle:
                handle.write(stdout)


def remove_insertions(alignment_temp_dir, fasta_temp_dir, ref_id, ref_seq, ref_genome_path ):
    print('recording insertions and removing insertions from sequences')
    
    ref = SeqIO.read(ref_genome_path, 'fasta')
    ref_seq = str(ref_seq)
    
    # prepare empty lists for data table
    seq_list = []
    start_list = []
    length_list = []
    ref_seq_list = []
    upstream_list = []
    downstream_list = []
    
    
    os.chdir(alignment_temp_dir)
    
    for file in glob.glob('*.alignment.fasta'):
        if not re.search('mod', file):
            alignment_name = file.split('.alignment')[0]
            for record in AlignIO.read(file, 'fasta'):
                seq_str = str(record.seq)

                if re.search('[-]+', seq_str) and record.id == ref_id:
                    insertions = re.findall('[-]+', seq_str)
                    print('')
                    print('found %d insertion(s) in %s' % (len(insertions), alignment_name))
                    print('will remove insertion(s) from %s and record in table' % alignment_name)
                 
                # we need the sequence of the sample:
                    for sample_record in AlignIO.read(file, 'fasta'):
                        if sample_record.id == alignment_name:
                            sample_seq_str= str(sample_record.seq)   
                    
                    moving_length = 0
                    k = 0 
                    for i in re.finditer('[-]+', seq_str):
                        k = k+1
                        print('remove insertion %d of %d' % (k , len(insertions)))
                                 
                        # locate first insert, record size and location, remove insert, 
                        # update the moving length repeat until no more insertions
                        seq_list.append(alignment_name)

                        start = i.start() - moving_length
                        length = i.end() - i.start()

                        start_list.append(start)
                        length_list.append(length)

                        # get the ref sequence
                        insert_nucleotides = ref_seq[start: start+length]
                        ref_seq_list.append('+%s' % insert_nucleotides)

                        # get the sequence around the insert
                        upstream = ref_seq[start-7: start]
                        downstream = ref_seq[start: start + 7]

                        upstream_list.append(upstream)
                        downstream_list.append(downstream)

                        # remove the insertion
                        first_half_seq =sample_seq_str[:start]
                        last_half_seq = sample_seq_str[start+length:]

                        new_seq = first_half_seq + last_half_seq
                        
                        # replaace the old sample_seq_str
                        sample_seq_str = new_seq
                        
                        # update the moving length
                        moving_length = moving_length + length 
                        
                        print('sample_seq is %d bases long' % len(sample_seq_str))
                        
                    # save the new sequence to temp records...
                    new_record = SeqRecord(Seq(sample_seq_str), id = alignment_name, name = alignment_name, description = '')
#                         records = [ref]
#                         records.append(new_record)
                    mod_fasta_path = os.path.join(alignment_temp_dir, '%s_mod.alignment.fasta' % alignment_name)
                    with open(mod_fasta_path, 'w') as handle:
                        SeqIO.write(new_record, handle, 'fasta')
                        


    # fill in pandas table
    df = pd.DataFrame()
    df['accession_id'] = seq_list
    df['indel'] = ref_seq_list
    df['start_pos'] = start_list
    df['length'] = length_list
    df['upstream_ref'] = upstream_list
    df['downstream_ref'] = downstream_list

    return {'insert_df': df, 'mod_seq_list': seq_list}

def record_deletions(alignment_temp_dir, ref_seq, insert_df, mod_seq_list, wd ) :
    print('recording deletions')
    
    # get date
    today_date = str(date.today())
    
    # prepare empty lists for data table
    ref_seq = str(ref_seq)
    seq_list = []
    start_list = []
    length_list = []
    ref_seq_list = []
    upstream_list = []
    downstream_list = []
    
    os.chdir(alignment_temp_dir)
     
    for file in glob.glob('*.alignment.fasta'):
        alignment_name = file.split('.alignment.fasta')[0]
        
        if alignment_name not in mod_seq_list:
            for record in AlignIO.read(file, 'fasta'):
                seq_str = str(record.seq)

                if re.search('[-]+', seq_str) and record.id == alignment_name:
                    deletions = re.findall('[-]+', seq_str)
                    print('found %d deletion(s) in %s' % (len(deletions), alignment_name))
                    print('will record deletions in table')
                    print('')

                    for i in re.finditer('[-]+', seq_str):

                            # locate first deletion, record size and location, repeat for next deletion
                            start = i.start()
                            if start == 0:
                                print('deletion at begining of sequence...will not record in table')
                            elif start != 0:
                                length = i.end() - i.start()

                                start_list.append(start)
                                length_list.append(length)
                                seq_list.append(alignment_name)

                                # get the ref sequence
                                insert_nucleotides = ref_seq[start: start+length]
                                ref_seq_list.append('-%s' % insert_nucleotides)

                                # get the sequence around the insert
                                upstream = ref_seq[start-7: start]
                                downstream = ref_seq[start: start + 7]

                                upstream_list.append(upstream)
                                downstream_list.append(downstream)

    df = pd.DataFrame()
    df['accession_id'] = seq_list
    df['indel'] = ref_seq_list
    df['start_pos'] = start_list
    df['length'] = length_list
    df['upstream_ref'] = upstream_list
    df['downstream_ref'] = downstream_list
    
    if insert_df.shape[0] > 0:
        dataframe_list = [insert_df, df]
        joined_df = pd.concat(dataframe_list)
    else:
        joined_df = df
    
    # write out joined dataframe
    path = os.path.join(wd, 'indel_results_table.csv')
    joined_df.to_csv(path, index = False)               
                
                
def create_dir_with_mod_fastas(alignment_temp_dir, mod_seq_list, wd):
    # get date
    today_date = str(date.today())
    
    if mod_seq_list:

        # create the directory
        mod_fasta_dir = os.path.join(wd, 'fasta_files_insertions_removed' )
        if os.path.exists(mod_fasta_dir):
            shutil.rmtree(mod_fasta_dir)
        os.mkdir(mod_fasta_dir)

        # create a copy of the mod fasta file in the directory
        for seq_name in mod_seq_list:
            file_name = os.path.join(alignment_temp_dir, '%s_mod.alignment.fasta' % seq_name)
            
            for record in AlignIO.read(file_name, 'fasta'):
                if record.id == seq_name:
                    
                    new_record = SeqRecord(record.seq, id = seq_name, description = '')
                    outpath = os.path.join(mod_fasta_dir, '%s.fasta' % seq_name)
                    
                    with open(outpath, 'w') as handle:
                        SeqIO.write(new_record, handle, 'fasta')
                
        # check length of all mod sequences
        print('checking length of all modified sequences')
        print('all sequences should be 29903 bp')

        os.chdir(mod_fasta_dir)
        for file in glob.glob('*.fasta'):
            record = SeqIO.read(file, 'fasta')
            print(record.id, 'length = %d' % len(record.seq))
            
    else:
        mod_fasta_dir = os.path.join(wd, 'fasta_files_insertions_removed' )
        print('no sequences with insertions')
        print('no fasta_files_with_insertions_removed directory is made')
      
    return mod_fasta_dir

def create_concatenated_seq_records(wd, temp_alignment_dir, mod_fasta_dir, mod_seq_list, ref_path):
    print('concatenating sequence records into single fasta alignment file')
    print('')
    
    # get date
    today_date = str(date.today())
        
    # create file that will append each record to
    concat_fasta_outfile = os.path.join(wd, 'aligned_sequences.alignment.fasta' )
    
    # first add the reference sequence
    ref_record = SeqIO.read(ref_path, 'fasta')
    with open(concat_fasta_outfile, 'w') as handle:
        SeqIO.write(ref_record, handle, 'fasta')
        
    # first append the non - modified records
    os.chdir(temp_alignment_dir)
    for file in glob.glob('*.fasta'):
        sample_name = file.split('.alignment.fasta')[0]
        if sample_name not in mod_seq_list:
            for record in AlignIO.read(file, 'fasta'):
                if record.id == sample_name:
                    with open(concat_fasta_outfile, 'a') as handle:
                        SeqIO.write(record, handle, 'fasta')
            
    # now append the records from the mod_seq_list
    if os.path.exists(mod_fasta_dir):
        os.chdir(mod_fasta_dir)
        for file in glob.glob('*.fasta'):
            record = SeqIO.read(file, 'fasta')
            with open(concat_fasta_outfile, 'a') as handle:
                SeqIO.write(record, handle, 'fasta')

    alignment = AlignIO.read(concat_fasta_outfile, 'fasta')
    print('******************')
    print('alignment has %d sequences' % len(alignment))
    print('sequence alignment length is %s' % alignment.get_alignment_length())
    print('*******************')



def main(wd, fasta_files_dir, ref_genome_path):

         
    ref = get_ref_seq_record_id(ref_genome_path= ref_genome_path)

    temp_dirs= create_temp_directories(wd)

    add_ref_genome(fasta_files_dir_path = fasta_files_dir, 
                                   wd = wd, 
                                   ref_genome_path = ref_genome_path, 
                                   fasta_temp_dir = temp_dirs['fasta_temp_dir'])  

    align_sequences(fasta_temp_dir = temp_dirs['fasta_temp_dir'], 
                    alignment_temp_dir = temp_dirs['alignment_temp_dir'], 
                    wd = wd) 

    insert = remove_insertions(alignment_temp_dir = temp_dirs['alignment_temp_dir'], 
                      fasta_temp_dir = temp_dirs['alignment_temp_dir'], 
                      ref_id = ref['ref_id'], 
                      ref_seq = ref['ref_seq'], 
                      ref_genome_path = ref_genome_path)

    record_deletions(alignment_temp_dir = temp_dirs['alignment_temp_dir'], 
                     ref_seq = ref['ref_seq'], 
                     insert_df= insert['insert_df'], 
                     mod_seq_list = insert['mod_seq_list'], 
                     wd = wd ) 

    mod_fasta_dir = create_dir_with_mod_fastas(alignment_temp_dir= temp_dirs['alignment_temp_dir'], 
                               mod_seq_list = insert['mod_seq_list'], 
                               wd = wd)
    
    create_concatenated_seq_records(wd= wd, 
                                    temp_alignment_dir =temp_dirs['alignment_temp_dir'], 
                                    mod_fasta_dir= mod_fasta_dir , 
                                    mod_seq_list = insert['mod_seq_list'],
                                   ref_path = ref_genome_path)

    delete_temp_directory(directory_name = temp_dirs['alignment_temp_dir'])
    delete_temp_directory(directory_name = temp_dirs['fasta_temp_dir'])
    delete_temp_directory(directory_name = mod_fasta_dir)
    os.chdir(wd)


    





