#!/usr/bin/env python
# coding: utf-8

# ### Parameters

# In[9]:


import os

REDUN_CUTOFF = 0.987
E_VALUE_THRESH = 0.04
BLASTN_DB = "16S_ribosomal_RNA"
#BLASTN_DB = "SILVA_138.1_SSU_tax_silva.fasta.blastdb"
DATA_DIR = "../data/result"
DATA_SUFFIX = ".seq"
OUT_DIR = os.path.join(DATA_DIR, 'output')
if not os.path.exists(OUT_DIR):
    os.mkdir(OUT_DIR)


# In[2]:


import Bio
from Bio import pairwise2
from Bio import SeqIO

def align_score_from_seqs(seq1, seq2):
    '''
    Global alignment of two sequences
    
    Return alignment score of aligned counts versus length of alignment
    '''
    alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)
    try:
        alignment = alignments[0]
    except IndexError:
        raise Exception('No alignment found between:\n{}\n{}'.format(seq1.seq, seq2.seq))
    return alignments[0].score / (alignments[0].end - alignments[0].start)

def align_score_fasta_files(file1, file2):
    '''
    Global alignment for two fasta files
    
    Return alignment score of aligned counts versus length of alignment
    '''
    try:
        seq1 = SeqIO.read(file1, "fasta")
    except IOError:
        raise Exception('No such file: {}'.format(file1))
    try:
        seq2 = SeqIO.read(file2, "fasta")
    except IOError:
        raise Exception('No such file: {}'.format(file2))
    return align_score_from_seqs(seq1, seq2)
    
def check_redundancy(scores, cutoff, seqnames):
    '''
    Given alignment scores to other sequences, return uniqueness of this input sequence:
    
    1: unique clone - score < cutoff to all other sequences
    2: similar clone - score >= cutoff to one or multiple other sequences
    3: repeat clone - score = 100% to one or multiple other sequences
    
    You can set the cutoff that determined to be similar or unique
    '''
    repeat_ids = []
    similar_ids = []
    unique_flag = True
    for i, x in enumerate(scores):
        if x == 1:
            repeat_ids.append(i)
            unique_flag = False
        elif x >= cutoff:
            similar_ids.append(i)
            unique_flag = False
    if unique_flag:
        return 'uniq'
    if repeat_ids:
        return 'repeat', [seqnames[i] for i in repeat_ids]
    return 'simi', [seqnames[i] for i in similar_ids]

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
def blastn_identity(file_query, db, evalue, file_out):
    '''
    Run BLAST and return the top-1 matched sequence with identity score
    '''
    blastn_cline = NcbiblastnCommandline(query=file_query, db=db, evalue=evalue, outfmt=5, out=file_out)
    print(blastn_cline)
    stdout, stderr = blastn_cline()
    
    blast_record = NCBIXML.read(open(file_out))
    alignment = blast_record.alignments[0]
    hsp = alignment.hsps[0]
    if hsp.expect < evalue:
        try:
            identity = hsp.identities / hsp.align_length
        except Exception as e:
            print('Error in calculating identity score for the top-1 alignment against DB: {}\n{}'.format(db, e))
            return
        return alignment.title, identity
    else:
        print('No identity score passed the E-VALUE threshold against DB: {}\n{}'.format(db))
        return
    


# ### Pairwise sequence similairty

# In[3]:


import time
import itertools
import numpy
import pandas

filenames = [x for x in os.listdir(DATA_DIR) if x.endswith(DATA_SUFFIX)]
files = [os.path.join(DATA_DIR, x) for x in filenames]
file_ids = range(len(files))
seqnames = [os.path.splitext(x)[0] for x in filenames]


# In[4]:


time0 = time.time()
align_score_df = pandas.DataFrame(data=0, columns=file_ids, index=file_ids, dtype=numpy.float64)
for i, j in itertools.combinations(file_ids, 2):
    align_score_df.loc[i,j] = align_score_df.loc[j,i] = align_score_fasta_files(files[i], files[j])
print('Elapsed time: ', time.time() - time0)


# ### Check redundancy: unique, similar or repeat sequence

# In[5]:


redun_series = align_score_df.apply(lambda x: check_redundancy(x, cutoff=REDUN_CUTOFF, seqnames=seqnames), axis=1)
redun_series.name = "Redundancy"


# ### BLAST for taxonomy

# In[10]:


tax_series = []
identity_series = []
for i in redun_series.index:
    tmp = blastn_identity(files[i], BLASTN_DB, E_VALUE_THRESH, os.path.join(OUT_DIR, os.path.split(files[i])[-1] + '.blastn.{}.xml').format(BLASTN_DB))
    if tmp is None:
        tax = identity = None
    else:
        tax, identity = tmp
    tax_series.append(tax)
    identity_series.append(identity)
tax_series = pandas.Series(data=tax_series, name='Taxonomy')
identity_series = pandas.Series(data=identity_series, name='Identity')


# ### Merge and output results

# In[11]:


seqname_series = pandas.Series(data=seqnames, name='Sequence')
result = pandas.concat([seqname_series, redun_series, tax_series, identity_series], axis=1)
result.to_excel(os.path.join(OUT_DIR, '16s_redun_tax_identity_{}.xlsx'.format(BLASTN_DB)))


# In[12]:


display(result)

