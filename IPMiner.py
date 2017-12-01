# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from sklearn import metrics
#from sklearn.metrics import confusion_matrix
import theano
#from pystruct.datasets import load_letters
#from pystruct.models import ChainCRF
#from pystruct.learners import OneSlackSSVM
from sklearn import svm, grid_search
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.cross_validation import train_test_split
from sklearn.calibration import CalibratedClassifierCV
from sklearn.cross_validation import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
import gzip
from random import randint
#import xgboost as xgb
import pandas as pd
import pdb
import os
import sys
import random
import argparse
from theano import tensor as T
# keras version is 0.1.2, please install this version of keras
#import keras
#keras_version =  keras.__version__
try:
    print 'please install Keras 0.1.2'
    sys.path.insert(0, '/usr/local/lib/python2.7/dist-packages/Keras-0.1.2-py2.7.egg')
except:
    print 'install keras 0.1.2'
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, AutoEncoder, Flatten, Merge
from keras.layers.normalization import BatchNormalization
from keras.layers.advanced_activations import PReLU
from keras.utils import np_utils, generic_utils
from keras.optimizers import SGD, RMSprop, Adadelta, Adagrad, Adam
from keras.layers import containers, normalization
#from keras.layers.convolutional import Convolution2D, MaxPooling2D
from keras.layers.recurrent import LSTM
from keras.layers.embeddings import Embedding
from keras import regularizers
from keras.constraints import maxnorm
from keras.optimizers import kl_divergence
#from sknn.mlp import Classifier, Layer



''' 
echo -e "3I5F\n2p4k\n2p4m" | while read I; do curl -s "http://www.rcsb.org/pdb/rest/customReport?pdbids=${I}&customReportColumns=structureId,chainId,entityId,sequence,db_id,
db_name&service=wsdisplay&format=csv"; done >result.csv

$ echo -e "3I5F\n2p4k\n2p4m" | while read I; do curl -s "http://www.rcsb.org/pdb/rest/customReport?pdbids=${I}
&customReportColumns=structureId,chainId,entityId,sequence,db_id,db_name&service=wsdisplay&format=text" 
| xsltproc stylesheet.xsl - ; done | fold -w 80   
def deep_learning_classifier(X_train, y_train):
    nn = Classifier(
    layers=[
        Layer("Rectifier", units=100),
        Layer("Linear")],
    learning_rate=0.02,
    n_iter=10)
    nn.fit(X_train, y_train)
    
    y_valid = nn.predict(X_valid)
    
    score = nn.score(X_test, y_test)
'''    
def get_uniq_pdb_protein_rna():
    protein_set = set()
    with open('ncRNA-protein/NegativePairs.csv', 'r') as fp:
        for line in fp:
            if 'Protein ID' in line:
                continue
            pro1, pro2 = line.rstrip().split('\t')
            protein_set.add(pro1.split('-')[0])
            protein_set.add(pro2.split('-')[0])
    
    with open('ncRNA-protein/PositivePairs.csv', 'r') as fp:
        for line in fp:
            if 'Protein ID' in line:
                continue
            pro1, pro2 = line.rstrip().split('\t')
            protein_set.add(pro1.split('-')[0])
            protein_set.add(pro2.split('-')[0])
    return protein_set  

def download_seq_from_PDB(protein_set, outfile_name):
    fw = open(outfile_name, 'w')
    for val in protein_set:
        cli_str = 'curl -s "http://www.rcsb.org/pdb/rest/customReport?pdbids='+ val +'&customReportColumns=structureId,chainId,sequence&service=wsdisplay&format=csv" >ncRNA-protein/tmpseq.csv'        
        cli_fp = os.popen(cli_str, 'r')
        cli_fp.close()
        #pdb.set_trace()
        f_in = open('ncRNA-protein/tmpseq.csv', 'r')
        for line in f_in:
            values = line.rstrip().split('<br />')
            for val in values:
                if 'structureId' in val:
                    continue
                if len(val) ==0:
                    continue
                pdbid, chainid, seq = val.split(',')
                fasa_name = pdbid[1:-1] + '-' + chainid[1:-1]
                fw.write('>' + fasa_name + '\n')
                fw.write(seq[1:-1] + '\n')
                #pdb.set_trace()
        f_in.close()
    fw.close()
    
def get_all_PDB_id():
    protein_set =  get_uniq_pdb_protein_rna()
    download_seq_from_PDB(protein_set, 'ncRNA-protein/all_seq.fa')

def get_protein_rna_id(inputfile):
    protein_set = set()
    with open(inputfile, 'r') as fp:
        for line in fp:
            if line[0] == '#':
                continue
            else:
                protein, rna = line.rstrip('\r\n').split()
                protein_set.add(protein.split('_')[0])
                protein_set.add(rna.split('_')[0])
    return protein_set

def judge_RNA_protein(seq):
    if 'U' in seq:
        return 'RNA'
    if all([c in 'AUGCTIN' for c in seq]):
        return 'RNA'
    else:
        return 'PROTEIN'

def generate_negative_samples_RPI2241_RPI369(seq_fasta, interaction_file, whole_file):
    seq_dict = read_fasta_file(seq_fasta)
    name_strand = {}
    type_dict = {}
    for key, tmpseq in seq_dict.iteritems():
        val, strand = key.split('-')
        seqtype = judge_RNA_protein(tmpseq)
        name_strand.setdefault(val, []).append(key)
        type_dict[key] = seqtype
        
    fw = open(whole_file, 'w')    
    existing_postive = set()
    with open(interaction_file, 'r') as fp:
        for line in fp:
            if line[0] == '#':
                continue
            pro1, pro2 = line.rstrip().split('\t')
            pro1 = pro1.replace('_', '-').upper()
            pro2 = pro2.replace('_', '-').upper()
            if type_dict[pro1] == 'RNA' and type_dict[pro2] == 'PROTEIN':
                existing_postive.add((pro2, pro1))
                fw.write(pro2 + '\t' + pro1 + '\t' + '1' + '\n')
            else:
                existing_postive.add((pro1, pro2))
                fw.write(pro1 + '\t' + pro2 + '\t' + '1' + '\n')
    
    #generate negative samples
    all_pairs = list(existing_postive)
    num_posi = len(all_pairs)
    nega_list  = []
    for val in all_pairs:
        pro1, pro2 = val
        for i in range(50):
            pro2_pare = pro2.split('-')[0]
            if name_strand.has_key(pro2_pare):
                for val in name_strand[pro2_pare]:
                    if type_dict[val] == 'PROTEIN' and val != pro1:
                        new_sele = (val, pro2)
                        if new_sele not in existing_postive:
                            #fw.write(val + '\t' + pro2 + '\t' + '0' + '\n')
                            nega_list.append(new_sele)
                            existing_postive.add(new_sele)
    
    random.shuffle(nega_list)
    for val in nega_list[:num_posi]:
        fw.write(val[0] + '\t' + val[1] + '\t' + '0' + '\n')
    fw.close()
        
    #protein_set.add(pro1.replace('_', '-').upper())
            #RNA_set.add(pro2.replace('_', '-').upper())


def get_RNA_protein_RPI2241_RPI369(seq_fasta, interaction_file, protein_file, rna_file):
    seq_dict = read_fasta_file(seq_fasta)
    RNA_set = set()
    protein_set = set()
    with open(interaction_file, 'r') as fp:
        for line in fp:
            if line[0] == '#':
                continue
            pro1, pro2, label = line.rstrip().split('\t')
            #protein_set.add(pro1.replace('_', '-').upper())
            #RNA_set.add(pro2.replace('_', '-').upper())
            protein_set.add(pro1)
            RNA_set.add(pro2)
            
    fw_pro = open(protein_file, 'w')
    for val in protein_set:
        if seq_dict.has_key(val):
            fw_pro.write('>' + val + '\n')
            fw_pro.write(seq_dict[val] + '\n')
        else:
            print val
    fw_pro.close()
    
    fw_pro = open(rna_file, 'w')
    for val in RNA_set:
        if seq_dict.has_key(val):
            fw_pro.write('>' + val + '\n')
            fw_pro.write(seq_dict[val].replace('N', '') + '\n')
        else:
            print val
    fw_pro.close()                 

def get_RPI2241_RPI369_seq():
    print 'downloading seqs'
    protein_set =  get_protein_rna_id('ncRNA-protein/RPI2241.txt')
    download_seq_from_PDB(protein_set, 'ncRNA-protein/RPI2241.fa')
    protein_set =  get_protein_rna_id('ncRNA-protein/RPI369.txt')
    download_seq_from_PDB(protein_set, 'ncRNA-protein/RPI369.fa')

def get_RPI2241_RPI369_ind_file():
    get_RNA_protein_RPI2241_RPI369('ncRNA-protein/RPI2241.fa', 'ncRNA-protein/RPI2241_all.txt', 'ncRNA-protein/RPI2241_protein.fa', 'ncRNA-protein/RPI2241_rna.fa')
    get_RNA_protein_RPI2241_RPI369('ncRNA-protein/RPI369.fa', 'ncRNA-protein/RPI369_all.txt', 'ncRNA-protein/RPI369_protein.fa', 'ncRNA-protein/RPI369_rna.fa')
    
def read_fasta_file(fasta_file):
    seq_dict = {}    
    fp = open(fasta_file, 'r')
    name = ''
    #pdb.set_trace()
    for line in fp:
        #let's discard the newline at the end (if any)
        line = line.rstrip()
        #distinguish header from sequence
        if line[0]=='>': #or line.startswith('>')
            #it is the header
            name = line[1:].upper() #discarding the initial >
            seq_dict[name] = ''
        else:
            #it is sequence
            seq_dict[name] = seq_dict[name] + line
    fp.close()
    
    return seq_dict
    
def get_RNA_protein():
    seq_dict = read_fasta_file('ncRNA-protein/all_seq.fa')
    RNA_set = set()
    protein_set = set()
    with open('ncRNA-protein/NegativePairs.csv', 'r') as fp:
        for line in fp:
            if 'Protein ID' in line:
                continue
            pro1, pro2 = line.rstrip().split('\t')
            protein_set.add(pro1)
            RNA_set.add(pro2)
    
    with open('ncRNA-protein/PositivePairs.csv', 'r') as fp:
        for line in fp:
            if 'Protein ID' in line:
                continue
            pro1, pro2 = line.rstrip().split('\t')
            protein_set.add(pro1)
            RNA_set.add(pro2) 
    
    fw_pro = open('ncRNA-protein/protein_seq.fa', 'w')
    for val in protein_set:
        if seq_dict.has_key(val):
            fw_pro.write('>' + val + '\n')
            fw_pro.write(seq_dict[val] + '\n')
        else:
            print val
    fw_pro.close()
    
    fw_pro = open('ncRNA-protein/RNA_seq.fa', 'w')
    for val in RNA_set:
        if seq_dict.has_key(val):
            fw_pro.write('>' + val + '\n')
            fw_pro.write(seq_dict[val] + '\n')
        else:
            print val
    fw_pro.close() 

def read_name_from_fasta(fasta_file):
    name_list = []
    fp = open(fasta_file, 'r')
    for line in fp:
        if line[0] == '>':
            name = line.rstrip('\r\n')[1:]
            name_list.append(name.upper())
    fp.close()
    return name_list

def get_noncode_seq():
    ncRNA_seq_dict = {}
    head  = True
    name = ''
    #pdb.set_trace()
    with open('ncRNA-protein/ncrna_NONCODE[v3.0].fasta', 'r') as fp:
        for line in fp:
            if head:
                head =False
                continue
            line = line.rstrip()
            if line == 'sequence':
                continue
            if line[0] == '>':
                name1 = line.split('|')
                name = name1[0][1:].strip()
                ncRNA_seq_dict[name] = ''
            else:
                #it is sequence
                ncRNA_seq_dict[name] = ncRNA_seq_dict[name] + line
            
    return ncRNA_seq_dict

def get_npinter_protein_seq():
    pro_dict = {}
    target_dir = 'ncRNA-protein/uniprot_seq/'
    files = os.listdir(target_dir)
    for file_name in files:
        protein_name = file_name.split('.')[0]
        with open(target_dir + file_name, 'r') as fp:
            for line in fp:
                line = line.rstrip()
                if line[0] == '>':
                    pro_dict[protein_name] = ''
                else:
                    pro_dict[protein_name] = pro_dict[protein_name] + line
                    
    return pro_dict   

def read_RNA_pseaac_fea(name_list, pseaac_file='ncRNA-protein/RNA_pse.csv'):
    print pseaac_file
    data = {}
    fp = open(pseaac_file, 'r')
    index = 0
    for line in fp:
        values = line.rstrip('\r\n').split(',')
        data[name_list[index]] = [float(val) for val in values]
        index = index + 1
    fp.close()
    return data      
    
def read_RNA_graph_feature(name_list, graph_file='ncRNA-protein/RNA_seq.gz.feature', fea_imp = None):
    print graph_file
    data = {}
    fea_len = 32768
    fp = open(graph_file, 'r')
    index = 0
    for line in fp:
        tmp_data = [0] * fea_len
        values = line.split()
        for value in values:
            val = value.split(':')
            tmp_data[int(val[0])] = float(val[1])
        if fea_imp is None:
            data[name_list[index]] = tmp_data
        else:
            data[name_list[index]] = [tmp_data[val] for val in fea_imp]
        index = index + 1
    fp.close()
    return data   

def read_protein_feature(protein_fea_file = 'ncRNA-protein/trainingSetFeatures.csv'):
    print protein_fea_file
    feature_dict = {}
    df = pd.read_csv(protein_fea_file)
    X = df.values.copy()
    for val in X:
        feature_dict[val[0].upper()] = val[2:].tolist()
    
    #pdb.set_trace()
    return feature_dict

def read_lncRNA_protein_feature(protein_fea_file = 'ncRNA-protein/trainingSetFeatures.csv'):
    print protein_fea_file
    feature_dict = {}
    df = pd.read_csv(protein_fea_file)
    X = df.values.copy()
    for val in X:
        feature_dict[val[0]] = val[2:].tolist()
    
    #pdb.set_trace()
    return feature_dict

def get_4_trids():
    nucle_com = []
    chars = ['A', 'C', 'G', 'U']
    base=len(chars)
    end=len(chars)**4
    for i in range(0,end):
        n=i
        ch0=chars[n%base]
        n=n/base
        ch1=chars[n%base]
        n=n/base
        ch2=chars[n%base]
        n=n/base
        ch3=chars[n%base]
        nucle_com.append(ch0 + ch1 + ch2 + ch3)
    return  nucle_com   

def get_3_trids():
    nucle_com = []
    chars = ['A', 'C', 'G', 'U']
    base=len(chars)
    end=len(chars)**3
    for i in range(0,end):
        n=i
        ch0=chars[n%base]
        n=n/base
        ch1=chars[n%base]
        n=n/base
        ch2=chars[n%base]
        nucle_com.append(ch0 + ch1 + ch2)
    return  nucle_com  
           
def get_4_nucleotide_composition(tris, seq, pythoncount = True):
    seq_len = len(seq)
    tri_feature = []
    
    if pythoncount:
        for val in tris:
            num = seq.count(val)
            tri_feature.append(float(num)/seq_len)
    else:
        k = len(tris[0])
        tmp_fea = [0] * len(tris)
        for x in range(len(seq) + 1- k):
            kmer = seq[x:x+k]
            if kmer in tris:
                ind = tris.index(kmer)
                tmp_fea[ind] = tmp_fea[ind] + 1
        tri_feature = [float(val)/seq_len for val in tmp_fea]
        #pdb.set_trace()        
    return tri_feature

def TransDict_from_list(groups):
    transDict = dict()
    tar_list = ['0', '1', '2', '3', '4', '5', '6']
    result = {}
    index = 0
    for group in groups:
        g_members = sorted(group) #Alphabetically sorted list
        for c in g_members:
            # print('c' + str(c))
            # print('g_members[0]' + str(g_members[0]))
            result[c] = str(tar_list[index]) #K:V map, use group's first letter as represent.
        index = index + 1
    return result

def translate_sequence (seq, TranslationDict):
    '''
    Given (seq) - a string/sequence to translate,
    Translates into a reduced alphabet, using a translation dict provided
    by the TransDict_from_list() method.
    Returns the string/sequence in the new, reduced alphabet.
    Remember - in Python string are immutable..

    '''
    import string
    from_list = []
    to_list = []
    for k,v in TranslationDict.items():
        from_list.append(k)
        to_list.append(v)
    # TRANS_seq = seq.translate(str.maketrans(zip(from_list,to_list)))
    TRANS_seq = seq.translate(string.maketrans(str(from_list), str(to_list)))
    #TRANS_seq = maketrans( TranslationDict, seq)
    return TRANS_seq

def get_protein_trids(seq, group_dict):

    #protein='MQNEEDACLEAGYCLGTTLSSWRLHFMEEQSQSTMLMGIGIGALLTLAFVGIFFFVYRR'
    tran_seq = translate_sequence (seq, group_dict)
    #pdb.set_trace()
    return tran_seq

def get_3_protein_trids():
    nucle_com = []
    chars = ['0', '1', '2', '3', '4', '5', '6']
    base=len(chars)
    end=len(chars)**3
    for i in range(0,end):
        n=i
        ch0=chars[n%base]
        n=n/base
        ch1=chars[n%base]
        n=n/base
        ch2=chars[n%base]
        nucle_com.append(ch0 + ch1 + ch2)
    return  nucle_com

def get_4_protein_trids():
    nucle_com = []
    chars = ['0', '1', '2', '3', '4', '5', '6']
    base=len(chars)
    end=len(chars)**4
    for i in range(0,end):
        n=i
        ch0=chars[n%base]
        n=n/base
        ch1=chars[n%base]
        n=n/base
        ch2=chars[n%base]
        n=n/base
        ch3=chars[n%base]
        nucle_com.append(ch0 + ch1 + ch2 + ch3)
    return  nucle_com

def get_NPinter_interaction():
    RNA_set = set()
    protein_set = set()
    with open('ncRNA-protein/NPInter10412_dataset.txt', 'r') as fp:
        head  = True
        for line in fp:
            if head:
                head = False
                continue
            pro1, pro1_len, pro2, pro2_len, org = line.rstrip().split('\t')
            protein_set.add(pro2)
            RNA_set.add(pro1)
    pro_dict = get_npinter_protein_seq()        
    fw_pro = open('ncRNA-protein/NPinter_protein_seq.fa', 'w')
    for val in protein_set:
        if pro_dict.has_key(val):
            fw_pro.write('>' + val + '\n')
            fw_pro.write(pro_dict[val] + '\n')
        else:
            print val
    fw_pro.close()
    ncRNA_dict = get_noncode_seq()
    fw_pro = open('ncRNA-protein/NPinter_RNA_seq.fa', 'w')
    for val in RNA_set:
        if ncRNA_dict.has_key(val):
            fw_pro.write('>' + val + '\n')
            seq = ncRNA_dict[val].replace('T', 'U')
            #seq = seq.replace('N', '')get_RPI2241_RPI369_seq()
            fw_pro.write( seq + '\n')
        else:
            print val
    fw_pro.close()     

def get_own_lncRNA_protein(datafile = 'ncRNA-protein/lncRNA-protein-488.txt'):
    protein_seq = {}
    RNA_seq = {}
    interaction_pair = {}
    with open(datafile, 'r') as fp:
        for line in fp:
            if line[0] == '>':
                values = line[1:].strip().split('|')
                label = values[1]
                name = values[0].split('_')
                protein = name[0] + '-' + name[1]
                RNA = name[0] + '-' + name[2]
                if label == 'interactive':
                    interaction_pair[(protein, RNA)] = 1
                else:
                    interaction_pair[(protein, RNA)] = 0
                index  = 0
            else:
                seq = line[:-1]
                if index == 0:
                    protein_seq[protein] = seq
                else:
                    RNA_seq[RNA] = seq
                index = index + 1
    #pdb.set_trace()
    fw = open('ncRNA-protein/lncRNA_protein.fa', 'w')
    for key, val in protein_seq.iteritems():
        fw.write('>' + key + '\n')
        fw.write(val + '\n')
        
    fw.close()
    
    fw = open('ncRNA-protein/lncRNA_RNA.fa', 'w')
    for key, val in RNA_seq.iteritems():
        fw.write('>' + key + '\n')
        fw.write(val.replace('N', '') + '\n')
        
    fw.close()
    '''
    cli_str = "python ProFET/ProFET/feat_extract/pipeline.py --trainingSetDir 'ncRNA-protein/lncRNA-protein/' \
    --trainFeatures True --resultsDir 'ncRNA-protein/lncRNA-protein/' --classType file"
    fcli = os.popen(cli_str, 'r')
    fcli.close()
    '''
def read_name_from_lncRNA_fasta(fasta_file):
    name_list = []
    fp = open(fasta_file, 'r')
    for line in fp:
        if line[0] == '>':
            name = line.rstrip('\r\n')[1:]
            name_list.append(name)
    fp.close()
    return name_list

def prepare_complex_feature(rna_file, protein_file, seperate = False):
    RNA_seq_dict = read_fasta_file(rna_file)
    protein_seq_dict = read_fasta_file(protein_file)
    groups = ['AGV', 'ILFP', 'YMTS', 'HNQW', 'RK', 'DE', 'C']
    group_dict = TransDict_from_list(groups)
    protein_tris = get_3_protein_trids()
    tris = get_4_trids()
    pairs = []
    train = []
    label = []
    for RNA, RNA_seq in RNA_seq_dict.iteritems():
        RNA_seq = RNA_seq.replace('T', 'U')
        for protein, protein_seq in protein_seq_dict.iteritems():
            pairs.append((RNA, protein))
            protein_seq1 = translate_sequence (protein_seq, group_dict)
            RNA_tri_fea = get_4_nucleotide_composition(tris, RNA_seq, pythoncount =False)
            protein_tri_fea = get_4_nucleotide_composition(protein_tris, protein_seq1, pythoncount =False)
            if seperate:
                tmp_fea = (protein_tri_fea, RNA_tri_fea)
            else:
                tmp_fea = protein_tri_fea + RNA_tri_fea
            train.append(tmp_fea)            
    
    return np.array(train), pairs

def prepare_RPI488_feature(extract_only_posi = False, 
                                   pseaac_file = None, deepmind = False, seperate = False, chem_fea = True):
    print 'RPI488 dataset'
    interaction_pair = {}
    RNA_seq_dict = {}
    protein_seq_dict = {}
    with open('ncRNA-protein/lncRNA-protein-488.txt', 'r') as fp:
        for line in fp:
            if line[0] == '>':
                values = line[1:].strip().split('|')
                label = values[1]
                name = values[0].split('_')
                protein = name[0] + '-' + name[1]
                RNA = name[0] + '-' + name[2]
                if label == 'interactive':
                    interaction_pair[(protein, RNA)] = 1
                else:
                    interaction_pair[(protein, RNA)] = 0
                index  = 0
            else:
                seq = line[:-1]
                if index == 0:
                    protein_seq_dict[protein] = seq
                else:
                    RNA_seq_dict[RNA] = seq
                index = index + 1
    #name_list = read_name_from_lncRNA_fasta('ncRNA-protein/lncRNA_RNA.fa')           
    groups = ['AGV', 'ILFP', 'YMTS', 'HNQW', 'RK', 'DE', 'C']
    group_dict = TransDict_from_list(groups)
    protein_tris = get_3_protein_trids()
    tris = get_4_trids()
    #tris3 = get_3_trids()
    train = []
    label = []
    chem_fea = []
    for key, val in interaction_pair.iteritems():
        protein, RNA = key[0], key[1]
        #pdb.set_trace()
        if RNA_seq_dict.has_key(RNA) and protein_seq_dict.has_key(protein): #and protein_fea_dict.has_key(protein) and RNA_fea_dict.has_key(RNA):
            label.append(val)
            RNA_seq = RNA_seq_dict[RNA]
            protein_seq = translate_sequence (protein_seq_dict[protein], group_dict)
            if deepmind:
                RNA_tri_fea = get_RNA_seq_concolutional_array(RNA_seq)
                protein_tri_fea = get_RNA_seq_concolutional_array(protein_seq) 
                train.append((RNA_tri_fea, protein_tri_fea))
            else:
                #pdb.set_trace()
                RNA_tri_fea = get_4_nucleotide_composition(tris, RNA_seq, pythoncount =False)
                protein_tri_fea = get_4_nucleotide_composition(protein_tris, protein_seq, pythoncount =False)
                #RNA_tri3_fea = get_4_nucleotide_composition(tris3, RNA_seq, pythoncount =False)
                #RNA_fea = [RNA_fea_dict[RNA][ind] for ind in fea_imp]
                #tmp_fea = protein_fea_dict[protein] + tri_fea #+ RNA_fea_dict[RNA]
                if seperate:
                    tmp_fea = (protein_tri_fea, RNA_tri_fea)
                    #chem_tmp_fea = (protein_fea_dict[protein], RNA_fea_dict[RNA])
                else:
                    tmp_fea = protein_tri_fea + RNA_tri_fea
                    #chem_tmp_fea = protein_fea_dict[protein] + RNA_fea_dict[RNA] 
                train.append(tmp_fea)
                #chem_fea.append(chem_tmp_fea)
        else:
            print RNA, protein   
    
    return np.array(train), label             

def prepare_RPIntDB_feature(extract_only_posi = False, 
                                   pseaac_file = None, deepmind = False, seperate = False, chem_fea = True):
    groups = ['AGV', 'ILFP', 'YMTS', 'HNQW', 'RK', 'DE', 'C']
    group_dict = TransDict_from_list(groups)
    protein_tris = get_3_protein_trids()
    tris = get_4_trids()
    #tris3 = get_3_trids()
    train = []
    label = []
    chem_fea = []
    head = True
    with open('ncRNA-protein/RPIntDB_interactions_new.txt', 'r') as fp:
        for line in fp:
            if head:
                head = False
                continue
            values = line.rstrip('\r\n').split('\t')
            protein = values[0]
            protein_seq = values[1]
            RNA = values[2]
            RNA_seq = values[3] 
            label.append(1)
            RNA_tri_fea = get_4_nucleotide_composition(tris, RNA_seq, pythoncount =False)
            protein_tri_fea = get_4_nucleotide_composition(protein_tris, protein_seq, pythoncount =False)
    
            if seperate:
                tmp_fea = (protein_tri_fea, RNA_tri_fea)
                #chem_tmp_fea = (protein_fea_dict[protein], RNA_fea_dict[RNA])
            else:
                tmp_fea = protein_tri_fea + RNA_tri_fea
                #chem_tmp_fea = protein_fea_dict[protein] + RNA_fea_dict[RNA] 
            train.append(tmp_fea) 
    
    return np.array(train), label
        
def prepare_RPI2241_369_feature(rna_fasta_file, data_file, protein_fasta_file, extract_only_posi = False, 
                                graph = False, deepmind = False, seperate = False, chem_fea = True):
    seq_dict = read_fasta_file(rna_fasta_file)
    protein_seq_dict = read_fasta_file(protein_fasta_file)
    groups = ['AGV', 'ILFP', 'YMTS', 'HNQW', 'RK', 'DE', 'C']
    group_dict = TransDict_from_list(groups)
    protein_tris = get_3_protein_trids()
    tris = get_4_trids()
    train = []
    label = []
    chem_fea = []
    #posi_set = set()
    #pro_set = set()
    with open(data_file, 'r') as fp:
        for line in fp:
            if line[0] == '#':
                continue
            protein, RNA, tmplabel = line.rstrip('\r\n').split('\t')
            if seq_dict.has_key(RNA) and protein_seq_dict.has_key(protein): 
                label.append(int(tmplabel))
                RNA_seq = seq_dict[RNA]
                protein_seq = translate_sequence (protein_seq_dict[protein], group_dict)
                if deepmind:
                    RNA_tri_fea = get_RNA_seq_concolutional_array(RNA_seq)
                    protein_tri_fea = get_RNA_seq_concolutional_array(protein_seq) 
                    train.append((RNA_tri_fea, protein_tri_fea))
                else:
                    RNA_tri_fea = get_4_nucleotide_composition(tris, RNA_seq, pythoncount =False)
                    protein_tri_fea = get_4_nucleotide_composition(protein_tris, protein_seq, pythoncount =False)
                    if seperate:
                        tmp_fea = (protein_tri_fea, RNA_tri_fea)
                    else:
                        tmp_fea = protein_tri_fea + RNA_tri_fea
                    train.append(tmp_fea)
            else:
                print RNA, protein  
            
    return np.array(train), label

def get_npinter_interaction(uniq_nid_dict = {}):
    pair_interaction = set()
    org_list  = []
    pair_list = []
    with open('ncRNA-protein/NPInter10412_dataset.txt', 'r') as fp:
        head  = True
        for line in fp:
            if head:
                head = False
                continue
            RNA, RNA_len, protein, protein_len, org = line.rstrip().split('\t')
            org_list.append(org)
            pair_list.append((RNA, protein))
            if uniq_nid_dict.has_key(RNA):
                uniq_name = uniq_nid_dict[RNA]
                pair_interaction.add((uniq_name, protein))
    return pair_list, org_list

def get_RPI367_interaction():
    org_list  = []
    with open('ncRNA-protein/RPI367.txt', 'r') as fp:
        head  = True
        for line in fp:
            if head:
                head = False
                continue
            values = line.rstrip().split('\t')
            org_list.append(values[2])
    return org_list
    
def prepare_RPI367_feature(extract_only_posi = True, graph = False, deepmind = False, seperate = False, chem_fea = False):
    print 'RPI367 data'
    seq_dict = {}
    uniq_nid_dict  = {}
    fzip = open('ncRNA-protein/ncrna_NONCODE[v1.0].fasta', 'r')
    for line in fzip:
        if line[0] == '#':
            continue
        if line[0] == '>':
            values = line.rstrip('\r\n').split(',')
            name = values[1]
            uniq_nid_dict[name] = values[0][1:]
        else:
            seq_dict[name] = line[:-1].replace('T', 'U')
    fzip.close()
    protein_seq_dict = get_npinter_protein_seq() 
    groups = ['AGV', 'ILFP', 'YMTS', 'HNQW', 'RK', 'DE', 'C']
    group_dict = TransDict_from_list(groups)
    protein_tris = get_3_protein_trids()
    #pdb.set_trace()
    train = []
    label = []
    posi_set = set()
    pro_set = set()
    tris = get_4_trids()
    with open('ncRNA-protein/RPI367.txt', 'r') as fp:
        head  = True
        for line in fp:
            if head:
                head = False
                continue
            #RNA, RNA_len, protein, protein_len, org = line.rstrip().split('\t')
            #RNA = RNA.upper()
            values = line.rstrip('\r\n').split('\t')
            protein = values[-1]
            RNA = values[-2]
            protein = protein.upper()
            posi_set.add((RNA, protein))
            pro_set.add(protein)
            if seq_dict.has_key(RNA) and protein_seq_dict.has_key(protein):
                label.append(1)
                RNA_seq = seq_dict[RNA]
                protein_seq = translate_sequence (protein_seq_dict[protein], group_dict)
                if deepmind:
                    RNA_tri_fea = get_RNA_seq_concolutional_array(RNA_seq)
                    protein_tri_fea = get_RNA_seq_concolutional_array(protein_seq) 
                    train.append((RNA_tri_fea, protein_tri_fea))
                else:
                    RNA_tri_fea = get_4_nucleotide_composition(tris, RNA_seq, pythoncount =False)
                    protein_tri_fea = get_4_nucleotide_composition(protein_tris, protein_seq, pythoncount =False)
                    if seperate:
                        tmp_fea = (protein_tri_fea, RNA_tri_fea)

                    else:
                        tmp_fea = protein_tri_fea + RNA_tri_fea

                    train.append(tmp_fea)

            else:
                print RNA, protein
        
    return np.array(train), label

def read_orf_seq(fasta_file, RNA = False):
    protein_seq_dict = {} 
    with open(fasta_file, 'r') as fp:
        for line in fp:
            line = line.rstrip()
            if line[0] == '>':
                name1 = line.split()
                name = name1[0][1:].strip()
                protein_seq_dict[name] = ''
            else:
                if RNA:
                    line = line.replace('T', 'U')
                protein_seq_dict[name] = protein_seq_dict[name] + line
    
    return protein_seq_dict

def read_orf_interaction(interaction_file):
    interacton_pair = []
    with open(interaction_file, 'r') as fp:
        head  = True
        for line in fp:
            if head:
                head = False
                continue
            values = line.rstrip().split()
            protein, RNA = values[0].split('_') 
            interacton_pair.append((protein, RNA))
    return interacton_pair

def prepare_RPI13254_feature(deepmind = False, seperate = False, extract_only_posi = False, indep_test = False):
    protein_seq_dict = read_orf_seq('ncRNA-protein/RPI13254_RNA_seq.fa')
    RNA_seq_dict = read_orf_seq('ncRNA-protein/RPI13254_protein_seq.fa', RNA = True)
    positive_pairs = read_orf_interaction('ncRNA-protein/RPI13254_positive.txt')
    negative_pairs = read_orf_interaction('ncRNA-protein/RPI13254_negative.txt')
    groups = ['AGV', 'ILFP', 'YMTS', 'HNQW', 'RK', 'DE', 'C']
    group_dict = TransDict_from_list(groups)
    protein_tris = get_3_protein_trids()
    tris = get_4_trids()
    train = []
    label = []
    if not indep_test:
        random.shuffle(positive_pairs)
        nega_num =len(negative_pairs)
        positive_pairs = positive_pairs[:nega_num]
    #pdb.set_trace()
    for val in positive_pairs:
        protein, RNA = val
        if RNA_seq_dict.has_key(RNA) and protein_seq_dict.has_key(protein):
            label.append(1)
            #RNA_fea = [RNA_fea_dict[RNA][ind] for ind in fea_imp]
            RNA_seq = RNA_seq_dict[RNA]
            protein_seq = translate_sequence (protein_seq_dict[protein], group_dict)
            if deepmind:
                RNA_tri_fea = get_RNA_seq_concolutional_array(RNA_seq)
                protein_tri_fea = get_RNA_seq_concolutional_array(protein_seq) 
                train.append((RNA_tri_fea, protein_tri_fea))
            else:
                RNA_tri_fea = get_4_nucleotide_composition(tris, RNA_seq, pythoncount =False)
                protein_tri_fea = get_4_nucleotide_composition(protein_tris, protein_seq, pythoncount =False)
                if seperate:
                    tmp_fea = (protein_tri_fea, RNA_tri_fea)

                else:
                    tmp_fea = protein_tri_fea + RNA_tri_fea

                train.append(tmp_fea)

        else:
            print RNA, protein
    if not extract_only_posi:        
        for val in negative_pairs:
            protein, RNA = val
            if RNA_seq_dict.has_key(RNA) and protein_seq_dict.has_key(protein):
                label.append(0)
                #RNA_fea = [RNA_fea_dict[RNA][ind] for ind in fea_imp]
                RNA_seq = RNA_seq_dict[RNA]
                protein_seq = translate_sequence (protein_seq_dict[protein], group_dict)
                if deepmind:
                    RNA_tri_fea = get_RNA_seq_concolutional_array(RNA_seq)
                    protein_tri_fea = get_RNA_seq_concolutional_array(protein_seq) 
                    train.append((RNA_tri_fea, protein_tri_fea))
                else:
                    RNA_tri_fea = get_4_nucleotide_composition(tris, RNA_seq, pythoncount =False)
                    protein_tri_fea = get_4_nucleotide_composition(protein_tris, protein_seq, pythoncount =False)
                    if seperate:
                        tmp_fea = (protein_tri_fea, RNA_tri_fea)
    
                    else:
                        tmp_fea = protein_tri_fea + RNA_tri_fea
    
                    train.append(tmp_fea)
    
            else:
                print RNA, protein
            
    return np.array(train), label    

def prepare_NPinter_feature(extract_only_posi = False, graph = False, deepmind = False, seperate = False, chem_fea = True):
    print 'NPinter data'
    name_list = read_name_from_fasta('ncRNA-protein/NPinter_RNA_seq.fa')
    seq_dict = read_fasta_file('ncRNA-protein/NPinter_RNA_seq.fa')
    protein_seq_dict = read_fasta_file('ncRNA-protein/NPinter_protein_seq.fa')
    groups = ['AGV', 'ILFP', 'YMTS', 'HNQW', 'RK', 'DE', 'C']
    group_dict = TransDict_from_list(groups)
    protein_tris = get_3_protein_trids()
    #pdb.set_trace()
    train = []
    label = []
    chem_fea = []
    posi_set = set()
    pro_set = set()
    tris = get_4_trids()
    with open('ncRNA-protein/NPInter10412_dataset.txt', 'r') as fp:
        head  = True
        for line in fp:
            if head:
                head = False
                continue
            RNA, RNA_len, protein, protein_len, org = line.rstrip().split('\t')
            RNA = RNA.upper()
            protein = protein.upper()
            posi_set.add((RNA, protein))
            pro_set.add(protein)
            if seq_dict.has_key(RNA) and protein_seq_dict.has_key(protein):
                label.append(1)
                #RNA_fea = [RNA_fea_dict[RNA][ind] for ind in fea_imp]
                RNA_seq = seq_dict[RNA]
                protein_seq = translate_sequence (protein_seq_dict[protein], group_dict)
                if deepmind:
                    RNA_tri_fea = get_RNA_seq_concolutional_array(RNA_seq)
                    protein_tri_fea = get_RNA_seq_concolutional_array(protein_seq) 
                    train.append((RNA_tri_fea, protein_tri_fea))
                else:
                    RNA_tri_fea = get_4_nucleotide_composition(tris, RNA_seq, pythoncount =False)
                    protein_tri_fea = get_4_nucleotide_composition(protein_tris, protein_seq, pythoncount =False)
                    if seperate:
                        tmp_fea = (protein_tri_fea, RNA_tri_fea)

                    else:
                        tmp_fea = protein_tri_fea + RNA_tri_fea

                    train.append(tmp_fea)

            else:
                print RNA, protein
    
    if not extract_only_posi:
        pro_list = list(pro_set)   
        total_pro_len = len(pro_list)       
        # get negative data
        with open('ncRNA-protein/NPInter10412_dataset.txt', 'r') as fp:
            head  = True
            for line in fp:
                if head:
                    head = False
                    continue
                RNA, RNA_len, protein, protein_len, org = line.rstrip().split('\t')
                RNA = RNA.upper()
                protein = protein.upper()
                for val in range(50):
                    random_choice = randint(0,total_pro_len-1)
                    select_pro = pro_list[random_choice]
                    selec_nega= (RNA, select_pro)
                    if selec_nega not in posi_set:
                        posi_set.add(selec_nega)
                        #print selec_nega
                        break
                        
                if seq_dict.has_key(RNA) and protein_seq_dict.has_key(select_pro): #and RNA_fea_dict.has_key(RNA) and protein_fea_dict.has_key(select_pro) :
                    label.append(0)
                    #RNA_fea = [RNA_fea_dict[RNA][ind] for ind in fea_imp]
                    RNA_seq = seq_dict[RNA]
                    protein_seq = translate_sequence (protein_seq_dict[select_pro], group_dict)
                    if deepmind:
                        RNA_tri_fea = get_RNA_seq_concolutional_array(RNA_seq)
                        protein_tri_fea = get_RNA_seq_concolutional_array(protein_seq) 
                        train.append((RNA_tri_fea, protein_tri_fea))
                    else:
                        RNA_tri_fea = get_4_nucleotide_composition(tris, RNA_seq, pythoncount =False)
                        protein_tri_fea = get_4_nucleotide_composition(protein_tris, protein_seq, pythoncount =False)
                    if seperate:
                        tmp_fea = (protein_tri_fea, RNA_tri_fea)
                    else:
                        tmp_fea = protein_tri_fea + RNA_tri_fea
                    train.append(tmp_fea)
                    #chem_fea.append(chem_tmp_fea)
                else:
                    print RNA, protein    
        #for key, val in RNA_fea_dict.iteritems():
            
            
    return np.array(train), label

def prepare_RPI1807_feature(graph = False, deepmind = False, seperate = False, chem_fea = True):
    print 'RPI-Pred data'
    #name_list = read_name_from_fasta('ncRNA-protein/RNA_seq.fa')
    seq_dict = read_fasta_file('ncRNA-protein/RPI1807_RNA_seq.fa')
    protein_seq_dict = read_fasta_file('ncRNA-protein/RPI1807_protein_seq.fa')
    groups = ['AGV', 'ILFP', 'YMTS', 'HNQW', 'RK', 'DE', 'C']
    group_dict = TransDict_from_list(groups)
    protein_tris = get_3_protein_trids()
    tris = get_4_trids()
    #pdb.set_trace()
    train = []
    label = []
    chem_fea = []
    #pdb.set_trace()
    with open('ncRNA-protein/RPI1807_PositivePairs.csv', 'r') as fp:
        for line in fp:
            if 'Protein ID' in line:
                continue
            pro1, pro2 = line.rstrip().split('\t')
            pro1 = pro1.upper()
            pro2 = pro2.upper()
            if seq_dict.has_key(pro2) and protein_seq_dict.has_key(pro1):#and protein_fea_dict.has_key(pro1) and RNA_fea_dict.has_key(pro2):
                label.append(1)
                RNA_seq = seq_dict[pro2]
                protein_seq = translate_sequence (protein_seq_dict[pro1], group_dict)
                if deepmind:
                    RNA_tri_fea = get_RNA_seq_concolutional_array(RNA_seq)
                    protein_tri_fea = get_RNA_seq_concolutional_array(protein_seq) 
                    train.append((RNA_tri_fea, protein_tri_fea))
                else:
                    RNA_tri_fea = get_4_nucleotide_composition(tris, RNA_seq, pythoncount =False)
                    protein_tri_fea = get_4_nucleotide_composition(protein_tris, protein_seq, pythoncount =False)
                    #RNA_fea = [RNA_fea_dict[RNA][ind] for ind in fea_imp]
                    #tmp_fea = protein_fea_dict[protein] + tri_fea #+ RNA_fea_dict[RNA]
                    if seperate:
                        tmp_fea = (protein_tri_fea, RNA_tri_fea)
                        #chem_tmp_fea = (protein_fea_dict[pro1], RNA_fea_dict[pro2])
                    else:
                        tmp_fea = protein_tri_fea + RNA_tri_fea
                        #chem_tmp_fea = protein_fea_dict[pro1] + RNA_fea_dict[pro2] 
                    train.append(tmp_fea)
                    #chem_fea.append(chem_tmp_fea)
            else:
                print pro1, pro2
    with open('ncRNA-protein/RPI1807_NegativePairs.csv', 'r') as fp:
        for line in fp:
            if 'Protein ID' in line:
                continue
            pro1, pro2 = line.rstrip().split('\t')
            pro1 = pro1.upper()
            pro2 = pro2.upper()            
            if seq_dict.has_key(pro2) and protein_seq_dict.has_key(pro1): #and protein_fea_dict.has_key(pro1) and RNA_fea_dict.has_key(pro2):
                label.append(0)
                RNA_seq = seq_dict[pro2]
                protein_seq = translate_sequence (protein_seq_dict[pro1], group_dict)
                if deepmind:
                    RNA_tri_fea = get_RNA_seq_concolutional_array(RNA_seq)
                    protein_tri_fea = get_RNA_seq_concolutional_array(protein_seq) 
                    train.append((RNA_tri_fea, protein_tri_fea))
                else:
                    RNA_tri_fea = get_4_nucleotide_composition(tris, RNA_seq, pythoncount =False)
                    protein_tri_fea = get_4_nucleotide_composition(protein_tris, protein_seq, pythoncount =False)
                    #RNA_fea = [RNA_fea_dict[RNA][ind] for ind in fea_imp]
                    #tmp_fea = protein_fea_dict[protein] + tri_fea #+ RNA_fea_dict[RNA]
                    if seperate:
                        tmp_fea = (protein_tri_fea, RNA_tri_fea)
                        #chem_tmp_fea = (protein_fea_dict[pro1], RNA_fea_dict[pro2])
                    else:
                        tmp_fea = protein_tri_fea + RNA_tri_fea
                        #chem_tmp_fea = protein_fea_dict[pro1] + RNA_fea_dict[pro2] 
                    train.append(tmp_fea)
                    #chem_fea.append(chem_tmp_fea)
            else:
                print pro1, pro2
    return np.array(train), label

def calculate_performace(test_num, pred_y,  labels):
    tp =0
    fp = 0
    tn = 0
    fn = 0
    for index in range(test_num):
        if labels[index] ==1:
            if labels[index] == pred_y[index]:
                tp = tp +1
            else:
                fn = fn + 1
        else:
            if labels[index] == pred_y[index]:
                tn = tn +1
            else:
                fp = fp + 1               
            
    acc = float(tp + tn)/test_num
    precision = float(tp)/(tp+ fp)
    sensitivity = float(tp)/ (tp+fn)
    specificity = float(tn)/(tn + fp)
    MCC = float(tp*tn-fp*fn)/(np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return acc, precision, sensitivity, specificity, MCC 

def plot_feature_importance(importance):
    df = pd.DataFrame(importance, columns=['feature', 'fscore'])
    df['fscore'] = df['fscore'] / df['fscore'].sum()
    #pdb.set_trace()
    plt.figure()
    df.plot()
    df.plot(kind='barh', x='feature', y='fscore', legend=False, figsize=(6, 10), fontsize='20')
    plt.title('XGBoost Feature Importance')
    plt.xlabel('relative importance')
    #plt.gcf().savefig('feature_importance_xgb.png')
    plt.show()           

def calculate_performace_without_MCC(test_num, pred_y,  labels):
    tp =0
    fp = 0
    tn = 0
    fn = 0
    for index in range(test_num):
        if labels[index] ==1:
            if labels[index] == pred_y[index]:
                tp = tp +1
            else:
                fn = fn + 1
        else:
            if labels[index] == pred_y[index]:
                tn = tn +1
            else:
                fp = fp + 1               
            
    acc = float(tp + tn)/test_num
    #precision = float(tp)/(tp+ fp)
    sensitivity = float(tp)/ (tp+fn)
    return acc, sensitivity

def get_protein_seq_concolutional_array(seq, motif_len = 4):
    #tar_list = ['0', '1', '2', '3', '4', '5', '6']
    #data = {}
    alpha = '0123456'
    #for key, seq in seq_dict.iteritems():
    row = (len(seq) + 2*motif_len - 2)
    new_array = np.zeros((row, 7))
    for i, val in enumerate(seq):
        if val not in alpha:
            new_array[i] = np.array([0.15]*7)
            continue
        if val == 'N' or i < motif_len or i > len(seq) - motif_len:
            new_array[i] = np.array([0.15]*7)
        else:
            index = alpha.index(val)
        new_array[i][index] = 1
        #data[key] = new_array
    return new_array

def get_RNA_seq_concolutional_array(seq, motif_len = 4):
    '''
        𝑆𝑖,𝑗 = {
    .25 if 𝑠𝑖−𝑚+1 = N or 𝑖 < 𝑚 or 𝑖 > 𝑛 − 𝑚
    1 if 𝑠𝑖−𝑚+1 = 𝑗th base in (A, C, G, T)
    0 otherwise
    '''
    #data = {}
    alpha = 'ACGT'
    #for seq in seqs:
    #for key, seq in seqs.iteritems():
    row = (len(seq) + 2*motif_len - 2)
    new_array = np.zeros((row, 4))
    for i, val in enumerate(seq):
        if val not in 'ACGTN':
            new_array[i] = np.array([0.25]*4)
            continue
        if val == 'N' or i < motif_len or i > len(seq) - motif_len:
            new_array[i] = np.array([0.25]*4)
        else:
            index = alpha.index(val)
            new_array[i][index] = 1
        #data[key] = new_array
    return new_array

def get_RNA_protein_concolutional_array(seq, motif_len = 4):
    '''
        𝑆𝑖,𝑗 = {
    .25 if 𝑠𝑖−𝑚+1 = N or 𝑖 < 𝑚 or 𝑖 > 𝑛 − 𝑚
    1 if 𝑠𝑖−𝑚+1 = 𝑗th base in (A, C, G, T)
    0 otherwise
    '''
    data = {}
    alpha = 'ACGT0123456'
    #for seq in seqs:
    #for key, seq in seqs.iteritems():
    row = (len(seq) + 2*motif_len - 2)
    new_array = np.zeros((row, 11))
    for i, val in enumerate(seq):
        if val not in alpha:
            new_array[i] = np.array([0.25]*11)
        if val == 'N' or i < motif_len or i > len(seq) - motif_len:
            new_array[i] = np.array([0.25]*11)
        else:
            index = alpha.index(val)
            new_array[i][index] = 1
        #data[key] = new_array
    return new_array

def pca_reduce_dimension(group_data, n_components = 50):
    print 'running PCA'
    pca = PCA(n_components=n_components)
    pca.fit(group_data)
    group_data = pca.transform(group_data)
    return group_data

def preprocess_data(X, scaler=None, stand = True):
    if not scaler:
        if stand:
            scaler = StandardScaler()
        else:
            scaler = MinMaxScaler()
        scaler.fit(X)
    X = scaler.transform(X)
    return X, scaler


def preprocess_labels(labels, encoder=None, categorical=True):
    if not encoder:
        encoder = LabelEncoder()
        encoder.fit(labels)
    y = encoder.transform(labels).astype(np.int32)
    if categorical:
        y = np_utils.to_categorical(y)
    return y, encoder

def load_data(path, train=True):
    df = pd.read_csv(path)
    X = df.values.copy()
    if train:
        np.random.shuffle(X)  # https://youtu.be/uyUXoap67N8
        X, labels = X[:, 1:-1].astype(np.float32), X[:, -1]
        return X, labels
    else:
        X, ids = X[:, 1:].astype(np.float32), X[:, 0].astype(str)
        return X, ids

def get_data(dataset):
    if dataset == 'RPI1807':
        X, labels = prepare_RPI1807_feature(graph = False, chem_fea = False) # load_data('train.csv', train=True)
    elif dataset == 'NPInter':
        X, labels = prepare_NPinter_feature(graph = False, chem_fea = False)
    elif dataset == 'RPI2241':
        X, labels = prepare_RPI2241_369_feature('ncRNA-protein/RPI2241_trainingSetFeatures.csv', 'ncRNA-protein/RPI2241.csv', 
                                                  'ncRNA-protein/RPI2241_rna.fa', 'ncRNA-protein/RPI2241_all.txt', 'ncRNA-protein/RPI2241_protein.fa', graph = False)
    elif dataset == 'RPI369':
        X, labels, chem_fea = prepare_RPI2241_369_feature('ncRNA-protein/RPI369_trainingSetFeatures.csv', 'ncRNA-protein/RPI369.csv', 
                                                  'ncRNA-protein/RPI369_rna.fa', 'ncRNA-protein/RPI369_all.txt', 'ncRNA-protein/RPI369_protein.fa', graph = False)
    elif dataset == 'RPI488':
        X, labels = prepare_RPI488_feature()
        
    print X.shape
    #pdb.set_trace()
    X, scaler = preprocess_data(X)
    
    #chem_fea, newscale = preprocess_data(chem_fea)
    
    dims = X.shape[1]
    print(dims, 'dims')
    
    return X, labels


def deep_classifier_keras(dataset = 'RPI-Pred'):
    X, labels, chem_fea = get_data(dataset)
    #y = np.array(labels)
    #X_train, X_test, b_train, b_test = train_test_split(X, y, test_size=0.25, random_state=42)
    #X_test, ids = load_data('test.csv', train=False)
    #X_test, _ = preprocess_data(X_test, scaler)
    y, encoder = preprocess_labels(labels)
    
    num_cross_val = 5
    all_performance = []
    for fold in range(num_cross_val):
        train = []
        test = []
        train = np.array([x for i, x in enumerate(X) if i % num_cross_val != fold])
        test = np.array([x for i, x in enumerate(X) if i % num_cross_val == fold])
        train_label = np.array([x for i, x in enumerate(y) if i % num_cross_val != fold])
        test_label = np.array([x for i, x in enumerate(y) if i % num_cross_val == fold])
        chem_train = np.array([x for i, x in enumerate(chem_fea) if i % num_cross_val != fold])
        chem_test = np.array([x for i, x in enumerate(chem_fea) if i % num_cross_val == fold])
        
        print("Building deep learning model...")
        
        model = Sequential()
        num_hidden = 128
        model.add(Dense(train.shape[1], num_hidden, init='uniform', activation='tanh'))
        #model.add(Activation('relu'))
        #model.add(PReLU((128,)))
        #model.add(BatchNormalization((128,)))
        model.add(Dropout(0.5))
        model.add(Dense(num_hidden, 128, init='uniform', activation='tanh'))
        #model.add(PReLU((128,)))
        #model.add(BatchNormalization((128,)))
        #model.add(Activation('relu'))
        model.add(Dropout(0.5))
        #model.add(Dense(128, 128, init='uniform', activation='tanh'))
        #model.add(Dropout(0.5))
        model.add(Dense(128, train_label.shape[1], init='uniform', activation='softmax'))

        sgd = SGD(lr=0.1, decay=1e-6, momentum=0.9, nesterov=True)
        #sgd = SGD()
        #sgd = RMSprop()
        model.compile(loss='categorical_crossentropy', optimizer=sgd) #"rmsprop")
        #model.fit(np.array(train), np.array(train_label), nb_epoch=20, batch_size=128, validation_split=0.15)
        print("Training model...")
        
        model.fit(train, train_label, nb_epoch=100, batch_size=100, verbose=0 )#, validation_split=0.15)
        #print model.get_weights()
        #model.set_weights(np.array([np.random.uniform(size=k.shape) for k in model.get_weights()]))
        #pdb.set_trace()
        print("Generating submission...")
        #pdb.set_trace()
        proba = model.predict_classes(test)
        
        #for pred, real in zip(proba, test_label):
        #    if real[0] == 1 and pred == 0:
        real_labels = []
        for val in test_label:
            if val[0] == 1:
                real_labels.append(0)
            else:
                real_labels.append(1)
        
        #pdb.set_trace()            
        acc, precision, sensitivity, specificity, MCC = calculate_performace(len(real_labels), proba,  real_labels)
        print acc, precision, sensitivity, specificity, MCC
        all_performance.append([acc, precision, sensitivity, specificity, MCC])
    print 'mean performance'
    print np.mean(np.array(all_performance), axis=0)

def build_deep_classical_autoencoder(autoencoder, input_dim, hidden_dim, activation, weight_reg = None, activity_reg = None):
    encoder = containers.Sequential([Dense(input_dim, hidden_dim, activation=activation, W_regularizer=weight_reg, activity_regularizer=activity_reg),
                            Dense(hidden_dim, hidden_dim/2, activation=activation)])
    decoder = containers.Sequential([Dense(hidden_dim/2, hidden_dim, activation=activation), 
                                     Dense(hidden_dim, input_dim, activation=activation, W_regularizer=weight_reg, activity_regularizer=activity_reg)])
    autoencoder.add(AutoEncoder(encoder=encoder, decoder=decoder, output_reconstruction=False))
    return autoencoder


def deep_autoencoder(dataset = 'RPI-Pred'):
    X, labels, chem_fea = get_data(dataset)
    y, encoder = preprocess_labels(labels)
    
    num_cross_val = 5
    batch_size  = 50
    nb_epoch = 100
    all_performance = []
    activation = 'linear' #'linear' #'relu, softmax, tanh'
    for fold in range(num_cross_val):
        train = []
        test = []
        train = np.array([x for i, x in enumerate(X) if i % num_cross_val != fold])
        test = np.array([x for i, x in enumerate(X) if i % num_cross_val == fold])
        train_label = np.array([x for i, x in enumerate(y) if i % num_cross_val != fold])
        test_label = np.array([x for i, x in enumerate(y) if i % num_cross_val == fold])
        chem_train = np.array([x for i, x in enumerate(chem_fea) if i % num_cross_val != fold])
        chem_test = np.array([x for i, x in enumerate(chem_fea) if i % num_cross_val == fold])
                
        autoencoder = Sequential()
        autoencoder = build_deep_classical_autoencoder(autoencoder, train.shape[1], 256, activation)
        autoencoder.get_config(verbose=0)
        #norm_m0 = normalization.BatchNormalization((599,))
        #autoencoder.add(norm_m0)
        autoencoder.add(Dropout(0.5))
        sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
        autoencoder.compile(loss='mean_squared_error', optimizer= sgd) #'adam')
        # Do NOT use validation data with return output_reconstruction=True
        autoencoder.fit(train, train, batch_size=batch_size, nb_epoch=nb_epoch, show_accuracy=False, verbose=0)
    
        prefilter_train = autoencoder.predict(train, verbose=0)
        prefilter_test = autoencoder.predict(test, verbose=0)
        prefilter_train = np.concatenate((prefilter_train, chem_train), axis = 1)
        prefilter_test = np.concatenate((prefilter_test, chem_test), axis = 1)
        
        print 'using random forest'
        train_label_new = []
        for val in train_label:
            if val[0] == 1:
                train_label_new.append(0)
            else:
                train_label_new.append(1)
        #parameters = {'kernel': ['linear', 'rbf'], 'C': [1, 2, 3, 4, 5, 6, 10], 'gamma': [0.5,1,2,4, 6, 8]}
        #svr = svm.SVC(probability = True)
        #clf = grid_search.GridSearchCV(svr, parameters, cv=3)        
        clf = RandomForestClassifier(n_estimators=50)
        clf.fit(prefilter_train, train_label_new)
        y_pred = clf.predict(prefilter_test)
        real_labels = []
        for val in test_label:
            if val[0] == 1:
                real_labels.append(0)
            else:
                real_labels.append(1)
        acc, precision, sensitivity, specificity, MCC = calculate_performace(len(real_labels), y_pred,  real_labels)
        print acc, precision, sensitivity, specificity, MCC
        all_performance.append([acc, precision, sensitivity, specificity, MCC])

        prefilter_train = []
        prefilter_test = []    
    print 'mean performance'
    print np.mean(np.array(all_performance), axis=0)

def construct_one_layer_network(X_train, X_test, input_dim, output_dim, activation = 'linear', batch_size = 100, nb_epoch = 100):
    print 'constructing one-layer network'
    autoencoder = Sequential()
    autoencoder = build_deep_classical_autoencoder(autoencoder, input_dim, output_dim, activation)
    autoencoder.get_config(verbose=0)
    #norm_m0 = normalization.BatchNormalization((599,))
    #autoencoder.add(norm_m0)
    autoencoder.add(Dropout(0.5))
    sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    autoencoder.compile(loss='mean_squared_error', optimizer= sgd) #'adam')
    # Do NOT use validation data with return output_reconstruction=True
    autoencoder.fit(X_train, X_train, batch_size=batch_size, nb_epoch=nb_epoch, show_accuracy=False, verbose=0)
    output_reconstruction=False
    first_train = autoencoder.predict(X_train, verbose=0)
    first_test = autoencoder.predict(X_test, verbose=0)
    
    return autoencoder, first_train, first_test

def autoencoder_fine_tuning(encoders, X_train, Y_train, X_test, batch_size, nb_epoch):
    print 'fine tunning'
    model = Sequential()
    for encoder in encoders:
        model.add(encoder.layers[0].encoder)
    model.add(Dense(128, 2, activation='softmax'))

    model.compile(loss='categorical_crossentropy', optimizer='rmsprop')
    model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=nb_epoch, show_accuracy=False, verbose=0)
    #pdb.set_trace()
    X_train_tmp = np.copy(X_train)
    X_test_tmp = np.copy(X_test)   
    model2 = Sequential() 
    for ae in model.layers[:-1]:
        model2.add(ae)
        
    model2.compile(loss='mean_squared_error', optimizer='adam')
    X_train_tmp = model2.predict(X_train_tmp)
    print X_train_tmp.shape
    X_test_tmp = model2.predict(X_test_tmp)
    
    return X_train_tmp, X_test_tmp        

def get_preds( score1, score2, score3, weights):
    new_score  = [weights[0]*val1 + weights[1]* val2 + weights[2]* val3  for val1, val2, val3 in zip(score1, score2, score3)]
    return new_score

def transfer_label_from_prob(proba):
    label = [1 if val>=0.5 else 0 for val in proba]
    return label

def get_blend_data(j, clf, skf, X_test, X_dev, Y_dev, blend_train, blend_test):
        print 'Training classifier [%s]' % (j)
        blend_test_j = np.zeros((X_test.shape[0], len(skf))) # Number of testing data x Number of folds , we will take the mean of the predictions later
        for i, (train_index, cv_index) in enumerate(skf):
            print 'Fold [%s]' % (i)
            
            # This is the training and validation set
            X_train = X_dev[train_index]
            Y_train = Y_dev[train_index]
            X_cv = X_dev[cv_index]
            Y_cv = Y_dev[cv_index]
            
            clf.fit(X_train, Y_train)
            
            # This output will be the basis for our blended classifier to train against,
            # which is also the output of our classifiers
            #blend_train[cv_index, j] = clf.predict(X_cv)
            #blend_test_j[:, i] = clf.predict(X_test)
            blend_train[cv_index, j] = clf.predict_proba(X_cv)[:,1]
            blend_test_j[:, i] = clf.predict_proba(X_test)[:,1]
        # Take the mean of the predictions of the cross validation set
        blend_test[:, j] = blend_test_j.mean(1)
    
        print 'Y_dev.shape = %s' % (Y_dev.shape)

def multiple_autoencoder_extract_feature(dataset = 'RPI-Pred'):
    X, labels = get_data(dataset)
    #X = pca_reduce_dimension(X, n_components = 300)
    y, encoder = preprocess_labels(labels)
    num_cross_val = 5
    batch_size  = 50
    nb_epoch = 100
    all_performance = []
    all_performance_rf = []
    all_performance_ensemb = []
    all_performance_ae = []
    all_performance_rf_seq = []
    all_performance_chem = []
    all_performance_blend = []
    activation = 'linear' #'linear' #'relu, softmax, tanh'
    for fold in range(num_cross_val):
        train = []
        test = []
        train = np.array([x for i, x in enumerate(X) if i % num_cross_val != fold])
        test = np.array([x for i, x in enumerate(X) if i % num_cross_val == fold])
        train_label = np.array([x for i, x in enumerate(y) if i % num_cross_val != fold])
        test_label = np.array([x for i, x in enumerate(y) if i % num_cross_val == fold])
        #chem_train = np.array([x for i, x in enumerate(chem_fea) if i % num_cross_val != fold])
        #chem_test = np.array([x for i, x in enumerate(chem_fea) if i % num_cross_val == fold])

        real_labels = []
        for val in test_label:
            if val[0] == 1:
                real_labels.append(0)
            else:
                real_labels.append(1)

        train_label_new = []
        for val in train_label:
            if val[0] == 1:
                train_label_new.append(0)
            else:
                train_label_new.append(1)
        blend_train = np.zeros((train.shape[0], 3)) # Number of training data x Number of classifiers
        blend_test = np.zeros((test.shape[0], 3)) # Number of testing data x Number of classifiers 
        skf = list(StratifiedKFold(train_label_new, 3))  
        #pdb.set_trace()     
        class_index = 0                   
        encoders = multiple_layer_autoencoder(train, test, activation = 'sigmoid', batch_size = 100, nb_epoch = 100, last_dim = 128)
        prefilter_train = np.copy(train)
        prefilter_test = np.copy(test) 
        for ae in encoders:
            prefilter_train = ae.predict(prefilter_train)
            print prefilter_train.shape
            prefilter_test = ae.predict(prefilter_test)
            
        clf = RandomForestClassifier(n_estimators=50)
        clf.fit(prefilter_train, train_label_new)
        y_pred_prob = clf.predict_proba(prefilter_test)[:,1]
        #pdb.set_trace()
        y_pred = transfer_label_from_prob(y_pred_prob)
        acc, precision, sensitivity, specificity, MCC = calculate_performace(len(real_labels), y_pred,  real_labels)
        print acc, precision, sensitivity, specificity, MCC
        all_performance_ae.append([acc, precision, sensitivity, specificity, MCC])
        print '---' * 50
        
        get_blend_data(class_index, RandomForestClassifier(n_estimators=50), skf, prefilter_test, prefilter_train, np.array(train_label_new), blend_train, blend_test)
        
        
        prefilter_train, prefilter_test = autoencoder_fine_tuning(encoders, train, train_label, test, 100, 100)
        #prefilter_train, new_scaler = preprocess_data(prefilter_train)
        #prefilter_test, new_scaler = preprocess_data(prefilter_test, scaler = new_scaler)
        #prefilter_train = np.concatenate((prefilter_train, chem_train), axis = 1)
        #prefilter_test = np.concatenate((prefilter_test, chem_test), axis = 1)        
        print 'using random forest after sequence autoencoder'
        class_index = class_index + 1
        #parameters = {'kernel': ['linear', 'rbf'], 'C': [1, 2, 3, 4, 5, 6, 10], 'gamma': [0.5,1,2,4, 6, 8]}
        #svr = svm.SVC(probability = True)
        #clf = grid_search.GridSearchCV(svr, parameters, cv=3)        
        clf = RandomForestClassifier(n_estimators=50)
        clf.fit(prefilter_train, train_label_new)
        y_pred_prob = clf.predict_proba(prefilter_test)[:,1]
        #pdb.set_trace()
        y_pred = transfer_label_from_prob(y_pred_prob)
        acc, precision, sensitivity, specificity, MCC = calculate_performace(len(real_labels), y_pred,  real_labels)
        print acc, precision, sensitivity, specificity, MCC
        all_performance.append([acc, precision, sensitivity, specificity, MCC])
        print '---' * 50
        
        get_blend_data(class_index, RandomForestClassifier(n_estimators=50), skf, prefilter_test, prefilter_train, np.array(train_label_new), blend_train, blend_test)
        
        print 'using RF using only sequence feature'
        class_index = class_index + 1
        clf = RandomForestClassifier(n_estimators=50)
        clf.fit(train, train_label_new)
        y_pred_rf_prob = clf.predict_proba(test)[:,1]
        y_pred_rf = transfer_label_from_prob(y_pred_rf_prob)
        acc, precision, sensitivity, specificity, MCC = calculate_performace(len(real_labels), y_pred_rf,  real_labels)
        print acc, precision, sensitivity, specificity, MCC
        all_performance_rf_seq.append([acc, precision, sensitivity, specificity, MCC]) 
        print '---' * 50
        get_blend_data(class_index, RandomForestClassifier(n_estimators=50), skf, test, train, np.array(train_label_new), blend_train, blend_test)
        
        #roc = metrics.roc_auc_score(y_valid, valid_preds)
        # Start blending!
        bclf = LogisticRegression()
        bclf.fit(blend_train, train_label_new)
        Y_test_predict = bclf.predict(blend_test)
        print 'blend result'
        acc, precision, sensitivity, specificity, MCC = calculate_performace(len(real_labels), Y_test_predict,  real_labels)
        print acc, precision, sensitivity, specificity, MCC   
        all_performance_blend.append([acc, precision, sensitivity, specificity, MCC])     
        print '---' * 50
        
        '''
        print 'ensemble deep learning and rf'
        ensemb_prob = get_preds( y_pred_prob, y_pred_rf_prob, ae_y_pred, [0.3, 0.30, 0.4])       
        ensemb_label = transfer_label_from_prob(ensemb_prob)
        acc, precision, sensitivity, specificity, MCC = calculate_performace(len(real_labels), ensemb_label,  real_labels)
        print acc, precision, sensitivity, specificity, MCC
        all_performance_ensemb.append([acc, precision, sensitivity, specificity, MCC]) 
        print '---' * 50
        '''
    print 'in summary'
    print 'mean performance of chem autoencoder without fine tunning'
    print np.mean(np.array(all_performance_ae), axis=0)  
    print '---' * 50
    print 'mean performance of sequence autoencoder'
    print np.mean(np.array(all_performance), axis=0)
    print '---' * 50   
    print 'mean performance of only chem using RF'
    print np.mean(np.array(all_performance_rf), axis=0)
    print '---' * 50     
    print 'mean performance of blend fusion'
    print np.mean(np.array(all_performance_blend), axis=0) 
    print '---' * 50 
           
def random_forest_classify(dataset = 'RPI-Pred', SVM = False):
    X, labels, chem_fea = get_data(dataset)
    #X, scaler = preprocess_data(chem_fea)
    #y, encoder = preprocess_labels(labels)
    y = np.array(labels)
    #X_train, X_test, b_train, b_test = train_test_split(X, y, test_size=0.25, random_state=42)
    #X_test, ids = load_data('test.csv', train=False)
    #X_test, _ = preprocess_data(X_test, scaler)
    
    dims = X.shape[1]
    print(dims, 'dims')
    num_cross_val = 5
    all_performance = []
    for fold in range(num_cross_val):
        train = []
        test = []
        train = [x for i, x in enumerate(X) if i % num_cross_val != fold]
        test = [x for i, x in enumerate(X) if i % num_cross_val == fold]
        train_label = [x for i, x in enumerate(y) if i % num_cross_val != fold]
        test_label = [x for i, x in enumerate(y) if i % num_cross_val == fold]
        if SVM:
            parameters = {'kernel': ['linear', 'rbf'], 'C': [1, 2, 3, 4, 5, 6, 10], 'gamma': [0.5,1,2,4, 6, 8]}
            svr = svm.SVC(probability = True)
            clf = grid_search.GridSearchCV(svr, parameters, cv=3)  
        else:  
            print 'using random forest'
            clf = RandomForestClassifier(n_estimators=50)
        clf.fit(np.array(train), train_label)
        y_pred = clf.predict(np.array(test))
        acc, precision, sensitivity, specificity, MCC = calculate_performace(len(test_label), y_pred,  test_label)
        #print acc, precision, sensitivity, specificity, MCC
        all_performance.append([acc, precision, sensitivity, specificity, MCC])
    print 'mean performance'
    print np.mean(np.array(all_performance), axis=0)

def indep_validation():
    X, labels = prepare_RPI1807_feature(graph = False)
    print X.shape
    test, test_label = prepare_NPinter_feature(extract_only_posi = False, graph = False)
    print test.shape
    X, scaler = preprocess_data(X)
    test, _ = preprocess_data(test, scaler=scaler)
    print 'using random forest'
    clf = RandomForestClassifier(n_estimators=50)
    #calibrated_clf = CalibratedClassifierCV(clf, method='isotonic', cv=5)
    #calibrated_clf.fit(np.array(X), np.array(labels))
    #y_pred = calibrated_clf.predict(test)
    
    clf.fit(np.array(X), np.array(labels))
    #pdb.set_trace()
    y_pred = clf.predict(test)
    #acc, sensitivity = calculate_performace_without_MCC(len(test_label), y_pred,  test_label)
    #print acc, sensitivity  
    acc, precision, sensitivity, specificity, MCC = calculate_performace(len(test_label), y_pred,  test_label)
    print acc, precision, sensitivity, specificity, MCC  


def get_data_deepmind(dataset, deepmind =False, seperate = True, chem_fea = False, extract_only_posi = False, indep_test = False):
    if dataset == 'RPI1807':
        X, labels = prepare_RPI1807_feature(graph = False, deepmind = deepmind, seperate = seperate, chem_fea = chem_fea) # load_data('train.csv', train=True)
    elif dataset == 'NPInter':
        X, labels = prepare_NPinter_feature(graph = False, deepmind = deepmind, seperate = seperate, chem_fea = chem_fea, extract_only_posi = extract_only_posi)
    elif dataset == 'RPI2241':
        X, labels = prepare_RPI2241_369_feature('ncRNA-protein/RPI2241_rna.fa', 'ncRNA-protein/RPI2241_all.txt', 
                                                  'ncRNA-protein/RPI2241_protein.fa', graph = False, deepmind = deepmind, seperate = seperate, chem_fea = chem_fea)
    elif dataset == 'RPI369':
        X, labels = prepare_RPI2241_369_feature('ncRNA-protein/RPI369_rna.fa', 'ncRNA-protein/RPI369_all.txt', 
                                                  'ncRNA-protein/RPI369_protein.fa', graph = False, deepmind = deepmind, seperate = seperate, chem_fea = chem_fea)
    elif dataset == 'RPI488':
        X, labels = prepare_RPI488_feature(deepmind = deepmind, seperate = seperate, chem_fea = chem_fea)
    elif dataset == 'RPIntDB':
        X, labels = prepare_RPIntDB_feature(seperate = seperate)
    elif dataset == 'RPI13254':
        X, labels = prepare_RPI13254_feature(seperate = seperate, extract_only_posi = extract_only_posi, indep_test = indep_test) 
    elif dataset == 'RPI367':
        X, labels = prepare_RPI367_feature(seperate = seperate, extract_only_posi = extract_only_posi)            
    return X, labels
    #return chem_fea, labels

def transfer_array_format(data):
    formated_matrix1 = []
    formated_matrix2 = []
    #pdb.set_trace()
    #pdb.set_trace()
    for val in data:
        #formated_matrix1.append(np.array([val[0]]))
        formated_matrix1.append(val[0])
        formated_matrix2.append(val[1])
        #formated_matrix1[0] = np.array([val[0]])
        #formated_matrix2.append(np.array([val[1]]))
        #formated_matrix2[0] = val[1]      
    
    return np.array(formated_matrix1), np.array(formated_matrix2)  

def get_rnn_fea(train, sec_num_hidden = 128):
    model = Sequential()
    num_hidden = 128
    #sec_num_hidden = 128
    #max_features = train.shape[1]
    #model.add(Embedding(max_features, 256))
    #model.add(LSTM(256, num_hidden, activation='sigmoid', inner_activation='hard_sigmoid'))
    model.add(Dense(train.shape[1], num_hidden, init='uniform', activation='sigmoid'))
    #model.add(PReLU((num_hidden,)))
    #model.add(BatchNormalization((num_hidden,)))
    model.add(Dropout(0.5))
    model.add(Dense(num_hidden, sec_num_hidden, init='uniform', activation='sigmoid'))
    model.add(PReLU((sec_num_hidden,)))
    model.add(BatchNormalization((sec_num_hidden,)))
    #model.add(Activation('relu'))
    model.add(Dropout(0.5))
    return model

def set_rnn_fea(train, sec_num_hidden = 128, weights1 = None):
    model = Sequential()
    num_hidden = 256
    #sec_num_hidden = 128
    #max_features = train.shape[1]
    #model.add(Embedding(max_features, 256))
    #model.add(LSTM(256, num_hidden, activation='sigmoid', inner_activation='hard_sigmoid'))
    model.add(Dense(train.shape[1], num_hidden, init='uniform', activation='sigmoid'))
    #model.add(PReLU((num_hidden,)))
    #model.add(BatchNormalization((num_hidden,)))
    model.add(Dropout(0.5))
    model.add(Dense(num_hidden, sec_num_hidden, init='uniform', activation='sigmoid'))
    model.add(PReLU((sec_num_hidden,)))
    model.add(BatchNormalization((sec_num_hidden,)))
    #model.add(Activation('relu'))
    model.add(Dropout(0.5))
    return model

def set_seperate_network(old_model, X_train1, X_train2):
    left_hid = 64
    right_hid = 64
    weight =  old_model.layers[0].get_weights()
    left = set_rnn_fea(X_train1, sec_num_hidden = left_hid)
    left.set_weights(weight[:7])
    right = set_rnn_fea(X_train2, sec_num_hidden = right_hid)
    right.set_weights(weight[7:])
    
    model = Sequential()
    model.add(Merge([left, right], mode='concat'))
    return model

def multiple_layer_autoencoder(X_train, X_test, activation = 'linear', batch_size = 100, nb_epoch = 100, last_dim = 64):
    nb_hidden_layers = [X_train.shape[1], 256, 128, last_dim]
    X_train_tmp = np.copy(X_train)
    #X_test_tmp = np.copy(X_test)
    encoders = []
    for i, (n_in, n_out) in enumerate(zip(nb_hidden_layers[:-1], nb_hidden_layers[1:]), start=1):
        print('Training the layer {}: Input {} -> Output {}'.format(i, n_in, n_out))
        # Create AE and training
        ae = Sequential()
        encoder = containers.Sequential([Dense(n_in, n_out, activation=activation)])
        decoder = containers.Sequential([Dense(n_out, n_in, activation=activation)])
        ae.add(AutoEncoder(encoder=encoder, decoder=decoder,
                           output_reconstruction=False))
        ae.add(Dropout(0.5))
        #sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
        ae.compile(loss='mean_squared_error', optimizer='adam')#'rmsprop')
        ae.fit(X_train_tmp, X_train_tmp, batch_size=batch_size, nb_epoch=nb_epoch, show_accuracy=False, verbose=0)
        # Store trainined weight and update training data
        #encoders.append(ae.layers[0].encoder)
        encoders.append(ae)
        X_train_tmp = ae.predict(X_train_tmp)
        print X_train_tmp.shape
        #X_test_tmp = ae.predict(X_test_tmp)
        
    #return encoders, X_train_tmp, X_test_tmp
    return encoders

def autoencoder_two_subnetwork_fine_tuning(X_train1, X_train2, Y_train, X_test1, X_test2, Y_test = None, batch_size =100, nb_epoch = 100):
    print 'autoencode learning'
    last_dim = 64
    encoders1 = multiple_layer_autoencoder(X_train1, X_test1, activation = 'sigmoid', batch_size = batch_size, nb_epoch = nb_epoch, last_dim = last_dim)
    encoders2 = multiple_layer_autoencoder(X_train2, X_test2, activation = 'sigmoid', batch_size = batch_size, nb_epoch = nb_epoch, last_dim = last_dim)
    #pdb.set_trace()
    
    X_train1_tmp_bef = np.copy(X_train1)
    X_test1_tmp_bef = np.copy(X_test1) 
    for ae in encoders1:
        X_train1_tmp_bef = ae.predict(X_train1_tmp_bef)
        print X_train1_tmp_bef.shape
        X_test1_tmp_bef = ae.predict(X_test1_tmp_bef)
    
    X_train2_tmp_bef = np.copy(X_train2)
    X_test2_tmp_bef = np.copy(X_test2) 
    for ae in encoders2:
        X_train2_tmp_bef = ae.predict(X_train2_tmp_bef)
        print X_train2_tmp_bef.shape
        X_test2_tmp_bef = ae.predict(X_test2_tmp_bef)
        
    prefilter_train_bef = np.concatenate((X_train1_tmp_bef, X_train2_tmp_bef), axis = 1)
    prefilter_test_bef = np.concatenate((X_test1_tmp_bef, X_test2_tmp_bef), axis = 1)
        
    print 'fine tunning'
    print 'number of layers:', len(encoders1)
    sec_num_hidden = last_dim
    model1 = Sequential()
    ind = 0
    for encoder in encoders1:
        model1.add(encoder.layers[0].encoder)
        if ind != len(encoders1)  - 1 :
            model1.add(Dropout(0.5)) 
            ind = ind + 1
    model1.add(PReLU((sec_num_hidden,)))
    model1.add(BatchNormalization((sec_num_hidden,)))
    model1.add(Dropout(0.5))
    

    model2 = Sequential()
    ind = 0
    for encoder in encoders2:
        model2.add(encoder.layers[0].encoder)
        if ind != len(encoders2)  - 1 :
            model2.add(Dropout(0.5)) 
            ind = ind + 1
    model2.add(PReLU((sec_num_hidden,)))
    model2.add(BatchNormalization((sec_num_hidden,)))   
    model2.add(Dropout(0.5))     
         
    model = Sequential()
    model.add(Merge([model1, model2], mode='concat'))
    total_hid = sec_num_hidden + sec_num_hidden
    
    model.add(Dense(total_hid, 2))
    model.add(Dropout(0.5))
    model.add(Activation('softmax'))
    #model.get_config(verbose=0)
    sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy', optimizer=sgd) #'rmsprop')
    model.fit([X_train1, X_train2], Y_train, batch_size=100, nb_epoch=100, verbose=0)
    #config = autoencoder.get_config(verbose=1)
    #autoencoder = model_from_config(config)
    #pdb.set_trace()
    X_train1_tmp = np.copy(X_train1)
    X_test1_tmp = np.copy(X_test1)  
    ae=model.layers[0].layers[0]  
    ae.compile(loss='mean_squared_error', optimizer='adam')
    X_train1_tmp = ae.predict(X_train1_tmp)
    X_test1_tmp = ae.predict(X_test1_tmp)

    X_train2_tmp = np.copy(X_train2)
    X_test2_tmp = np.copy(X_test2)  
    ae=model.layers[0].layers[1]  
    ae.compile(loss='mean_squared_error', optimizer='adam')
    X_train2_tmp = ae.predict(X_train2_tmp)
    X_test2_tmp = ae.predict(X_test2_tmp)
    
    prefilter_train = np.concatenate((X_train1_tmp, X_train2_tmp), axis = 1)
    prefilter_test = np.concatenate((X_test1_tmp, X_test2_tmp), axis = 1)
    #return X_train1_tmp, X_test1_tmp, X_train2_tmp, X_test2_tmp, model
    return prefilter_train, prefilter_test, prefilter_train_bef, prefilter_test_bef
    #return model

def merge_seperate_network(X_train1, X_train2, Y_train):
    left_hid = 128
    right_hid = 64
    left = get_rnn_fea(X_train1, sec_num_hidden = left_hid)
    right = get_rnn_fea(X_train2, sec_num_hidden = right_hid)
    
    model = Sequential()
    model.add(Merge([left, right], mode='concat'))
    total_hid = left_hid + right_hid
    
    model.add(Dense(total_hid, 2))
    model.add(Dropout(0.3))
    model.add(Activation('softmax'))
    
    sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy', optimizer=sgd) #'rmsprop')
    
    model.fit([X_train1, X_train2], Y_train, batch_size=100, nb_epoch=100, verbose=0)
    
    return model

def merge_seperate_network_with_chem(X_train1, X_train2, chem_train1, chem_train2, Y_train):
    left_hid = 128
    right_hid = 64
    chem_left_hid = 128
    chem_right_hid = 20
    left = get_rnn_fea(X_train1, sec_num_hidden = left_hid)
    right = get_rnn_fea(X_train2, sec_num_hidden = right_hid)
    chem_left = get_rnn_fea(chem_train1, sec_num_hidden = chem_left_hid)
    chem_right = get_rnn_fea(chem_train2, sec_num_hidden = chem_right_hid)
    
    model = Sequential()
    model.add(Merge([left, right, chem_left, chem_right], mode='concat'))
    total_hid = left_hid + right_hid + chem_left_hid + chem_right_hid
    
    model.add(Dense(total_hid, 2))
    model.add(Activation('softmax'))
    
    sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy', optimizer=sgd) #'rmsprop')
    
    model.fit([X_train1, X_train2, chem_train1, chem_train2], Y_train, batch_size=100, nb_epoch=100, verbose=0)
    
    return model
    
def shuffle_two_array(X, Y, labels):
    new_array = range(len(labels))
    np.random.shuffle(new_array)
    #pdb.set_trace()
    new_X = []
    new_Y = []
    new_labels = []
    for val in new_array:
        new_X.append(X[val])
        new_Y.append(Y[val])
        new_labels.append(labels[val])
    #X, labels = new_array[:, 1:-1].astype(np.float32), new_array[:, -1]
    
    return np.array(new_X), np.array(new_Y), np.array(new_labels)

def plot_roc_curve(labels, probality, legend_text, auc_tag = True):
    #fpr2, tpr2, thresholds = roc_curve(labels, pred_y)
    fpr, tpr, thresholds = roc_curve(labels, probality) #probas_[:, 1])
    roc_auc = auc(fpr, tpr)
    if auc_tag:
        rects1 = plt.plot(fpr, tpr, label=legend_text +' (AUC=%6.3f) ' %roc_auc)
    else:
        rects1 = plt.plot(fpr, tpr, label=legend_text )

def IPMiner(dataset):
    X, labels = get_data_deepmind(dataset, seperate = True)
    #pdb.set_trace()
    #X, labels = shuffle_two_array(X, np.array(labels))
    X_data1, X_data2 = transfer_array_format(X)
    print X_data1.shape, X_data2.shape
    X_data1, scaler1 = preprocess_data(X_data1)
    X_data2, scaler2 = preprocess_data(X_data2)
    #X_data1 = X_data1 +np.random.normal(0, 1.0, size=(X_data1.shape))
    #X_data2 = X_data2 +np.random.normal(0, 1.0, size=X_data2.shape)
    y, encoder = preprocess_labels(labels)
    
    #param = {'bst:max_depth':2, 'bst:eta':1, 'silent':1, 'objective':'binary:logistic' }
    #param['nthread'] = 4
    #plst = param.items()
    
    num_cross_val = 5
    all_performance = []
    all_performance_rf = []
    all_performance_bef = []
    all_performance_blend = []
    all_labels = []
    all_prob = {}
    num_classifier = 3
    all_prob[0] = []
    all_prob[1] = []
    all_prob[2] = []
    all_prob[3] = []
    all_averrage = []
    for fold in range(num_cross_val):
        train1 = np.array([x for i, x in enumerate(X_data1) if i % num_cross_val != fold])
        test1 = np.array([x for i, x in enumerate(X_data1) if i % num_cross_val == fold])
        train2 = np.array([x for i, x in enumerate(X_data2) if i % num_cross_val != fold])
        test2 = np.array([x for i, x in enumerate(X_data2) if i % num_cross_val == fold])
        train_label = np.array([x for i, x in enumerate(y) if i % num_cross_val != fold])
        test_label = np.array([x for i, x in enumerate(y) if i % num_cross_val == fold])
  
          
        real_labels = []
        for val in test_label:
            if val[0] == 1:
                real_labels.append(0)
            else:
                real_labels.append(1)

        train_label_new = []
        for val in train_label:
            if val[0] == 1:
                train_label_new.append(0)
            else:
                train_label_new.append(1)
        
        blend_train = np.zeros((train1.shape[0], num_classifier)) # Number of training data x Number of classifiers
        blend_test = np.zeros((test1.shape[0], num_classifier)) # Number of testing data x Number of classifiers 
        skf = list(StratifiedKFold(train_label_new, num_classifier))  
        class_index = 0
        prefilter_train, prefilter_test, prefilter_train_bef, prefilter_test_bef = autoencoder_two_subnetwork_fine_tuning(train1, train2, train_label, test1, test2, test_label)
        #X_train1_tmp, X_test1_tmp, X_train2_tmp, X_test2_tmp, model = autoencoder_two_subnetwork_fine_tuning(train1, train2, train_label, test1, test2, test_label)
        #model = autoencoder_two_subnetwork_fine_tuning(train1, train2, train_label, test1, test2, test_label)
        #model = merge_seperate_network(train1, train2, train_label)
        #proba = model.predict_proba([test1, test2])[:1]
        
        
        real_labels = []
        for val in test_label:
            if val[0] == 1:
                real_labels.append(0)
            else:
                real_labels.append(1)
                
        all_labels = all_labels + real_labels
        #prefilter_train, new_scaler = preprocess_data(prefilter_train, stand =False)
        #prefilter_test, new_scaler = preprocess_data(prefilter_test, scaler = new_scaler, stand = False)
        '''
        prefilter_train1 = xgb.DMatrix( prefilter_train, label=train_label_new)
        evallist  = [(prefilter_train1, 'train')]
        num_round = 10
        clf = xgb.train( plst, prefilter_train1, num_round, evallist )
        prefilter_test1 = xgb.DMatrix( prefilter_test)
        ae_y_pred_prob = clf.predict(prefilter_test1)
        '''
        tmp_aver = [0] * len(real_labels)
        print 'deep autoencoder'
        clf = RandomForestClassifier(n_estimators=50)
        clf.fit(prefilter_train, train_label_new)
        ae_y_pred_prob = clf.predict_proba(prefilter_test)[:,1]
        all_prob[class_index] = all_prob[class_index] + [val for val in ae_y_pred_prob]
        tmp_aver = [val1 + val2/3 for val1, val2 in zip(ae_y_pred_prob, tmp_aver)]
        proba = transfer_label_from_prob(ae_y_pred_prob)
        #pdb.set_trace()            
        acc, precision, sensitivity, specificity, MCC = calculate_performace(len(real_labels), proba,  real_labels)
        print acc, precision, sensitivity, specificity, MCC
        all_performance.append([acc, precision, sensitivity, specificity, MCC])
        get_blend_data(class_index, RandomForestClassifier(n_estimators=50), skf, prefilter_test, prefilter_train, np.array(train_label_new), blend_train, blend_test)
        
        print 'deep autoencoder without fine tunning'
        class_index = class_index + 1
        clf = RandomForestClassifier(n_estimators=50)
        clf.fit(prefilter_train_bef, train_label_new)
        ae_y_pred_prob_bef = clf.predict_proba(prefilter_test_bef)[:,1]
        all_prob[class_index] = all_prob[class_index] + [val for val in ae_y_pred_prob_bef]
        tmp_aver = [val1 + val2/3 for val1, val2 in zip(ae_y_pred_prob_bef, tmp_aver)]
        proba = transfer_label_from_prob(ae_y_pred_prob_bef)
        #pdb.set_trace()            
        acc, precision, sensitivity, specificity, MCC = calculate_performace(len(real_labels), proba,  real_labels)
        print acc, precision, sensitivity, specificity, MCC
        all_performance_bef.append([acc, precision, sensitivity, specificity, MCC])
        get_blend_data(class_index, RandomForestClassifier(n_estimators=50), skf, prefilter_test_bef, prefilter_train_bef, np.array(train_label_new), blend_train, blend_test)
        
        print 'random forest using raw feature'
        class_index = class_index + 1
        prefilter_train = np.concatenate((train1, train2), axis = 1)
        prefilter_test = np.concatenate((test1, test2), axis = 1)
        
        clf = RandomForestClassifier(n_estimators=50)
        clf.fit(prefilter_train, train_label_new)
        ae_y_pred_prob = clf.predict_proba(prefilter_test)[:,1]
        all_prob[class_index] = all_prob[class_index] + [val for val in ae_y_pred_prob]
        tmp_aver = [val1 + val2/3 for val1, val2 in zip(ae_y_pred_prob, tmp_aver)]
        proba = transfer_label_from_prob(ae_y_pred_prob)
                 
        acc, precision, sensitivity, specificity, MCC = calculate_performace(len(real_labels), proba,  real_labels)
        print acc, precision, sensitivity, specificity, MCC
        all_performance_rf.append([acc, precision, sensitivity, specificity, MCC])
        get_blend_data(class_index, RandomForestClassifier(n_estimators=50), skf, prefilter_test, prefilter_train, np.array(train_label_new), blend_train, blend_test)
        
        class_index = class_index + 1
        bclf = LogisticRegression()
        bclf.fit(blend_train, train_label_new)
        #print bclf.coef_, bclf.intercept_
        stack_y_prob = bclf.predict_proba(blend_test)[:,1]
        all_prob[class_index] = all_prob[class_index] + [val for val in stack_y_prob]
        Y_test_predict = bclf.predict(blend_test)
        print 'stacked ensembling'
        acc, precision, sensitivity, specificity, MCC = calculate_performace(len(real_labels), Y_test_predict,  real_labels)
        print acc, precision, sensitivity, specificity, MCC   
        all_performance_blend.append([acc, precision, sensitivity, specificity, MCC])     
        print '---' * 50
        all_averrage = all_averrage + tmp_aver
        
    print 'mean performance of deep autoencoder'
    print np.mean(np.array(all_performance), axis=0)
    print '---' * 50 
    print 'mean performance of deep autoencoder without fine tunning'
    print np.mean(np.array(all_performance_bef), axis=0)
    print '---' * 50 
    print 'mean performance of random forest using raw feature'
    print np.mean(np.array(all_performance_rf), axis=0)
    print '---' * 50    
    print 'mean performance of stacked ensembling'
    print np.mean(np.array(all_performance_blend), axis=0)
    print '---' * 50
    
def predict_new_samples(RNA_file, protein_file):
    train, trainlabel = get_data_deepmind('RPI488', seperate = True)
    #test, testlabel = get_data_deepmind(datatype, seperate = True, extract_only_posi = extract_only_posi, indep_test = test_indep_dataset)
    test, pairs = prepare_complex_feature(RNA_file, protein_file, seperate = True)
    #pdb.set_trace()
    #X, labels = shuffle_two_array(X, np.array(labels))
    train_data1, train_data2 = transfer_array_format(train)
    print train_data1.shape, train_data2.shape
    train_data1, scaler1 = preprocess_data(train_data1)
    train_data2, scaler2 = preprocess_data(train_data2)
    
    test_data1, test_data2 = transfer_array_format(test)
    print test_data1.shape, test_data2.shape
    test_data1, scaler1 = preprocess_data(test_data1, scaler = scaler1)
    test_data2, scaler2 = preprocess_data(test_data2, scaler = scaler2)
    trainlabel_new, encoder = preprocess_labels(trainlabel)
    num_class = 3
    blend_train = np.zeros((train.shape[0], num_class)) # Number of training data x Number of classifiers
    blend_test = np.zeros((test.shape[0], num_class)) # Number of testing data x Number of classifiers 
    skf = list(StratifiedKFold(trainlabel, num_class))  
    class_index = 0
    prefilter_train, prefilter_test, prefilter_train_bef, prefilter_test_bef = autoencoder_two_subnetwork_fine_tuning(train_data1, train_data2, 
                                                                            trainlabel_new, test_data1, test_data2)

    print 'deep autoencoder'
    clf = RandomForestClassifier(n_estimators=50)
    clf.fit(prefilter_train, trainlabel)
    ae_y_pred_prob = clf.predict_proba(prefilter_test)[:,1]
    
    proba = transfer_label_from_prob(ae_y_pred_prob)
    

    #all_performance.append([acc, precision, sensitivity, specificity, MCC])
    get_blend_data(class_index, RandomForestClassifier(n_estimators=50), skf, prefilter_test, prefilter_train, np.array(trainlabel), blend_train, blend_test)
    
    print 'deep autoencoder without fine tunning'
    class_index = class_index + 1
    clf = RandomForestClassifier(n_estimators=50)
    clf.fit(prefilter_train_bef, trainlabel)
    ae_y_pred_prob = clf.predict_proba(prefilter_test_bef)[:,1]
    
    proba = transfer_label_from_prob(ae_y_pred_prob)


    #all_performance_bef.append([acc, precision, sensitivity, specificity, MCC])
    get_blend_data(class_index, RandomForestClassifier(n_estimators=50), skf, prefilter_test_bef, prefilter_train_bef, np.array(trainlabel), blend_train, blend_test)
    
    print 'random forest'
    class_index = class_index + 1
    prefilter_train = np.concatenate((train_data1, train_data2), axis = 1)
    prefilter_test = np.concatenate((test_data1, test_data2), axis = 1)
    
    clf = RandomForestClassifier(n_estimators=50)
    clf.fit(prefilter_train, trainlabel)
    ae_y_pred_prob = clf.predict_proba(prefilter_test)[:,1]
    proba = transfer_label_from_prob(ae_y_pred_prob)  

    #all_performance_rf.append([acc, precision, sensitivity, specificity, MCC])
    
    #if not (datatype == 'MSL' or datatype == 'PRC2' or datatype == 'MRP'):
    get_blend_data(class_index, RandomForestClassifier(n_estimators=50), skf, prefilter_test, prefilter_train, np.array(trainlabel), blend_train, blend_test)
    
        
    blend_train, scaler1 = preprocess_data(blend_train, stand = False)
    blend_test, scaler1 = preprocess_data(blend_test, scaler = scaler1, stand = False)
    
    bclf = LogisticRegression()
    bclf.fit(blend_train, trainlabel)
    proba_pred = bclf.predict_proba(blend_test)[:, 1]
    #proba = bclf.predict(blend_test)
    print 'stacked ensembling'
    for pro, pair in zip(proba_pred, pairs):
        print pair, pro


def calculate_acc(label_len, pred, label):
    num =0
    for val1, val2 in zip(pred, label):
        if val1 == val2:
            num = num + 1
    return float(num)/label_len


parser = argparse.ArgumentParser(description="""IPMiner: Hidden ncRNA-protein interaction sequential pattern mining with stacked autoencoder for accurate computational prediction""")

parser.add_argument('-dataset',
                    type=str, help='which dataset you want to do 5-fold cross-validation')

parser.add_argument('-r',
                    type=str, help='RNA fasta file to store RNAs')

parser.add_argument('-p',
                    type=str, help='protein fasta file to store proteins')

args = parser.parse_args()
dataset = args.dataset
if dataset is not None:
    IPMiner(dataset)
else:
    RNA_file = args.r
    protein_file = args.p
    if RNA_file is None or protein_file is None:
        print 'you must input RNA and protein fasta file'
    predict_new_samples(RNA_file, protein_file)
