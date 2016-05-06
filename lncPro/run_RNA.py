import sys
import os
import pdb
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

rna_all = sys.argv[1]
seq_dict = read_fasta_file(rna_all)
fw_all = open('RNA_fea_all.txt', 'w')
for key, val in seq_dict.iteritems():
    fw = open('tmp.fa', 'w')
    fw.write('>' + key + '\n')
    fw.write(val + '\n')
    fw.close()
    #pdb.set_trace()
    cli_str = './RNAScore_rna tmp.fa'
    fex = os.popen(cli_str, 'r')
    fex.close()
    if os.stat('RNA_fea.txt').st_size==0:
        print key
        feas = '\t'.join(map(str, [0] * 30))
        fw_all.write(key + '\t' + feas + '\n')
    else:
        fp = open('RNA_fea.txt', 'r')
        for line in fp:
            fw_all.write(key + '\t' + line)
        fp.close()
fw_all.close()
