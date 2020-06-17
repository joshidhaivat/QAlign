import numpy as np
import math
import pdb
import time
from random import choice

def edit_distance(string_x=None, string_y=None):
    matrix_a = np.zeros([2, len(string_y) + 1], dtype = int)
    matrix_a[0, :] = np.arange(0,len(string_y)+1,1)
    col = np.arange(0,len(string_x)+1,1)

    for i in range(1,len(col)):
        matrix_a[1, 0] = col[i]
        for j in range(1,len(matrix_a[1,:])):
            distVer = matrix_a[0, j] + 1
            distHor = matrix_a[1, j - 1] + 1
            if (string_x[i - 1] == string_y[j - 1]):
                distDiag = matrix_a[0, j - 1]
            else:
                distDiag = matrix_a[0, j - 1] + 1
            matrix_a[1, j] = min(([distVer, distHor, distDiag]))
        matrix_a[0, :] = matrix_a[1, :]

    editDistance = matrix_a[-1, -1]
    return editDistance


def edit_distance_gpu(string_x=None, string_y=None):
    import tensorflow as tf

    USE_GPU = False
    def select_device(use_gpu=True):
        from tensorflow.python.client import device_lib
        device = '/device:GPU:0' if use_gpu else '/CPU:0'
        return device

    device = select_device(use_gpu=USE_GPU)
    sess = tf.Session()
    # Create input data

    def create_sparse_vec(word_list):
        num_words = len(word_list)
        indices = [[xi, 0, yi] for xi,x in enumerate(word_list) for yi,y in enumerate(x)]
        chars = list(''.join(word_list))
        return tf.SparseTensorValue(indices, chars, [num_words,1,1])

    string_x_sparse = create_sparse_vec([string_x])
    string_y_sparse = create_sparse_vec([string_y])

    ed = sess.run(tf.edit_distance(string_x_sparse, string_y_sparse, normalize=False))
    return ed[0][0]


def complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)


def revcom(s):
    return complement(s[::-1])

def mash_distance(seq1=None,seq2=None):
        fid = open("mash_function.fa","w")
        fid.write(">read_1\n")
        fid.write("%s\n" % seq1)
        fid.write(">read_2\n")
        fid.write(seq2)
        fid.close()
        import os
        os.system("sh run.sh")
        fid = open("mash_result.dist","r")
        for aline in fid:
            a = aline.split()
            #print(a)
            if(a[0]=='read_2'):
                distance = float(a[2])
            else:
                distance = 0
        fid.close()
        #print('distance = '+str(distance))
        return distance
    
def get_kmer_map(filename):
    start=0
    myDict = {}
    filename = str(filename)
    fid = open(filename,'r')
    for aline in fid:
        a = aline.split()
        if(a[0]=='kmer'):
            start=1
            continue
        if(start):
            kmer = str(a[0])
            value = [float(a[1]), float(a[2])]
            myDict[kmer] = value
    fid.close()
    return myDict

def get_random_map(filename,quant_level):
    start=0
    myDict = {}
    filename = str(filename)
    fid = open(filename,'r')
    dict = {1:'A',2:'C',3:'G',4:'T'}
    for aline in fid:
        a = aline.split()
        if(a[0]=='kmer'):
            start=1
            continue
        if(start):
            kmer = str(a[0])
            import numpy as np
            value = np.random.random_integers(quant_level)
            value = dict[value]
            myDict[kmer] = value
    fid.close()
    return myDict

def get_reads_from_fastq(filename):
    reads = []
    name = []
    num_reads=0
    filename = str(filename)
    fid = open(filename,'r')
    seq_indicator=0
    for aline in fid:
        if(aline[0]=='@'):
            seq_indicator=1
            name.append(aline)
            continue
        if(seq_indicator==1):
            seq=aline
            num_reads+=1
            seq_indicator=0
            reads+=[seq[:-1]]
    print('Total_reads = '+str(num_reads))
    fid.close()
    return reads,num_reads,name

def get_reads_from_fasta(filename):
    reads = []
    read_name = []
    num_reads=0
    filename = str(filename)
    fid = open(filename,'r')
    seq_indicator = 0
    for aline in fid:
        if(aline[0]=='>'):
            seq_indicator = 1
            read_name.append(aline)
            continue
        if(seq_indicator==1):
            seq = aline
            num_reads+=1
            seq_indicator = 0
            reads.append(seq[:-1])
    print('Total_reads = '+str(num_reads))
    fid.close()
    return reads,num_reads,read_name

def get_genome_from_fasta(filename):
    genome = []
    num_genome = 0
    filename = str(filename)
    fid = open(filename,'r')
    seq=''
    name=[]
    for aline in fid:
        if(aline[0]=='>'):
            #pdb.set_trace()
            num_genome+=1
            name.append(aline)
            if not seq:
                continue
            else:
                genome.append(seq)
            seq=''
        else:
            seq+=aline[:-1]
    genome.append(seq)
    fid.close()
    print('Total genomes = '+str(num_genome))
    return genome,num_genome,name

def print_reads_to_fasta(reads,filename,name=None):
    filename = str(filename)
    fid = open(filename,'w')
    for i in range(0,len(reads)):
        seq = reads[i]
        if name == None:
            fid.write('>read_%d\n' % (i+1))
        else:
            fid.write('%s' % name[i])
        fid.write('%s\n\n' % seq)
    fid.close()
    print('Done!')   

def print_ref_to_fasta(ref,filename,name=None):
    filename = str(filename)
    fid = open(filename,'w')
    for i in range(0,len(ref)):
        seq = ref[i]
        if(name==None):
            if(len(ref)>1):
                fid.write(">ref_%d\n" % (i+1))
            else:
                fid.write(">ref\n")
        else:
            fid.write(name[i])
        fid.write('%s\n\n' % seq)
    fid.close()
    print('Done!')

def get_quantize_level(current_seq,level_num):
    qlevels = np.zeros(level_num+1)
    current_seq = np.sort(current_seq)
    a = np.where(current_seq>0)
    #pdb.set_trace()
    current_seq = current_seq[a]
    #pdb.set_trace
    interval = len(current_seq)/level_num
    #pdb.set_trace
    for i in range(0,len(qlevels)):
        if i==0:
            qlevels[i] = current_seq[0]
        elif i==len(qlevels)-1:
            qlevels[i] = current_seq[-1]
        else:
            qlevels[i] = current_seq[math.ceil(i*interval)]
    return qlevels

def get_quantized_sequences(current_seq,qlevels):
    quantized_seq = [None]*len(current_seq)
    quantized_value = 'ACGT'
    quantize_num = len(qlevels)-1
    for i in range(0,len(current_seq)):
        write_status = 0
        for j in range(1,quantize_num+1):
            if((current_seq[i]<=qlevels[j]) and (current_seq[i]>=qlevels[j-1])):
                quantized_seq[i] = quantized_value[j-1]
                write_status = 1
                break
        if(write_status==0):
            quantized_seq[i]='N'#quantized_value[quantize_num-1]
    quantized_seq = ''.join(quantized_seq)
    #pdb.set_trace()
    return quantized_seq

def get_quantized_seq(current_seq,qlevels):
    # Faster computation than get_quantized_sequences by 4 times
    quant_seq = -1*np.ones(len(current_seq),dtype=int)
    quant_value = {-1:'N',0:'A',1:'C',2:'G',3:'T'}
    quant_num = len(qlevels)-1
    for j in range(1,quant_num+1):
        quant_seq[(current_seq>=qlevels[j-1]) & (current_seq<=qlevels[j])] = j-1
    quantized_seq = [quant_value[x] for x in quant_seq]
    quantized_seq = ''.join(quantized_seq)
    return quantized_seq

def get_index(i,j,n):
    pair = np.array([i,j])
    pair.sort(axis=0)
    index=0
    for jj in range(0,pair[0]-1):
        index += n-jj
    index += pair[1]-pair[0]+1
    return index


def extract_from_paf(filename,cigar_string=False):
    start_time = time.time()
    fid = open(filename,'r')
    out = []
    out_ed = []
    cigar = []
    for aline in fid:
        import pdb
        #pdb.set_trace()
        temp = -1*np.ones((16),dtype=int)
        cols = aline.split('\t')
        seq1 = cols[0] # seq1_name
        if seq1[:5] == 'Read_' or seq1[:5] == 'read_':
            try:
                temp[0] = int(seq1[5:])
            except:
                print('Read number is not available\n')
        else:
            seq1 = seq1.split('.')
            try:
                temp[0] = int(seq1[1])
            except:
                print('Read number is not available\n')
        #pdb.set_trace()
        temp[1] = int(cols[1]) # seq1_length
        temp[2] = int(cols[2]) # seq1_start
        temp[3] = int(cols[3]) # seq1_end
        strand = cols[4] # Forward or reverse strand alignment
        if strand == "+":
            temp[4] = 1
        else:
            temp[4] = 0
        seq2 = cols[5] # seq2_name
        if seq2[:8] == 'ref_chr_':
            if seq2[8:] == 'x' or seq2[8:]=='X':
                temp[5] = 23
            elif seq2[8:] == 'y' or seq2[8:] == 'Y':
                temp[5] = 24
            elif seq2[8:] == 'm' or seq2[8:] == 'M':
                temp[5] = 25
            else:
                try:
                    temp[5] = int(seq2[8:])
                except:
                    print('Reference number is not available\n')
        elif seq2[:4] == 'ref_':
            try:
                temp[5] = int(seq2[4:])
            except:
                print('Reference number is not available\n')
        elif seq2 == 'ref':
            temp[5] = 1
        elif seq2[:5] == 'read_' or seq2[:5] == 'Read_':
            try:
                temp[5] = int(seq2[5:])
            except:
                print('Reference number is not available\n')
        else:
            print('Unidentified form of ref name')
        temp[6] = int(cols[6]) # Seq2_length
        temp[7] = int(cols[7]) # Seq2_start
        temp[8] = int(cols[8]) # Seq2_end
        temp[9] = int(cols[9]) # Match bases
        temp[10] = int(cols[10]) # total_match
        temp[11] = int(cols[11]) # mapping quality
        out_ed.append((temp[10] - temp[9])/temp[10])
        for ii in range(12,len(cols)):
            #pdb.set_trace()
            if cigar_string==True and len(cols[ii])>5 and cols[ii][:5]=='cg:Z:':
                cigar_seq = cols[ii][5:-1]
                cigar.append(cigar_seq)
            elif len(cols[ii])>5 and cols[ii][:5] == 'tp:A:':
                map_type = cols[ii]
                if map_type[5] == 'P':
                    temp[12] = 1
                elif map_type[5] == 'S':
                    temp[12] = 0
                else:
                    print('Unknown mapping type: {}\n'.format(cols[ii]))
            elif len(cols[ii])>5 and cols[ii][:5] == 'cm:i:':
                temp[13] = int(cols[ii][5:]) # Num of anchors
            elif len(cols[ii])>5 and cols[ii][:5] == 's1:i:':
                temp[14] = int(cols[ii][5:]) # Chaining score of Primary alignment
            elif len(cols[ii])>5 and cols[ii][:5] == 's2:i:':
                temp[15] = int(cols[ii][5:]) # Chaining score of Secondary alignment
        out.append(temp)
    fid.close()
    output = np.matrix(out)
    output_ed = np.array(out_ed)
    print('The dimension of the extracted PAF file is: {}'.format(output.shape))
    print('Total time taken: {:.2f} seconds'.format(time.time()-start_time))
    if cigar_string==False:
        return output,output_ed
    else:
        return output,output_ed,cigar

def rand_seq(length,alphabet_size):
    # Generates a random seq of the provided length and the provided alphabet size
    DNA=""
    for i in range(length):
        if alphabet_size == 4:
            DNA+=choice("ACGT")
        elif alphabet_size == 3:
            DNA+=choice("ACG")
        elif alphabet_size == 2:
            DNA+=choice("AC")
    return DNA

def extract_from_vcf(filename):
    fid = open(filename,'r')
    out = []
    out_seq = []
    sniffle_version = 11
    for aline in fid:
        precise = -1
        ref_start = -1
        index = -1
        ref_end = -1
        sv_length = 0
        read_support = -1
        ref_chr = 'NONE'
        sv_type = 'NONE'
        if aline[0]!='#':
            cols = aline.split('\t')
            ref_chr = cols[0]
            seq = cols[3]
            temp = cols[7].split(';')
            ref_start = int(cols[1])
            index = int(cols[2])
            for element in temp:
                if element=='PRECISE':
                    precise = 1
                    continue
                elif element == 'IMPRECISE':
                    precise = 0
                    continue
#                 elif len(element)>19 and element[:17] == 'SVMETHOD=Sniffles':
#                     sniffle_version = int(element[-2:])
#                     break
                elif len(element)>7 and element[:7] == 'SVTYPE=':
                    sv_type = element[7:]
                    continue
                elif len(element)>6 and element[:6] == 'SVLEN=':
                    sv_length = int(element[6:])
                    continue
                elif len(element)>4 and element[:4] == 'END=':
                    ref_end = int(element[4:])
                    continue
                elif len(element)>3 and element[:3] == 'RE=':
                    read_support = int(element[3:])
                    continue
            out.append([index,precise,ref_start,ref_end,sv_length,read_support,ref_chr,sv_type])
            out_seq.append(seq)
    fid.close()
    return out,out_seq