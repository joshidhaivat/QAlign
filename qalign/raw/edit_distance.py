import numpy as np
import pdb
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

	USE_GPU = True
	def select_device(use_gpu=True):
		from tensorflow.python.client import device_lib
		device = '/device:GPU:0' if use_gpu else '/CPU:0'
		return device

	device = select_device(use_gpu=USE_GPU)
	sess = tf.Session()
	# Create input data
	#pdb.set_trace()
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
