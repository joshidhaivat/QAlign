import numpy as np
import sys
sys.path.append('/home/djjoshi')
import scipy.io as sio
import edit_distance as ed
import time
from multiprocessing import Pool
import os
import pdb
import Levenshtein as L
import all_functions as all_func
from functools import partial
import multiprocessing as mp

#align = []
reads = []
ref = []
percent = []
avg_time = []
reverse = 0
MAX_PROCESS = 78#int(0.95*mp.cpu_count())
MAX_TPP = 100

def func(i,align,align_ed,cigar):
	global percent,avg_time,ref,reads,MAX_PROCESS
	percent.append(1)
	progress = np.sum(np.array(percent))
	start = time.time()
	alignments = np.zeros([11], dtype = 'f')
	read_index = int(align[i,0])
	read_length = align[i,1]
	read_start = int(align[i,2])
	if read_start < 1:
		read_start = 1
	read_end = int(align[i,3])
	if read_end > read_length-1:
		read_end = read_length
	strand = align[i,4]
	ref_index = int(align[i,5])
	ref_length = int(align[i,6])
	match_start = int(align[i,7])
	if match_start < 1:
		match_start = 1
	match_end = int(align[i,8])
	if match_end > ref_length-1:
		match_end = ref_length
	match_bases = align[i,9]
	total_bases = align[i,10]
	map_quality = align[i,11]
	read = reads[read_index-1]
	temp_score = align_ed[i]
	if(map_quality>10):
		#ed_start = time.time()
		#dist = 0
		if(strand==0 or reverse==1):
				dist = L.distance(ref[ref_index-1][match_start-1:match_end-1].upper(),ed.revcom(read.upper()))
				if(reverse==1):
					dist1 = L.distance(ref[ref_index-1][match_start-1:match_end-1].upper(),ed.revcom(read.upper())[read_start-1:read_end-1])
				elif(strand==0):
					dist1 = L.distance(ref[ref_index-1][match_start-1:match_end-1].upper(),ed.revcom(read[read_start-1:read_end-1].upper()))
		else:
				dist = L.distance(ref[ref_index-1][match_start-1:match_end-1].upper(),read.upper())
				dist1 = L.distance(ref[ref_index-1][match_start-1:match_end-1].upper(),read[read_start-1:read_end-1].upper())
		#ed_time = time.time()-ed_start
		score = dist/read_length
		score1 = dist1/max(match_end-match_start,read_end-read_start)
		if(total_bases/read_length >= 0.9):
			alignments[0] = 1
		alignments[1] = read_index
		alignments[2] = read_length
		alignments[3] = match_start
		alignments[4] = match_end
		alignments[5] = strand
		alignments[6] = score
		alignments[7] = ref_index
		alignments[8] = temp_score
		alignments[9] = len(ref[ref_index-1])
		alignments[10] = score1
		avg_time.append(time.time()-start)
		if progress%20==0:
			pass#print('Progress = {:.2f}%; Score = {:.4f}; Avg. Time per iteration = {:.4f}'.format((progress*MAX_PROCESS*100)/(align.shape[0]),score,np.sum(np.array(avg_time))/len(avg_time)))
	return alignments

if (__name__ == "__main__"):
	t1 = time.time()
	import pickle
	#with open('reads.txt','rb') as f:
	#	reads = pickle.load(f)
	reads,_,_ = all_func.get_reads_from_fasta('reads.fasta')
	ref,_,_ = all_func.get_genome_from_fasta('GRCh38_ref.fasta')
	reads_num = len(reads)
	print('reads = '+str(reads_num))
	align_results = np.zeros([reads_num,4],dtype='i')
	align_analysis = np.zeros([reads_num,12],dtype='i')
	score = np.ones([reads_num,8],dtype='f')
	filename = str(sys.argv[1])#str(input("Please provide the file name to analyse: "))
	output_file = sys.argv[2]#str(input("Please provide the output file name: "))
	col = int(sys.argv[3])#int(input("Please enter Quant_level[1 2 3 4]: "))
	col -= 1
	#matfile = sio.loadmat(name+'.mat')
	#align = matfile[name]
	align,align_ed,cigar = all_func.extract_from_paf(filename,1)
	if(col==0):
		align_results = np.zeros([reads_num,4],dtype='i')
		align_analysis = np.zeros([reads_num,12],dtype='i')
		score = np.ones([reads_num,8],dtype='f')
		minimap_score = np.zeros([reads_num,4],dtype='f')
		ref_align = np.zeros([reads_num,8],dtype='i')
		loop=1
	else:
		rc_filename = 'rc_'+filename
		#matfile = sio.loadmat(rc_name+'.mat')
		#rc_align = matfile[rc_name]
		rc_align,rc_align_ed,rc_cigar = all_func.extract_from_paf(rc_filename,1)
		loop=2
		#align = np.concatenate((align,rc_align),axis=0)
		acgt_name = sys.argv[4]#str(input("ACGT analysis file name: "))
		acgt_file = sio.loadmat('align_results_'+acgt_name+'_v2.mat')
		align_results = acgt_file['align_results']
		align_analysis = acgt_file['align_analysis']
		score = acgt_file['score']
		minimap_score = acgt_file['minimap_score']
		ref_align = acgt_file['ref_align']
	#with open('ref.txt','rb') as f:
	#	ref = pickle.load(f)
	print('number of Chromosomes = '+str(len(ref)))
	print('Loading complete!')
	print('Total iterations = '+str(align.shape[0]))
	print('Total reads = '+str(reads_num))
	
	alignments = []
	for j in range(0,loop):
		if(j==1):
			reverse=1
			align = rc_align
			cigar = rc_cigar
			align_ed = rc_align_ed
		#total_loops = int(align.shape[0]/(MAX_PROCESS*MAX_TPP))+1
		#for i in range(0,total_loops):
		#print('Starting loop '+str(i)+' of '+str(total_loops))
		#for i in range(align.shape[0]):
		#	alignments.append(func(i,align=align,align_ed = align_ed, cigar = cigar, ref = ref, reads = reads))
		func1 = partial(func,align=align,align_ed = align_ed, cigar = cigar)
		print('Starting the loop with multiprocessing')
		p = Pool(processes = MAX_PROCESS)
		#i_start = i*MAX_PROCESS*MAX_TPP
		#if(i != total_loops-1):
		#	i_end = (i+1)*MAX_PROCESS*MAX_TPP
		#else:
		#	i_end = align.shape[0]
		alignments += p.map(func1,range(align.shape[0]))#range(i_start,i_end))
		p.close()
		p.join()
	#alignments=[]
	#for i in range(0,align.shape[0]):
	#	alignments+=[func(i)]
	alignments = np.array(alignments)
	#print(alignments.shape)
	#sio.savemat(output_file+'.mat',{'read_alignments':alignments})

	for i in range(alignments.shape[0]):
		idx = int(alignments[i,1])-1
		temp_idx = alignments[i,0]
		col2 = 3*col
		if(alignments[i,1]!=0):
			write=0
			read_len = alignments[i,2]
			scr = alignments[i,6]
			temp_scr = score[idx,2*col]
			if(temp_scr==1):
				write=1
			else:
				if(scr<=temp_scr and scr!=0):
					write=1
			if(write==1):
				align_results[idx,col] = alignments[i,0]
				score[idx,2*col] = scr
				score[idx,2*col+1] = alignments[i,10]
				minimap_score[idx,col] = alignments[i,8]
				ref_align[idx,col*2] = alignments[i,7]
				ref_align[idx,col*2+1] = alignments[i,9]
				align_analysis[idx,col2+1] = alignments[i,3]
				align_analysis[idx,col2+2] = alignments[i,4]
				if(col==0):
					align_analysis[idx,col2] = alignments[i,7]
				else:
					if(align_analysis[idx,0]==0):
						align_analysis[idx,col2] = 2
					else:
						if(((abs(align_analysis[idx,1]-align_analysis[idx,col2+1]) <= 0.5* read_len) or (abs(align_analysis[idx,2]-align_analysis[idx,col2+2]) <= 0.5*read_len)) and align_analysis[idx,0]==alignments[i,7]):
							align_analysis[idx,col2] = 1
						else:
							align_analysis[idx,col2] = -1
	print('Total time = '+str(time.time()-t1))
	sio.savemat('align_results_'+output_file+'_v2.mat',{'align_results':align_results,'align_analysis':align_analysis,'score':score,'minimap_score':minimap_score,'ref_align':ref_align})
	#Provide a bit of analysis
	s_h = np.zeros_like(align_results)
	for i in range(align_results.shape[1]):
		num_alignments = np.sum(align_results[:,i])
		if i==0:
			print('Percentage of well aligned reads in ACGT = {} out of {} ({:.2f}%)'.format(num_alignments,align_results.shape[0],num_alignments*100/align_results.shape[0]))
		else:
			print('Percentage of well aligned reads in Q{} = {} out of {} ({:.2f}%)'.format(i+1,num_alignments,align_results.shape[0],num_alignments*100/align_results.shape[0]))
	for i in range(align_results.shape[1]):
		s_h[score[:,2*i+1]<=0.48,i] = 1
		num_align = np.sum(s_h[:,i])
		if i==0:
			print('Percentage of aligned reads in ACGT = {} out of {} ({:.2f}%)'.format(num_align,s_h.shape[0],num_align*100/s_h.shape[0]))
		else:
			print('Percentage of aligned reads in Q{} = {} out of {} ({:.2f}%)'.format(i+1,num_align,s_h.shape[0],num_align*100/s_h.shape[0]))
	print('All Done!')
