from multiprocessing import Pool
import sys
sys.path.append('/home/djjoshi')
import all_functions as all_func
import numpy as np
import scipy.io as sio
import time
import multiprocessing as mp
#import edit_distance as ed
import os

tol = []#int(sys.argv[2])
ed_score_thd = float(sys.argv[4])
mapping_thd = 0
overlap = []
reads_overlap = []
overlap_ed = []
thd = []#[int(sys.argv[3])]
n = []
reads = []
percent = 0
MAX_PROCESS = 80#int(mp.cpu_count())
MAX_TPP = 1000
total_loops = []
reverse = 0
overlap_factor = 0.9

def func(i):
	#print(f'Process {os.getpid()} started')
	global overlap,overlap_ed,n,thd,reads,percent,MAX_PROCESS,tol,reverse,reads_overlap,ed_score_thd,overlap_factor,mapping_thd
	percent += 1
	t1 = time.time()
	ovp_set = np.zeros([4,1],dtype='f')
	first_index = int(overlap[i,0])
	read_length1 = int(overlap[i,1])
	match_start1 = int(overlap[i,2])
	if(match_start1<1):
		match_start1 = 1
	match_end1 = int(overlap[i,3])
	strand = int(overlap[i,4])
	second_index = int(overlap[i,5])
	read_length2 = int(overlap[i,6])
	match_start2 = int(overlap[i,7])
	if(match_start2<1):
		match_start2 = 1
	match_end2 = int(overlap[i,8])
	matched_base = int(overlap[i,9])
	matched_total = int(overlap[i,10])
	ed_score = overlap_ed[i]
	mapping_quality = int(overlap[i,11])
	pair = np.array([first_index,second_index])
	pair.sort(axis=0)
	temp_index = 0
	for j in range(0,pair[0]-1):
		temp_index += n-j
	temp_index += pair[1]-pair[0]+1
	#print('tmp_idx='+str(temp_index))
	g1 = match_end1-match_start1
	g2 = match_end2-match_start2
	l = reads_overlap[temp_index-1,2]
	write = 0
	if((matched_total/min([read_length1,read_length2])) >= 0.85):
		tol1=0
		tol2=0
		diff1 = max(0,g1+tol1-l)
		diff2 = max(0,g2+tol2-l)
		r_l = min(read_length1,read_length2)
		if(g1 >= overlap_factor*r_l and g2 >= overlap_factor*r_l and (ed_score <= ed_score_thd or mapping_quality >= mapping_thd)):
			write=1
	else:
		tol1 = min(match_start1,read_length1-match_end1)
		tol2 = min(match_start2,read_length2-match_end2)
		diff1 = max(0,g1+tol1-l)
		diff2 = max(0,g2+tol2-l)
		if((tol1+tol2 <= (1-overlap_factor)*g1) and (tol1+tol2 <= (1-overlap_factor)*g2) and ed_score <= ed_score_thd):
			write = 1
	if(l!=0):
		score = (g1+g2-diff1-diff2)/(l+l+diff1+diff2)
	else:
		score = 0
	if(write==1):
		ovp_set[0] = 1
		ovp_set[1] = temp_index
		ovp_set[2] = ed_score
		ovp_set[3] = score
		#if(i%1000 == 0):
		#print('i='+str(percent)+'/'+str(MAX_TPP)+'; score='+str(score)+'; t='+str(time.time()-t1))
	return ovp_set

if __name__ == "__main__":
	t1 = time.time()
	foldername = str(sys.argv[1])#str(input("Enter the number : "))
	filename = str(sys.argv[2])
	out_filename = str(sys.argv[3])
	#matfile = sio.loadmat('chosen_reads.mat')
	#reads_list = matfile['chosen_reads']
	#for i in range(reads_list.shape[1]):
	#	reads += [reads_list[0,i][0]]
	reads,_,reads_name = all_func.get_reads_from_fasta(foldername+'reads.fasta')
	overlap,overlap_ed,cigar = all_func.extract_from_paf(foldername+filename,1)
	reads_overlap = sio.loadmat(foldername+'ground_truth.mat')
	reads_overlap = reads_overlap['reads_overlap']
	#reads_overlap = reads_overlap[:,:,1]
	n = int(len(reads))
	print('n='+str(n))
	print(len(thd))
	if 'q' in filename:
		loops = 2
	else:
		loops = 1
	for ll in range(0,loops):
		if(ll==0):
			print('Loading file : {}'.format(filename))
		else:
			overlap,overlap_ed,cigar = all_func.extract_from_paf(foldername+'rc_'+filename,1)
			print('Loading file : {}'.format('rc_'+filename))
			reverse = 1
		#ovp_set = np.zeros([overlap.shape[0],2,len(thd)],dtype='float')
		print('Loading complete')
		ovp_set = []
		print('# of CPUs = '+str(MAX_PROCESS))
		time.sleep(1)
		total_loops = int(overlap.shape[0]/(MAX_PROCESS*MAX_TPP))+1
		for i in range(0,total_loops):
			print('Starting loop '+str(i)+' of '+str(total_loops))
			p = Pool(processes = MAX_PROCESS, maxtasksperchild = 10)
			i_start = i*MAX_PROCESS*MAX_TPP
			if(i != total_loops-1):
				i_end = (i+1)*MAX_PROCESS*MAX_TPP
			else:
				i_end = overlap.shape[0]
			ovp_set += p.map(func,range(i_start,i_end))
			p.close()
			p.join()
		ovp_set = np.array(ovp_set)
		#print(ovp_set)
		#print(ovp_set.shape)
		#sio.savemat('ovlp_'+name+'.mat',{'ovp_set':ovp_set})
	
		if(ll==0):
			overlap_set = np.concatenate((reads_overlap, np.zeros([reads_overlap.shape[0],3],dtype='f')),axis=1)
		for i in range(0,ovp_set.shape[0]):
			if(ovp_set[i,0]!=0):
				write = 0
				index = int(ovp_set[i,1])
				score = ovp_set[i,2]
				temp_score = overlap_set[index-1,-2]
				if(temp_score!=0):
					if(score<=temp_score):
						write = 1
				else:
					write = 1
				if(write):
					overlap_set[index-1,-2] = score
					overlap_set[index-1,-3] = ovp_set[i,0]
					overlap_set[index-1,-1] = ovp_set[i,3]
					#print('k='+str(k)+' i='+str(i))
	ovp = overlap_set
	a1 = ovp[:,3]==1
	a = np.sum(a1)
	print('Ground truth = '+str(a))
	a2 = ovp[:,-3]==1
	a = (a1 & a2)
	a = np.sum(a)
	print('True positives = '+str(a))
	a1 = ovp[:,3]==0
	a2 = ovp[:,-3]==1
	a = (a1 & a2)
	a = np.sum(a)
	print('False positives = '+str(a))
	sio.savemat(foldername+'overlap_analysis_'+out_filename+'.mat',{'overlap_set':overlap_set})
	print('Done! '+foldername+filename+'_'+str(tol)+'  Time taken = '+str(time.time()-t1))
