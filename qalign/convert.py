import qalign.all_functions as func
import numpy as np
import os
import subprocess
import pdb

def run_cmd(cmd, shell=True, show_cmd=True):
  if show_cmd:
    print(cmd)
  subprocess.call(cmd, shell=shell)

def get_kmermap_default():
  path = os.path.dirname(os.path.abspath(__file__))
  path = os.path.join(path, 'kmermap_models',
    'r9.4_6mer_nucleotide_template_model.txt')
  kmermap = func.get_kmer_map(path)
  return kmermap

def get_kmermap(kmermap_path=None):
  """
  Input:
    kmermap_path: kmer model file absolute path. If None (default),
      use qaline/kmermap_models/r9.4_6mer_nucleotide_template_model.txt

  Output:
    kmermap: a dic with key as kmer and val as [mean, std]
  """
  if kmermap_path is not None:
    try:
      kmermap = func.get_kmer_map(kmermap_path)
    except:
      kmermap = get_kmermap_default()
  else:
    kmermap = get_kmermap_default()

  return kmermap

def get_qlevels(kmermap, num_levels):
  """
  Input:
    kmermap: a dic with key as kmer and val as [mean, std]
    num_levels: 2 (for Q2), 3 (for Q3), 4 (for ACGT)
  Output:
    qlevels: an array of num_levels+1 elements for quantization
  
  Note:
    use kmermap to construct a current seq to get qlevels,
    in order to avoid calculating qlevels that was originally calculated
    using long reference genome
  """
  current_seq = [
    kmermap[k][0] for k in kmermap.keys()]
  current_seq = np.array(current_seq)
  qlevels = func.get_quantize_level(current_seq, num_levels)

  return qlevels

def convert(input_fasta, output_dir, qlevel, rc=False, kmermap_path=None):
  """
  Input:
    input_fasta: a fasta file (could be read, ref genome etc) in ACGT format
    qlevel: 2 for Q2 and 3 for Q3
    rc: reverse complement
      True: typically for reads
      False: typically for ref genome
    kmermap_path: absolute path for kmer model, or default None
  Output:
    if rc is False:
      a quantized file will be stored at output_dir/res.fasta
    if rc is True:
      two quantized files will be stored at output_dir/res.fasta and rc_res.fasta
  """
  reads, num_reads, name = func.get_reads_from_fasta(input_fasta)

  kmermap = get_kmermap(kmermap_path)
  kmer_k = len(list(kmermap.items())[0][0])

  qlevels = get_qlevels(kmermap, qlevel)

  # current_reads = [None]*num_reads
  # rc_current_reads = [None]*num_reads
  reads_q = []
  rc_reads_q = []

  num_k = 2 if rc is True else 1

  for i in range(0, num_reads):
    for k in range(0, num_k):
      seq = reads[i]
      if(k==1):
        seq = func.revcom(seq)
      current_mean = np.zeros(len(seq)-kmer_k+1)
      for j in range(0,len(seq)-kmer_k+1):
        try:
          # exception when read/ref has unknown kmer, or ref has 'N'/'M'/'R' in kmer
          current_mean[j] = kmermap[seq[j:j+kmer_k]][0]
        except:
          continue
      if(k==0):
        # current_reads[i] = current_mean
        reads_q += [func.get_quantized_seq(current_mean, qlevels)]
      else:
        # rc_current_reads[i] = current_mean
        rc_reads_q += [func.get_quantized_seq(current_mean, qlevels)]
      current_mean = []
    if(i%1000==0):
      print(i)

  run_cmd('mkdir -p %s'%output_dir)
  func.print_reads_to_fasta(reads_q, '%s/res.fasta'%output_dir)
  if rc is True:
    func.print_reads_to_fasta(rc_reads_q, '%s/rc_res.fasta'%output_dir)

  return