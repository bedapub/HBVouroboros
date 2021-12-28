import os,sys
import re

aa_ref = 'ARNDCQEGHILKMFPSTWYV*-X'
nt_ref = 'ACGT-N'

def mapped(flag):
	return False if flag & 0x4 else True

def primary(flag):
	return False if flag & 0x100 else True

def parse_cigar_string(cigar):
	len_list = [ int(l) for l in re.split('M|N|I|D|S|H|X', cigar) if l != '']
	code_list = [ c for c in re.split('1|2|3|4|5|6|7|8|9|0', cigar) if c != '']
	return len_list, code_list

def cigar_seq_to_aligned_seq(cigar, seq):
	len_list, code_list = parse_cigar_string(cigar)
	aligned_seq = ''
	start_loc_in_read_seq = 1
	for tmp in range(len(code_list)):
		tmp_len = len_list[tmp]
		tmp_code = code_list[tmp]
		if tmp_code == 'S':
			start_loc_in_read_seq += tmp_len
		elif tmp_code == 'H':
			print(f"tmp_code == 'H'! Current code cannot handle hard clipping case!")
			exit(1)
		elif tmp_code in 'MX':
			aligned_seq += seq[(start_loc_in_read_seq-1): (start_loc_in_read_seq-1+tmp_len)]
			start_loc_in_read_seq += tmp_len
		elif tmp_code in 'DN':
			aligned_seq += ('-' * tmp_len)
		elif tmp_code == 'I':
			start_loc_in_read_seq += tmp_len
		else:
			print(f'Invalid code: {tmp_code}!')
			exit(1)
	return aligned_seq	

def get_end_pos_from_start_pos_cigar(start, cigar):
	len_list, code_list = parse_cigar_string(cigar)
	start -= 1
	for tmp in range(len(len_list)):
		if code_list[tmp] in 'MNDX':
			start += len_list[tmp]
	return start

def build_nt2indexMap(uniq_nt_seq):
	nt2indexMap = {}
	tmp_nt_index = 0
	for tmp in uniq_nt_seq:
		nt2indexMap[tmp] = tmp_nt_index
		tmp_nt_index += 1
	return nt2indexMap

def build_aa2indexMap(uniq_aa_seq):
	aa2indexMap = {}
	tmp_aa_index = 0
	for tmp in uniq_aa_seq:
		aa2indexMap[tmp] = tmp_aa_index
		tmp_aa_index += 1
	return aa2indexMap

def ntSeq2aaSeq(ntSeq):
	nt_num = len(ntSeq)
	if nt_num%3 != 0:
		print(f"Error,nt_num%3 != 0")
		exit(1)
	aa_seq = ''
	aa_num = nt_num//3
	for tmp in range(aa_num):
		tmp_codon = (ntSeq[tmp*3:tmp*3+3]).upper()
		tmp_aa = nt2aa(tmp_codon)
		aa_seq += tmp_aa
	return aa_seq

def nt2aa(codon):
	if codon[0] == 'A':
		if codon[1] == 'A':
			if codon[2] in 'AG':
				return 'K'
			elif codon[2] in 'CT':
				return 'N'
			else:
				return 'X'
		elif codon[1] == 'C':
			return 'T'
		elif codon[1] == 'G':
			if codon[2] in 'CT':
				return 'S'
			elif codon[2] in 'AG':
				return 'R'
			else:
				return 'X'
		elif codon[1] == 'T':
			if codon[2] in 'ACT':
				return 'I'
			elif codon[2] == 'G':
				return 'M'
			else:
				return 'X'
		else:
			return 'X'
	elif codon[0] == 'C':
		if codon[1] == 'A':
			if codon[2] in 'CT':
				return 'H'
			elif codon[2] in 'AG':
				return 'Q'
			else:
				return 'X'
		elif codon[1] == 'C':
			return 'P'
		elif codon[1] == 'G':
			return 'R'
		elif codon[1] == 'T':
			return 'L'
		else:
			return 'X'
	elif codon[0] == 'G':
		if codon[1] == 'A':
			if codon[2] in 'AG':
				return 'E'
			elif codon[2] in 'CT':
				return 'D'
			else:
				return 'X'
		elif codon[1] == 'C':
			return 'A'
		elif codon[1] == 'G':
			return 'G'
		elif codon[1] == 'T':
			return 'V'
		else:
			return 'X'
	elif codon[0] == 'T':
		if codon[1] == 'A':
			if codon[2] in 'AG':
				return '*'
			elif codon[2] in 'CT':
				return 'Y'
			else:
				return 'X'
		elif codon[1] == 'C':
			return 'S'
		elif codon[1] == 'G':
			if codon[2] == 'A':
				return '*'
			elif codon[2] == 'G':
				return 'W'
			elif codon[2] in 'CT':
				return 'C'
			else:
				return 'X'
		elif codon[1] == 'T':
			if codon[2] in 'AG':
				return 'L'
			elif codon[2] in 'CT':
				return 'F'
			else:
				return 'X'
		else:
			return 'X'
	else:
		return 'X'

def build_seqPos2alignedSeqPosMap(aligned_seq):
	seqPos2alignedSeqPosMap = {}
	pos = 0
	for tmp in range(len(aligned_seq)):
		if aligned_seq[tmp] != '-':
			pos += 1
			seqPos2alignedSeqPosMap[pos] = tmp + 1
	return seqPos2alignedSeqPosMap

def removeDel(aligned_seq):
	seq = ''
	for tmp in aligned_seq:
		if tmp != '-':
			seq += tmp
	return seq

def parse_multi_seq_fa(multi_seq_fa):
	id_list = []
	seq_list = []
	tmp_seq = ''
	fa_lines = open(multi_seq_fa, 'r').read().splitlines()
	id_list.append(fa_lines[0])
	for tmp in fa_lines[1:]:
		if tmp[0] == '>':
			id_list.append(tmp) # next sequence
			seq_list.append(tmp_seq) # last sequence
			tmp_seq = ''
		else:
			tmp_seq += tmp
	seq_list.append(tmp_seq)
	return id_list, seq_list

def mafft_cmd(in_file, out_file):
	return mafft_binary_path + ' ' + in_file + ' > ' + out_file

def alignMultiSeq2ref(ref_id, ref_seq, to_process_id_list, to_process_seq_list, tmp_file_in, tmp_file_out):
	to_process_seq_num = len(to_process_id_list)
	res_id_list_list = []
	res_seq_list_list = []
	for tmp in range(to_process_seq_num):
		#print(f'*********************** tmp:{tmp}')
		tmp_id = to_process_id_list[tmp]
		#print(f'*********************** tmp id: {tmp_id}')
		tmp_seq = to_process_seq_list[tmp]
		#print(f'tmp_id:{tmp_id}')
		#print(f'tmp_seq:{tmp_seq}')
		ofs = open(tmp_file_in, 'w')
		ofs.write(f'{ref_id}\n{ref_seq}\n{tmp_id}\n{tmp_seq}\n')
		ofs.close()
		#exit(1)
		os.system(mafft_cmd(tmp_file_in, tmp_file_out))
		res_id_list, res_seq_list = parse_multi_seq_fa(tmp_file_out)
		res_id_list_list.append(res_id_list)
		res_seq_list_list.append(res_seq_list)
	return res_id_list_list, res_seq_list_list
