import os,sys
from seq_utils import *
import numpy as np

mutation_cutoff_perc_list = [1,2,5,15,50]

def extractOverlapAASeq(gene_start, gene_end, read_start, read_end, read_seq): # return overlap_or_not, aa_seq, start_aa_in_gene
	if gene_start >= read_end or read_start >= gene_end:
		return False, '', -1
	nt_start = max(gene_start, read_start)
	nt_end = min(gene_end, read_end)
	shift_start = (nt_start - gene_start)%3
	if shift_start != 0:
		shift_start = 3 - shift_start
	shift_end = (gene_end - nt_end)%3
	if shift_end != 0:
		shift_end = 3 - shift_end
	nt_start = nt_start + shift_start
	nt_end = nt_end - shift_end
	if nt_end - nt_start + 1 < 3:
		return False, '', -1
	nt_seq = read_seq[shift_start:(shift_start + nt_end - nt_start + 1)]
	aa_seq = ntSeq2aaSeq(nt_seq)
	aa_num = len(aa_seq)
	start_aa_in_gene = (nt_start - gene_start)//3 + 1
	return True, aa_seq, start_aa_in_gene

def pass_mut_cutoff(uniq_nt_or_aa_seq, tmp_freq_list, tmp_total, tmp_ref, tmp_cutoff):
	tmp_mut_str = ""
	if tmp_total == 0:
		return tmp_mut_str
	for tmp_freq_index in range(len(tmp_freq_list)):
		if uniq_nt_or_aa_seq[tmp_freq_index].upper() != tmp_ref.upper() and tmp_freq_list[tmp_freq_index] * 100/tmp_total >= tmp_cutoff:
			tmp_mut_str += f"{uniq_nt_or_aa_seq[tmp_freq_index].upper()};"
	return tmp_mut_str
	
uniq_nt_seq = 'ACGTN-'
nt2indexDic = build_nt2indexMap(uniq_nt_seq)
uniq_nt_num = len(uniq_nt_seq)
uniq_aa_seq = 'ARNDCQEGHILKMFPSTWYV*-X'
aa2indexDic = build_aa2indexMap(uniq_aa_seq)
uniq_aa_num = len(uniq_aa_seq)

# read reference genome
ref_file = 'results/infref/infref_strain_dup.fasta'
ref_id_list, ref_seq_list = parse_multi_seq_fa(ref_file)
ref_id = ref_id_list[0]
ref_seq = ref_seq_list[0]
ref_seq_len = len(ref_seq)
nt_pos_freq_list_list = []
for _ in range(ref_seq_len):
	tmp_freq_list = [0] * uniq_nt_num
	nt_pos_freq_list_list.append(tmp_freq_list)
nt_pos_total_list = [0] * ref_seq_len

gene_start_end_list = [] # [['S', 155, 835]]
gene_start_end_list.append(['S', 155, 835])
gene_start_end_list.append(['X', 1374, 1838])
gene_start_end_list.append(['preC', 1814, 1900])
gene_start_end_list.append(['C', 1901, 2452])
gene_num = len(gene_start_end_list)
gene_aa_ref_list = []
for tmp_gene_index in range(len(gene_start_end_list)):
	tmp_gene_nt_start = gene_start_end_list[tmp_gene_index][1]
	tmp_gene_nt_end = gene_start_end_list[tmp_gene_index][2]
	tmp_gene_nt_ref = ref_seq[(tmp_gene_nt_start-1):tmp_gene_nt_end]
	tmp_gene_aa_ref = ntSeq2aaSeq(tmp_gene_nt_ref)
	gene_aa_ref_list.append(tmp_gene_aa_ref)

pos_aa_freq_list_list_list = []
pos_freq_list_list = []
for tmp in gene_start_end_list:
	tmp_start = tmp[1]
	tmp_end = tmp[2]
	tmp_aa_len = (tmp_end - tmp_start + 1)//3
	tmp_pos_aa_freq_list_list = []
	for tmp_aa_index in range(tmp_aa_len):
		tmp_aa_list = [0] * uniq_aa_num
		tmp_pos_aa_freq_list_list.append(tmp_aa_list)
	tmp_pos_freq_list = [0] * tmp_aa_len
	pos_aa_freq_list_list_list.append(tmp_pos_aa_freq_list_list)
	pos_freq_list_list.append(tmp_pos_freq_list)

sam_file = sys.argv[1]
out = open(f'{sam_file[:-4]}.aaFreq.csv', 'w')
out_nt = open(f'{sam_file[:-4]}.ntFreq.csv', 'w')
sam_lines = open(sam_file, 'r').read().splitlines()
processed_read_num = 0
for tmp in sam_lines:
	if tmp[0] == '@':
		continue
	processed_read_num += 1
	if processed_read_num%10000 == 0:
		print(f"Processed {processed_read_num} reads!")
	#print(f'sam:\n{tmp}')
	tmp_vec = tmp.split('\t')
	tmp_start = int(tmp_vec[3])
	tmp_cigar = tmp_vec[5]
	if 'I' in tmp_cigar or 'D' in tmp_cigar:
		continue
	tmp_seq = tmp_vec[9]
	tmp_aligned_seq = cigar_seq_to_aligned_seq(tmp_cigar, tmp_seq)
	tmp_end = get_end_pos_from_start_pos_cigar(tmp_start, tmp_cigar)
	for tmp_nt_pos in range(tmp_start, tmp_end + 1):
		tmp_nt = tmp_aligned_seq[tmp_nt_pos-tmp_start]
		tmp_nt_index = nt2indexDic[tmp_nt]
		nt_pos_freq_list_list[tmp_nt_pos-1][tmp_nt_index] += 1
		nt_pos_total_list[tmp_nt_pos-1] += 1
	#print(f'start:\n{tmp_start}\ncigar:{tmp_cigar}\nend:{tmp_end}\n')
	for tmp_gene_index in range(gene_num):
		tmp_gene_start_end_tuple = gene_start_end_list[tmp_gene_index]
		tmp_gene = tmp_gene_start_end_tuple[0]
		tmp_gene_start = tmp_gene_start_end_tuple[1]
		tmp_gene_end = tmp_gene_start_end_tuple[2]
		#print(f'Gene:{tmp_gene}\nGene_start:{tmp_gene_start}\nGene_end:{tmp_gene_end}\n')
		overlap_or_not, tmp_aligned_seq_aa, tmp_start_aa_in_gene = extractOverlapAASeq(tmp_gene_start, tmp_gene_end, tmp_start, tmp_end, tmp_aligned_seq)
		if not overlap_or_not:
			continue
		for tmp_aa_pos_index in range(len(tmp_aligned_seq_aa)):
			tmp_aa = tmp_aligned_seq_aa[tmp_aa_pos_index]
			tmp_aa_pos = tmp_start_aa_in_gene + tmp_aa_pos_index
			tmp_aa_index = aa2indexDic[tmp_aa]
			pos_aa_freq_list_list_list[tmp_gene_index][tmp_aa_pos-1][tmp_aa_index] += 1
			pos_freq_list_list[tmp_gene_index][tmp_aa_pos-1] += 1

# print header for ntFreq
out_nt.write(f'NT_pos,Ref,Cons,Total')
for tmp_nt in uniq_nt_seq:
	out_nt.write(f',{tmp_nt}')
for tmp_cutoff in mutation_cutoff_perc_list:
	out_nt.write(f',{tmp_cutoff}%Mut')
out_nt.write('\n')
for tmp_nt_pos_index in range(ref_seq_len):
	tmp_nt_pos = tmp_nt_pos_index + 1
	tmp_nt_ref = ref_seq[tmp_nt_pos_index].upper()
	tmp_cons = uniq_nt_seq[np.array(nt_pos_freq_list_list[tmp_nt_pos_index]).argmax()].upper()
	#tmp_mutateOrNot_str = 'Y' if tmp_nt_ref != tmp_cons else 'N'
	tmp_total = nt_pos_total_list[tmp_nt_pos_index]
	out_nt.write(f'{tmp_nt_pos},{tmp_nt_ref},{tmp_cons},{tmp_total}')
	for tmp_nt_index in range(uniq_nt_num):
		tmp_nt_freq = nt_pos_freq_list_list[tmp_nt_pos_index][tmp_nt_index]
		tmp_nt_perc = round(100*tmp_nt_freq/tmp_total, 2) if tmp_total > 0 else 0
		out_nt.write(f',{tmp_nt_freq}')
	for tmp_cutoff in mutation_cutoff_perc_list:
		tmp_pass_mut_cutoff_str = pass_mut_cutoff(uniq_nt_seq, nt_pos_freq_list_list[tmp_nt_pos_index], tmp_total, tmp_nt_ref, tmp_cutoff)
		out_nt.write(f',{tmp_pass_mut_cutoff_str}')
	out_nt.write(f'\n')
out_nt.close()

# print header for aaFreq
out.write(f'Gene,AA_pos,NT_codon_pos_start,Ref,Cons,Total')
for tmp_aa in uniq_aa_seq:
	out.write(f',{tmp_aa}')
for tmp_cutoff in mutation_cutoff_perc_list:
	out.write(f',{tmp_cutoff}%Mut')
out.write('\n')

for tmp_gene_index in range(gene_num):
	tmp_gene_start_end_tuple = gene_start_end_list[tmp_gene_index]
	tmp_gene = tmp_gene_start_end_tuple[0]
	tmp_gene_start = tmp_gene_start_end_tuple[1]
	tmp_gene_end = tmp_gene_start_end_tuple[2]
	tmp_gene_aa_len = (tmp_gene_end - tmp_gene_start + 1)//3	
	for tmp_aa_pos_index in range(tmp_gene_aa_len):
		tmp_gene_nt_codon_pos_start = tmp_gene_start + tmp_aa_pos_index * 3
		tmp_gene_aa_ref = gene_aa_ref_list[tmp_gene_index][tmp_aa_pos_index].upper()
		tmp_gene_aa_cons = uniq_aa_seq[np.array(pos_aa_freq_list_list_list[tmp_gene_index][tmp_aa_pos_index]).argmax()].upper()
		#tmp_gene_aa_mutateOrNot_str = 'Y' if tmp_gene_aa_ref != tmp_gene_aa_cons else 'N'
		tmp_gene_aa_total_freq = pos_freq_list_list[tmp_gene_index][tmp_aa_pos_index]
		out.write(f'{tmp_gene},{tmp_aa_pos_index+1},{tmp_gene_nt_codon_pos_start},{tmp_gene_aa_ref},{tmp_gene_aa_cons},{tmp_gene_aa_total_freq}')
		for tmp_aa_index in range(uniq_aa_num):
			tmp_aa = uniq_aa_seq[tmp_aa_index]
			tmp_gene_aa_freq = pos_aa_freq_list_list_list[tmp_gene_index][tmp_aa_pos_index][tmp_aa_index]
			tmp_gene_aa_perc = round(100*tmp_gene_aa_freq/tmp_gene_aa_total_freq, 2) if tmp_gene_aa_total_freq > 0 else 0
			out.write(f',{tmp_gene_aa_perc}')
		for tmp_cutoff in mutation_cutoff_perc_list:
			tmp_pass_mut_cutoff_str = pass_mut_cutoff(uniq_aa_seq, pos_aa_freq_list_list_list[tmp_gene_index][tmp_aa_pos_index], tmp_gene_aa_total_freq, tmp_gene_aa_ref, tmp_cutoff)
			out.write(f',{tmp_pass_mut_cutoff_str}')
		out.write(f'\n')
out.close()
