import h5py
import config
import collections
import time
import re
import os
import numpy as np
import pysam
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import glob
import utilities
import subprocess
import shutil
from Bio.Seq import Seq
from cigar import Cigar


# Obsolete don't use! Just for reference.

FNULL = open(os.devnull, 'w')


def get_ref_offset(ref):
    with open('temp_files/temp_index.fa', 'w') as fw:
        fw.write('>temp_index\n' + standard_ref_seq[0:200])
    subprocess.run('bwa mem -M -x ont2d -t 7 ' + ref + 
                   ' temp_files/temp_index.fa > temp_files/'
                   'temp.sam',
                   shell=True, stdout=FNULL,
                   stderr=subprocess.STDOUT)
    with open('temp_files/temp.sam') as f:
        return int(f.readlines()[2].split()[3]) - 1


def calc_qual(qual_str):
    qual_vals = []
    for s in qual_str:
        qual_vals.append(ord(s))
    qual_actual_prob = [10**((-i+33)/10) for i in qual_vals]
    smalls = [i for i in qual_vals if i < 37]
    return np.median(qual_actual_prob)


class Read(object):
    def _is_this_healthy_rDNA(self):
        """Check if the read is completely inside rDNA or at the end.
        """
        if self.length < 3000 or self.direction == 0:
            return 0
        mapping_state = []
        for item in self.sam_summary:
            if item[1] != '0':
                mapping_state.append(1)
            else:
                mapping_state.append(0)
        # sufficient mapping fraction of split mapped reads
        threshold = 0.8
        if sum(mapping_state)/len(mapping_state) > threshold:
            return 1
        else:
            for i in range(1, len(mapping_state) - 50):
                if sum(mapping_state[i:])/len(mapping_state[i:]) > threshold or \
                   sum(mapping_state[:-i])/len(mapping_state[:-i]) > threshold:
                    healthy = 2
        return 0
 
    def _get_direction(self):
        direction_array = self.sam_summary[:, 1]
        minus = 0.1
        plus = 0.1
        for i in direction_array:
            if i == '+':
                plus += 1
            if i == '-':
                minus += 1
        if plus/minus > 10:
            return '+'
        elif minus/plus > 10:
            return '-'
        else:
            return 0

    def _find_basecall_with_mod(self, h5file):
        for i in range(3):
            try:
                bc_template = h5file['Analyses/Basecall_1D_00' + str(i) + '/BaseCalled_template']
            except:
                continue
            if 'ModBaseProbs' in bc_template.keys():
                return str(i)

    def __init__(self, fast5file, ref, mod_loc=0):
        #Use fast5s that are basecalled with modification!
        self.fast5path = fast5file
        with h5py.File(fast5file) as f:
            self.ref_path = ref
            read_name = list(f['Raw/Reads'].keys())[0]
            self.read_id = f['Raw/Reads/' + read_name].attrs['read_id'].decode()
            self.raw_duration = f['Raw/Reads/' + read_name].attrs['duration']
            #self.raw = f['Raw/Reads/'+ read_name + '/Signal'][:]
            #mod_loc = self._find_basecall_with_mod(f)
            basecalled_template = f['Analyses/Basecall_1D_00{}/BaseCalled_template'.format(mod_loc)]
            fastq = basecalled_template['Fastq'][()].decode().strip().split('\n')
            self.seq = fastq[1]
            self.quality = fastq[3]
            self.ave_quality = np.mean([ord(i) - 33 for i in self.quality])
            self.guppy_mbt = basecalled_template['ModBaseProbs'][:]
            self.fastq_str = basecalled_template['Fastq'][()].decode().strip().split('\n')
            #fastq length! not raw signal length
            self.split_length = 300
            self.length = len(self.seq)
            self.sam_summary = np.array(utilities.split_mapping_and_sam_analysis(self.split_length, 'temp', 
                    self.seq, self.quality, ref))
            self.direction = self._get_direction()
            self.long_side_len = 9500
            self.short_side_len = 3300
            self.raw_start = int(f['Analyses/Segmentation_00{}/Summary/segmentation'.format(mod_loc)].attrs['first_sample_template'])
            self.raw_step = int(f['Analyses/Basecall_1D_00{}/Summary/basecall_1d_template'.format(mod_loc)].attrs['block_stride'])
            self.guppy_trace = basecalled_template['Trace'][:]
            self.guppy_move = basecalled_template['Move'][:]
            self.readp2refplist = self.read_pos2ref_pos()

    def refp_from_readp(self, readpos):
        block = readpos // self.split_length
        if block >= len(self.readp2refplist):
            pos = -1
        elif self.readp2refplist[block][1] == '+':
            pos = self.readp2refplist[block][0] + (readpos % self.split_length)
        elif self.readp2refplist[block][1] == '-':
            pos = self.readp2refplist[block][0] - (readpos % self.split_length)
        else:
            pos = -1
        return pos

    def read_pos2ref_pos(self):
        readp2refplist = []
        for n, item in enumerate(self.sam_summary):
            if item[1] == '+':
                readp2refplist.append((int(item[2]), item[1]))
            elif item[1] == '-':
                readp2refplist.append((int(item[2]) + int(item[5]), item[1]))
        return readp2refplist

    def get_coord2raw(self):
        """Gets raw signal position from fastq sequence coordinate.

        Args:
            coord (int): coordinate in fastq
        Returns:
            coord2raw (dict):
        """
        coord2raw = {}
        m = 0
        for n, move in enumerate(self.guppy_move):
            if move:
                coord2raw[m] = self.raw_start + n * self.raw_step
                m += 1
        return coord2raw

    def get_methylated_bases(self, threshold, guppy_v):
        if guppy_v == 5:
            cpg = self.guppy_mbt[:,2]
        elif guppy_v == 4:
            cpg = self.guppy_mbt[:,3]
        dam = self.guppy_mbt[:,1]
        cpg_pos = np.where(cpg > threshold)
        dam_pos = np.where(dam > threshold)
        """
        Since the outputs of np.where is tuple we need to slice by [0].
        The reason the output is tuple is so that it can work on higher
        dimensional arrays.
        """
        return cpg_pos[0], dam_pos[0]

    def get_coding_methylations(self):
        cpg_pos = self._find_cpg()
        coding_cos = self.get_coding_coordinates()
        head_score = -1
        tail_score = -1
        for c in coding_cos:
            if c[0] < 15000:
                head_cpg_pos = [i for i in cpg_pos if c[0] < i < c[1]]
                head_score = np.mean(np.array(self.guppy_mbt[:,2])[head_cpg_pos])
            elif c[0] > 25000:
                tail_cpg_pos = [i for i in cpg_pos if c[0] < i < c[1]]
                tail_score = np.mean(np.array(self.guppy_mbt[:,2])[tail_cpg_pos])
        return head_score, tail_score

    def get_coordinate_methylations(self, start, end):
        cpg_pos = self._find_cpg()
        target_pos = [i for i in cpg_pos if start < i <end]
        return np.mean(np.array(self.guppy_mbt[:,2])[target_pos])

    def get_coding_coordinates(self):
        res = []
        former = []
        latter = []
        ret = []
        for n, ms in enumerate(self.sam_summary):
            pos = int(ms[2])
            if 0 < pos < self.long_side_len:
                former.append(n)
            elif self.long_side_len < pos < 13300:
                latter.append(n)
        if len(former) > 7:
            res.append(former)
        if len(latter) > 7:
            res.append(latter)
        for i in res:
            ret.append((min(i)*300, max(i)*300))
        return ret

    def find_short_deletions(self):
        in_codings = []
        for n, ms in enumerate(self.sam_summary[:,2].astype(np.int32)):
            if 0 < ms < 13100:
                in_codings.append((n, ms))
        if len(in_codings) < 30:
            return 0
        start = in_codings[0][0]
        end = in_codings[-1][0]
        if end - start > 50:
            return 0
        in_codings = in_codings[1:-1]
        prev_pos = in_codings[0][0]
        prev_coord  = in_codings[0][1]
        for pair in in_codings[1:]:
            posdiff = pair[0] - prev_pos
            coodiff = pair[1] - prev_coord
            prev_pos = pair[0]
            prev_coord = pair[1]
            if 300 * posdiff + 100 < abs(coodiff):
                return (pair[0], pair[1], posdiff, abs(coodiff), self.read_id)

    def get_IGS_coordinates(self):
        temp = []
        ret = ()
        for n, ms in enumerate(self.sam_summary):
            pos = int(ms[2])
            if 13500 < pos < 43500:
                temp.append(n)
        if len(temp) > 20:
            ret = (min(temp)*300, max(temp)*300)
        return ret

    def get_unmapped_section(self):
        print(self.read_id)
        print(calc_qual(self.quality))
        pos = self.sam_summary[:,2].astype(np.uint32)
        for n, p in enumerate(pos):
            if p == 0:
                um_qual = self.quality[n*self.split_length:(n+1)*self.split_length]
                if um_qual:
                    print(calc_qual(um_qual))

    def _find_cpg(self):
        cpg_sites = []
        upper_seq = self.seq.upper()
        for n in range(len(upper_seq) - 1):
            if upper_seq[n] == 'C':
                if upper_seq[n+1] == 'G':
                    cpg_sites.append(n)
        return cpg_sites

    def cpg_met_average(self, start, end):
        cpg_sites = self._find_cpg()
        in_cpgs = [i for i in cpg_sites if start < i < end]
        mod_scores = np.array(self.guppy_mbt[:,2])[in_cpgs]
        return np.mean(mod_scores) / 255

    def cpg_met_proportion(self, start, end):
        cpg_sites = self._find_cpg()
        in_cpgs = [i for i in cpg_sites if start < i < end]
        mod_scores = np.array(self.guppy_mbt[:,2])[in_cpgs]
        if len(mod_scores) == 0:
            met_pro = -1
        else:
            met_pro = len([i for i in mod_scores if i > 200])/len(mod_scores)
        return met_pro

    def cpg_scores(self):
        cpg_sites = self._find_cpg()
        mod_scores = np.array(self.guppy_mbt[:,2])[cpg_sites]
        return (cpg_sites, mod_scores)

    def get_coding_met_stats(self):
        cpg_func = self.cpg_met_average
        coding_coords = self.get_coding_coordinates()
        cpg_stats = []
        for start, end in coding_coords:
            met_stat = cpg_func(start, end)
            if met_stat < 0:
                continue
            cpg_stats.append(met_stat)
        return cpg_stats

    def get_IGS_met_stats(self):
        cpg_func = self.cpg_met_average
        start, end = self.get_IGS_coordinates()
        met_stat = cpg_func(start, end)
        return met_stat

    def __coding_coordinate(self):
        """
        get coding coordinate of read, not reference
        """
        region1 = self.long_side_len
        region2 = self.short_side_len
        length = len(self.seq)
        if self.direction == '+':
            a_s = 0
            a_e = region2
            b_s = self.length - region1
            b_e = self.length - 1
        elif self.direction == '-':
            a_s = 0
            a_e = region1
            b_s = self.length - region2
            b_e = self.length - 1
        return (a_s, a_e, b_s, b_e)

    def _truncate_fast5_file(self, outname, start, end, read_id=''):
        filename = self.fast5path
        length = end - start
        with h5py.File(filename, 'r') as f:
            with h5py.File(outname, 'w') as fw:
                f.copy('Raw', fw)
                f.copy('UniqueGlobalKey', fw)
                read_name = list(f['Raw/Reads'].keys())[0]
                if read_id:
                    fw['Raw/Reads/'+ read_name].attrs.modify('read_id', read_id.encode())
                fw['/'].attrs.create('file_version', 2, dtype='double')
                fw['Raw/Reads/' + read_name].attrs['duration'] = length
                signal = f['Raw/Reads/'+ read_name + '/Signal']
                del fw['Raw/Reads/'+ read_name + '/Signal']
                temp = fw['Raw/Reads/' + read_name].create_dataset('Signal', \
                        (length,), dtype='<i2', maxshape=(None,), \
                        chunks=(length,), compression=1)
                temp[:] = signal[start:end]
                fw.flush()

    def _get_raw_pos(self, coordinate):
        signal_length = len(self.raw)
        return int(signal_length/self.length * coordinate)

    def get_position_in_ref(self, pos):
        around_seq = self.seq[pos-50:pos+50]
        around_qual = self.quality[pos-50:pos+50]
        with open('temp_files/temp_fastq.fastq', 'w') as fw:
            fw.write('@temp\n' + around_seq + '\n+\n' + 
                     around_qual)
        subprocess.run('bwa mem -M -x ont2d -t 7 ' + self.ref_path + 
                       ' temp_files/temp_fastq.fastq > temp_files/'
                       'temp.sam',
                       shell=True, stdout=FNULL,
                       stderr=subprocess.STDOUT)
        with open('temp_files/temp.sam') as f:
            row = f.readlines()[2].strip().split()
            if row[2] == '*':
                return None
            cigar = row[5]
            ref_pos = int(row[3])
            c = Cigar(cigar)
            split_cigar = ''
            for i in c.items():
                split_cigar += i[0] * i[1]
            shift = 0
            current = 0
            for l in split_cigar:
                current += 1
                if l == 'I':
                    shift -= 1
                elif l == 'D':
                    shift += 1
                    current -= 1
                if current == 50:
                    ref_coordinate = ref_pos + 49 + shift
                    break
            return ref_coordinate

    def _find_coordinates(self, coords, temp_fastq_length=500):
        """Find exact rDNA coordinate in the given nanopore reads.

        Args:
            read (str): read sequence
            coord (int): target coordinate of rDNA

        Returns:
            result: list of [mapped coordinate, direction of mapping]
        """
        result = []
        reference_seq = ''
        with open(self.ref_path) as f:
            lines = f.readlines()[1:]
        for line in lines:
            reference_seq += line.strip()
        with open('temp_index/temp_index.fasta', 'w') as fw:
            fw.write('>{}\n{}'.format(self.read_id, self.seq))
        subprocess.run('bwa index temp_index/temp_index.fasta', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        for coord in coords:
            with open('temp_index/coordinate_rDNA.fastq', 'w') as fw:
                fw.write('>temp\n{}\n+\n{}\n'.format(reference_seq[coord-1:coord+temp_fastq_length-1], 'J' * temp_fastq_length))
            # with -a option, multiple hits are more clearly shown
            utilities.bwa_mapping('temp_index/temp_index.fasta', 'temp_index/coordinate_rDNA.fastq', 'temp_index/temp_sam4coord.sam', multi=True)
            with open('temp_index/temp_sam4coord.sam') as samf:
                map_result = samf.readlines()[2:]
            for mapping in map_result:
                row = mapping.strip().split()
                AS = int(mapping.strip().split('AS:i:')[1].split()[0])
                flag = int(row[1])
                if utilities.easy_flag(flag, 16) != 1:
                    direction = '+'
                else:
                    direction = '-'
                mapped_coord = int(row[3])
                if AS > 0.3 * temp_fastq_length:
                    result.append([coord, mapped_coord, direction])
        return result

    def large_scale_by_mapping(self):
        mapped_results = self._find_coordinates([2000, 3500, 5000, 6500, 8000, 9500, 11000, 18000, 31500, 33500, 35000, 38000, 41000])
        mapped_sorted = sorted(mapped_results, key=lambda res: res[1])
        for n in range(len(mapped_sorted) - 1):
            if mapped_sorted[n][2] == '+':
                read_diff = mapped_sorted[n+1][1] - mapped_sorted[n][1]
                ref_diff = (mapped_sorted[n+1][0] - mapped_sorted[n][0]) % 42999
                diff_from_ref = abs(ref_diff - read_diff)
                if mapped_sorted[n+1][0] > 25000 > mapped_sorted[n][0]:
                    if diff_from_ref > 10000:
                        return 1
                elif mapped_sorted[n+1][0] > 14000 > mapped_sorted[n][0]:
                    if diff_from_ref > 6000:
                        return 1
                else:
                    if diff_from_ref > 3000:
                        return 1
            else:
                read_diff = mapped_sorted[n+1][1] - mapped_sorted[n][1]
                ref_diff = (mapped_sorted[n][0] - mapped_sorted[n+1][0]) % 42999
                diff_from_ref = abs(ref_diff - read_diff)
                if mapped_sorted[n][0] > 25000 > mapped_sorted[n+1][0]:
                    if diff_from_ref > 10000:
                        return 1
                elif mapped_sorted[n][0] > 14000 > mapped_sorted[n+1][0]:
                    if diff_from_ref > 6000:
                        return 1
                else:
                    if diff_from_ref > 3000:
                        return 1
        return 0

    def truncate_and_separate_by_dam(self, start, end, fast5_dir, control_dir, offset):
        #separate reads by the dam status within (start, end) reference coordinate
        coding_coords = self.__coding_coordinate()
        if end < 9200 or start > 20000:
            if self.direction == '+':
                target_side = 'tail'
            else:
                target_side = 'head'
        elif start > 10000:
            if self.direction == '+':
                target_side = 'head'
            else:
                target_side = 'tail'
        else:
            return None
        raw_coords = [self._get_raw_pos(i) for i in coding_coords]
        if target_side == 'head':
            raw_pos = (raw_coords[0], raw_coords[1])
        else:
            raw_pos = (raw_coords[2], raw_coords[3])
        dam_met_pos = self.get_methylated_bases(180)[1]
        ref_dam_pos = [self.get_position_in_ref(i) for i in dam_met_pos]
        ref_dam_pos = [i for i in ref_dam_pos if i != None]
        out_file = '{}/{}.fast5'
        if not ref_dam_pos:
            self._truncate_fast5_file(out_file.format(control_dir, self.read_id),
                                      raw_pos[0], raw_pos[1])
        elif any(((start + offset) % 42999 < i < (end + offset) % 42999 for i in ref_dam_pos)):
            self._truncate_fast5_file(out_file.format(fast5_dir, self.read_id),
                                      raw_pos[0], raw_pos[1])
        else:
            self._truncate_fast5_file(out_file.format(control_dir, self.read_id),
                                      raw_pos[0], raw_pos[1])

    def truncate_and_separate_by_cpg(self, guppy_v, dir_basename='', only_head='no'):
        dirs = [dir_basename + '_head_met',
                dir_basename + '_head_unmet',
                dir_basename + '_tail_met',
                dir_basename + '_tail_unmet']
        for direc in dirs:
            if not os.path.exists(direc):
                if only_head == 'yes':
                    if 'head' not in direc:
                        continue
                os.mkdir(direc)
        coding_coords = self.__coding_coordinate()
        raw_pos = [self._get_raw_pos(i) for i in coding_coords]
        cpg_table = self.guppy_mbt[:,3]
        head_met_ratio, tail_met_ratio = self.get_coding_methylation_pros('cpg', guppy_v)
        out_file = '{}/{}.fast5'
        threshold = 0.15
        if head_met_ratio < threshold:
            self._truncate_fast5_file(out_file.format(dirs[1], self.read_id),
                    raw_pos[0], raw_pos[1], read_id=self.read_id[:-2]+'01')
        elif head_met_ratio > threshold:
            self._truncate_fast5_file(out_file.format(dirs[0], self.read_id),
                                      raw_pos[0], raw_pos[1], read_id=self.read_id[:-2]+'01')
        if only_head == 'yes':
            return head_met_ratio, -1 
        if tail_met_ratio < threshold:
            self._truncate_fast5_file(out_file.format(dirs[3], self.read_id),
                                      raw_pos[2], raw_pos[3], read_id=self.read_id[:-2]+'02')
        elif tail_met_ratio > threshold:
            self._truncate_fast5_file(out_file.format(dirs[2], self.read_id),
                                      raw_pos[2], raw_pos[3], read_id=self.read_id[:-2]+'02')
        return head_met_ratio, tail_met_ratio

    def cpg_methylation_average(self, window, caller='default'):
        """
        Divide the region by window, calculate the average number of methylated CpGs.
        """
        x = []
        y = []
        cpg_sites = self._find_cpg()
        if caller == 'default':
            mbt_site = 2
        elif caller == 'rerio':
            mbt_site = 3
        for n in range(0, int(self.length/window)):
            x.append((n + 0.5) * window)
            in_cpgs = [i for i in cpg_sites if n * window < i < (n+1) * window]
            mod_scores = np.array(self.guppy_mbt[:,mbt_site])[in_cpgs]
            y.append(np.mean(mod_scores) * 40)
        return x, y

    def _find_a(self):
        a_sites = []
        upper_seq = self.seq.upper()
        for n in range(len(upper_seq) - 1):
            if upper_seq[n] == 'A':
                a_sites.append(n)
        return a_sites

    def m6a_methylation_average(self, window):
        x = []
        y = []
        a_sites = self._find_a()
        for n in range(0, int(self.length/window)):
            x.append((n + 0.5) * window)
            in_a = [i for i in a_sites if n * window < i < (n+1) * window]
            mod_scores = np.array(self.guppy_mbt[:,1])[in_a]
            y.append(np.mean(mod_scores) * 40)
        return x, y

    def is_this_end_to_end(self, gRNA='28S'):
        poss = [int(n) for n in self.sam_summary[:,2]]
        start, end = 0, 0
        for pos in poss[:3]:
            if pos != 0:
                start = pos
                break
        for pos in reversed(poss[-3:]):
            if pos != 0:
                end = pos
                break
        if not (start and end):
            return 0
        if gRNA == '28S':
            left = 7000
            right = 11000
        elif gRNA == '18S':
            left = 3000
            right = 7000
        if left < start < right and left < end < right:
            return 1
        else:
            return 0

    def plot_structure(self, savedir, mouse='no', met='yes', eps='no', adenine='no'):
        if not os.path.exists(savedir):
            os.mkdir(savedir)
        lc = utilities.plot_read_structure(self.read_id, self.split_length, self.sam_summary, mouse=mouse)
        fig = plt.figure()
        plt.subplots_adjust(left=0.2)
        ax = fig.add_subplot()
        met_window = 200
        """
        for n, item in enumerate(self.guppy_mbt[:,1]):
            if item > 120:
                ax.bar(n, item * 10000/255, width=100, color='red', zorder=3)
        """
        if met == 'yes':
            x, y = self.cpg_methylation_average(met_window)
            ax.bar(x, y, width=met_window, color = 'magenta', zorder=0)
        elif met == 'rerio':
            x, y = self.cpg_methylation_average(met_window, caller='rerio')
            ax.bar(x, y, width=met_window, color = 'magenta', zorder=0, alpha=0.3)
        if adenine == 'yes':
            x, y = self.m6a_methylation_average(met_window)
            ax.bar(x, y, width=met_window, color = 'blue', zorder=0, alpha=0.3)
        bin1 = 200
        x = []
        y = []
        for n in range(0, self.length - bin1, bin1):
            x.append(n)
            y.append(np.mean([ord(i) - 33 for i in self.quality[n:n+bin1]]) * 200 - 10000)
        ax.plot(x, y)
        ax.add_collection(lc)
        if mouse == 'yes':
            ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000, 50000))
            ax.set_yticklabels(('UM', '0', '10k', '20k', '30k', '40k', '50k'), fontsize=20)
        else:
            ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000))
            ax.set_yticklabels(('UM', '0', '10k', '20k', '30k', '40k'), fontsize=20)
        ax.set_xticks((0, 10000, 20000, 30000, 40000, 50000))
        ax.set_xticklabels(('0', '10k', '20k', '30k', '40k', '50k'), fontsize=20)
        ax.set_title('Mean quality: {:.1f}'.format(self.ave_quality))
        ax.autoscale()
        ax.set_ylim([-12000, 47000])
        if savedir[-1] == '/':
            savedir = savedir[:-1]
        if eps == 'yes':
            plt.savefig(savedir + '/' + self.read_id + '.eps')
        else:
            plt.savefig(savedir + '/' + self.read_id + '.png', dpi=300)
        plt.close()


def met_level_per_coding(savename='temp.png'):
    ref = 'rDNA_index/musRibosomal.fa'
    met_pros = []
    for n, f5 in enumerate(glob.glob('B6BM_rDNA_newbc/**/*.fast5')):
        r = Read(f5, ref, mod_loc=2)
        if r._is_this_healthy_rDNA() and r.length > 15000:
            mets = r.get_coding_met_stats()
            met_pros += mets
            print(mets)
    """
    if len(met_pros):
        print(f'{sample}\tnonmet: {len([i for i in met_pros if i < 50])/len(met_pros)}\tmet: {len([i for i in met_pros if i > 200])/len(met_pros)}')
    """
    savename = f'B6BM.png'
    #weights = np.ones_like(met_pros)/float(len(met_pros))
    plt.xlim(0, 1)
    plt.hist(met_pros, bins=20, range=(0, 1))
    plt.savefig(savename, dpi=300)
    plt.close()


def plot_dam_positions(reads):
    pos_in_ref = []
    for r in reads:
        cpg_pos, dam_pos = r.get_methylated_bases(180)
        for pos in dam_pos:
            ref_pos = r.get_position_in_ref(pos)
            if ref_pos:
                pos_in_ref.append(ref_pos)
    with open('dam_distribution.txt', 'w') as fw:
        fw.write(str(sorted(pos_in_ref)))
    plt.hist(pos_in_ref, bins=1000)
    plt.savefig('dam_distribution.png', dpi=500)
    return pos_in_ref


def basecalling(dirname, out_name='', ver=5):
    if dirname[-1] == '/':
        dirname = dirname[:-1]
    if not out_name:
        out_name = dirname + '_bc'
    if ver == 5:
        string = '~/Softwares/ont-guppy_5.0.16_linux64/ont-guppy/bin/guppy_basecaller -i {} -s {} -c dna_r9.4.1_450bps_modbases_5mc_hac_prom.cfg --device cuda:0 --fast5_out --recursive'.format(dirname, out_name)
    elif ver == 4:
        string = '~/Softwares/ont-guppy_4.2.2_linux64/ont-guppy/bin/guppy_basecaller -i {} -s {} -c dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac_prom.cfg --device cuda:0 --fast5_out --recursive'.format(dirname, out_name)
    elif ver == 6:
        string = '~/Softwares/ont-guppy_6.4.6_linux64/ont-guppy/bin/guppy_basecaller -i {} -s {} -c dna_r9.4.1_450bps_modbases_5mc_cg_hac.cfg --device cuda:0 --recursive'.format(dirname, out_name)
    subprocess.run(string, shell=True)


def resquiggle(bc_dirname, ref):
    string = 'tombo resquiggle {} {} --processes 10 --num-most-common-errors 5 --overwrite'.format(bc_dirname, ref)
    subprocess.run(string, shell=True)


def tombo_plot(coordinates, offset, bc_dirname, bc_alt_dirname, pdf_filename):
    """
    coordinates, bc_dirname, bc_alt_dirname should be list
    """
    coord_str = ''
    for co in (np.array(coordinates) + offset) % 42999:
        coord_str += ' gi\|555853\|gb\|U13369.1\|HSU13369:' + str(co)
    bc_str = ''
    for dirname in bc_dirname:
        bc_str += dirname + ' '
    bc_alt_str = ''
    for dirname in bc_alt_dirname:
        bc_alt_str += dirname + ' '
        string = 'tombo plot genome_locations --plot-standard-model --fast5-basedirs {} --control-fast5-basedirs {} --genome-locations {} --pdf-filename {}'.format(bc_str, bc_alt_str, coord_str, pdf_filename)
    subprocess.run(string, shell=True)


def separate_by_dam_and_plot_main(reads, coordinate, offset, ref):
    fast5_dir = 'dam_fast5s'
    bc_dir = fast5_dir + '_basecalled'
    control_dir = 'dam_control_fast5s'
    ctl_bc_dir = control_dir + '_basecalled'
    if os.path.exists(fast5_dir):
        shutil.rmtree(fast5_dir)
    if os.path.exists(control_dir):
        shutil.rmtree(control_dir)
    if os.path.exists(bc_dir):
        shutil.rmtree(bc_dir)
    if os.path.exists(ctl_bc_dir):
        shutil.rmtree(ctl_bc_dir)
    os.mkdir(fast5_dir)
    os.mkdir(control_dir)
    for r in reads:
        r.truncate_and_separate_by_dam(coordinate - 5, coordinate + 5, fast5_dir, control_dir, offset)
    basecalling(fast5_dir, bc_dir)
    basecalling(control_dir, ctl_bc_dir)
    shutil.rmtree(fast5_dir)
    shutil.rmtree(control_dir)
    resquiggle(bc_dir, ref)
    resquiggle(ctl_bc_dir, ref)
    tombo_plot(coordinate, offset, bc_dir, ctl_bc_dir, 
            pdf_filename='dam_' + str(coordinate) + '.pdf')


def separate_by_cpg_align_and_squiggle(reads, dir_basename='', guppy_v=5, only_head='no'):
    dirs = [dir_basename + 'head_met',
            dir_basename + 'head_unmet',
            dir_basename + 'tail_met',
            dir_basename + 'tail_unmet']
    called_dirs = [i + '_basecalled' for i in dirs]
    for directory in dirs:
        if os.path.exists(directory):
            shutil.rmtree(directory)
    for directory in called_dirs:
        if os.path.exists(directory):
            shutil.rmtree(directory)
    for directory in dirs:
        if only_head == 'yes':
            if 'head' not in directory:
                continue
        os.mkdir(directory)
    hs = []
    for r in reads:
        head, tail = r.truncate_and_separate_by_cpg(guppy_v, dir_basename, only_head)
        hs.append(head)
    plt.hist(hs, bins=10)
    plt.show()
    for in_dir, out_dir in zip(dirs, called_dirs):
        if only_head == 'yes':
            if 'head' not in in_dir:
                continue
        basecalling(in_dir, out_dir)
        shutil.rmtree(in_dir)
    for called_dir in called_dirs:
        if only_head == 'yes':
            if 'head' not in called_dir:
                continue
        resquiggle(called_dir, ref)


def truncate_head(reads, strand, basename):
    if strand == '-':
        dirname = 'head_minus_{}'.format(basename)
        called_dir = 'head_minus_{}_bc'.format(basename)
    elif strand == '+':
        dirname = 'head_plus_{}'.format(basename)
        called_dir = 'head_plus_{}_bc'.format(basename)
    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    os.mkdir(dirname)
    if strand == '+':
        length = 3000
    elif strand == '-':
        length = 11000
    for r in reads:
        if r.direction == strand and r.length > length:
            start = 0
            end = r.get_coord2raw()[length]
            out_file = '{}/{}.fast5'
            r._truncate_fast5_file(out_file.format(dirname, r.read_id), start, end)
    basecalling(dirname, called_dir)
    shutil.rmtree(dirname)
    resquiggle(called_dir, 'rDNA_index/rDNA_for_cas9_2.fasta')


def plot_by_cpg(reads, coordinates, offset, ref, make_fast5s, dir_basename=''):
    dirs = [dir_basename + 'head_met',
            dir_basename + 'head_unmet',
            dir_basename + 'tail_met',
            dir_basename + 'tail_unmet']
    called_dirs = [i + '_basecalled' for i in dirs]
    if make_fast5s == 'yes':
        separate_by_cpg_align_and_squiggle(reads, dir_basename)
    bc_dirs = [i for i in called_dirs if '_met' in i]
    ctl_bc_dirs = [i for i in called_dirs if '_unmet' in i]
    coordinates = np.array(coordinates)
    tombo_plot(coordinates, bc_dirs, ctl_bc_dirs, 
            pdf_filename='1020cpg_{}-{}.pdf'.format(coordinates[0], coordinates[1]))


def extract_rDNA_and_make_reads(rDNA_name_file, fast5_dir):
    if fast5_dir[-1] != '/':
        fast5_dir = fast5_dir + '/'
    rDNA_read_ids = []
    with open(rDNA_name_file) as f:
        for line in f:
            rDNA_read_ids.append(line.split()[0])
    fast5files = glob.glob(fast5_dir)
    rDNA_files = [fast5_dir+ i + '.fast5' for i in rDNA_read_ids]
    reads = []
    n = 0
    for rfile in rDNA_files:
        reads.append(Read(rfile, ref))
    return reads


def nanopore_reads_to_fig(in_dir, base_name, ref):
    if os.path.exists(base_name + '_plot'):
        shutil.rmtree(base_name + '_plot')
    os.mkdir(base_name + '_plot')
    subprocess.run('multi_to_single_fast5 -t 6 -i ' + in_dir + ' -s ' + base_name + '_single', shell=True)
    basecalling(base_name + '_single', base_name + '_bc_fast5s')
    shutil.rmtree(base_name + '_single')
    fast5s = glob.glob(base_name + '_bc_fast5s/workspace/**/*.fast5')
    ret_str = ''
    for f in fast5s:
        read = Read(f, ref)
        if read._is_this_healthy_rDNA():
            if read.is_this_end_to_end():
                ret_str += '\n'.join(read.seq) + '\n'
                read.plot_structure(base_name + '_plot')
    """
    with open(base_name + '_rDNAs.fastq', 'w') as fw:
        fw.write(ret_str)
    """

def find_abnormal_to_fig(in_dir, ref):
    fast5s = glob.glob(in_dir + 'workspace/**/*.fast5')
    basename = in_dir.split('/')[-2]
    out_dir = 'abnormal_reads2_nomet/' + basename
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.mkdir(out_dir)
    print(out_dir)
    for f in fast5s:
        read = Read(f, ref)
        if read._is_this_healthy_rDNA():
            if read.is_this_end_to_end():
                if read.large_scale_by_mapping():
                    read.plot_structure(out_dir, met='no')


def extract_rDNA_fq_from_f5():
    for bc_dir in glob.glob('/mnt/data/nanopore_cas9_bc/230214_Balbc_young_bc/'):
        basename = bc_dir.split('/')[-2].split('_bc')[0]
        fast5s = glob.glob(bc_dir + 'workspace/**/*.fast5')
        ret_str = ''
        for f in fast5s:
            read = Read(f, 'rDNA_index/musRibosomal.fa')
            if read.is_this_end_to_end() and read._is_this_healthy_rDNA():
                ret_str += '\n'.join(read.fastq_str) + '\n'
        with open('fastqs/' + basename + '_rDNAs.fastq', 'w') as fw:
            fw.write(ret_str)


def met_stats_hist():
    ref = 'rDNA_index/humRibosomal.fa'
    mets = []
    for fast5 in glob.glob('210215_iPS_bc/**/*.fast5', recursive=True):
        r = Read(fast5, ref)
        if r.length > 60000:
            continue
        codings = r.get_coding_coordinates()
        if codings:
            met_stats = r.get_coding_met_stats()
            mets += met_stats
    plt.hist(mets, bins=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig('B7iPS_met.eps')


def met_pickle():
    ref = 'rDNA_index/musRibosomal.fa'
    n = 0
    mets = []
    for fast5 in glob.glob('/media/owner/809f403f-5b66-4d70-be53-a585528402c5/nanopore_cas9_bc/*Balbc_young_bc/**/*.fast5', recursive=True):
        print(fast5)
        r = Read(fast5, ref, mod_loc=0)
        if not r._is_this_healthy_rDNA():
            continue
        if r.length > 70000:
            continue
        n += 1
        codings = r.get_coding_coordinates()
        if codings:
            met_stats = r.get_coding_met_stats()
            mets += met_stats
        if len(mets) > 300:
            break
    plt.hist(mets, bins=20, range=(0, 1.0))
    plt.xlabel('proportion of 45S methylation')
    plt.show()


def to_single_f5():
    for f5 in glob.glob('/mnt/data/nanopore_cas9_bc/210225_B6BM/workspace/B6BM/20210225_0758_MN32877_AFT144_44417472/'):
        subprocess.run('multi_to_single_fast5 -s {} -i {} -t 8 --recursive'.format('/mnt/data/nanopore_cas9_bc/210225_B6BM_single', f5, shell=True))


def rerio_analysis_28SPCR():
    ref = 'rDNA_index/humRibosomal2.fa'
    g1 = []
    g2 = []
    for n, f5 in enumerate(glob.glob('g1_3_bc/**/*.fast5', recursive=True)):
        r = Read(f5, ref)
        mbt = np.array(r.guppy_mbt)
        Y = mbt[:, 1]
        g1.append(np.mean(Y))
    for n, f5 in enumerate(glob.glob('g2_3_bc/**/*.fast5', recursive=True)):
        r = Read(f5, ref)
        mbt = np.array(r.guppy_mbt)
        Y = mbt[:, 1]
        g2.append(np.mean(Y))
    plt.hist(g1, color='red', alpha=0.3)
    plt.hist(g2, color='blue', alpha=0.3)
    plt.show()


def rerio_analysis_RPA43():
    ref = 'rDNA_index/humRibosomal2.fa'
    for n, f5 in enumerate(glob.glob('293_EcoGII_bc2/**/*.fast5', recursive=True)):
        r = Read(f5, ref, mod_loc='2')
        if r._is_this_healthy_rDNA() and r.length > 5000:
            r.plot_structure('temp_figs/293_ME', met='rerio', adenine='yes')


def prom_coding_met_correlation():
    ref = 'rDNA_index/musRibosomal.fa'
    n = 0
    prom_coding_mets = []
    for fast5 in glob.glob('/media/owner/809f403f-5b66-4d70-be53-a585528402c5/nanopore_cas9_bc/230216_Balbc_old_bc/workspace/**/*.fast5', recursive=True):
        try:
            r = Read(fast5, ref)
            if not r._is_this_healthy_rDNA():
                continue
            if r.length > 60000:
                continue
            coding_coords = r.get_coding_coordinates()
            if not coding_coords:
                continue
            """
            prom_coords = r._find_coordinates([44530], temp_fastq_length=300)
            for prom in prom_coords:
                if prom[2] == '+':
                    for coding in coding_coords:
                        if prom[1] - 1500 < coding[0] < prom[1] + 1500:
                            prom_coding_mets.append((r.cpg_met_average(prom[1] + 120, prom[1] + 340), r.cpg_met_average(coding[0], coding[1])))
                else:
                    for coding in coding_coords:
                        if prom[1] - 1500 < coding[1] < prom[1] + 1500:
                            prom_coding_mets.append((r.cpg_met_average(prom[1] - 340, prom[1] - 120), r.cpg_met_average(coding[0], coding[1])))
            """
        except:
            continue
    #prom_coding_mets = pd.read_pickle('Y10026_prom_coding.pkl')
    #pd.to_pickle(prom_coding_mets, 'Y10026_prom_coding.pkl')
    np_mets = np.array(prom_coding_mets)
    x = np_mets[:,0]
    y = np_mets[:,1]
    plt.scatter(x, y)
    plt.xlabel('promoter methylation', fontsize=15)
    plt.ylabel('gene body methylation', fontsize=15)
    plt.xlim((0, 1))
    plt.ylim((0, 1))
    plt.show()
    result = scipy.stats.linregress(np_mets[:,0], np_mets[:,1])
    print(result)


def f5_to_plot(date):
    dirs = glob.glob('/var/lib/minknow/data/*')
    ref = 'rDNA_index/humRibosomal2.fa'
    """
    for d in dirs:
        if date in d:
            data_d = d
            break
    fast5_dir = glob.glob('{}/*/*/fast5_pass'.format(data_d))[0]
    print(fast5_dir)
    subprocess.run('multi_to_single_fast5 -t 8 -i {} -s {}'.format(fast5_dir, date), shell=True)
    basecalling(date)
    """
    if os.path.exists('temp_figs/{}'.format(date)):
        shutil.rmtree('temp_figs/{}'.format(date))
    if os.path.exists('temp_figs/{}_full'.format(date)):
        shutil.rmtree('temp_figs/{}_full'.format(date))
    for n, f5 in enumerate(glob.glob('{}_bc/**/*.fast5'.format(date), recursive=True)):
        if n % 500 == 0:
            print(n)
        r = Read(f5, ref, mod_loc=0)
        if r._is_this_healthy_rDNA():
            if r.is_this_end_to_end(gRNA='18S'):
                r.plot_structure('temp_figs/{}_full'.format(date))
            else:
                r.plot_structure('temp_figs/{}'.format(date))


def multi_to_single(dir_loc):
    if dir_loc[-1] == '/':
        dir_loc = dir_loc[:-1]
    out_dir = '/mnt/data/nanopore_cas9_bc/{}'.format(dir_loc.split('/')[-1])
    subprocess.run('multi_to_single_fast5 --recursive -i {} -s {} -t 8'.format(dir_loc, out_dir), shell=True)
    return out_dir


def just_visualize(dir_loc, fig_dir='temp_figs', quality_threshold=15, cr=False, bc=False, mouse=False, already_bc=False):
    if dir_loc[-1] == '/':
        dir_loc = dir_loc[:-1]
    dir_name = dir_loc.split('/')[-1]
    ref = 'rDNA_index/humRibosomal2.fa'
    if mouse:
        ref = 'rDNA_index/musRibosomal.fa'
    out_dir = '/mnt/data/nanopore_cas9_bc/{}'.format(dir_loc.split('/')[-1])
    if already_bc == False:
        if bc:
            bc_dir = out_dir + '_bc'
        else:
            bc_dir = out_dir
        if not os.path.exists(bc_dir):
            multi_to_single(dir_loc)
            if bc:
                basecalling(out_dir)
                shutil.rmtree(out_dir)
    else:
        bc_dir = dir_loc
    if cr:
        for f5 in glob.glob('{}/**/*.fast5'.format(bc_dir), recursive=True):
            r = Read(f5, ref, mod_loc=0)
            mean_coord = np.mean(r.sam_summary[:, 2][10:-10].astype("i4"))
            if r._is_this_healthy_rDNA() and 19000 > r.length > 13000 and mean_coord < 15000 and r.ave_quality > quality_threshold:
                r.plot_structure('{}/{}'.format(fig_dir, dir_name))
    else:
        for f5 in glob.glob('{}/**/*.fast5'.format(bc_dir), recursive=True):
            r = Read(f5, ref, mod_loc=0)
            if r._is_this_healthy_rDNA():
                r.plot_structure('{}/{}'.format(fig_dir, dir_name))


def find_abnormal_reads(sample_name):
    ref = 'rDNA_index/humRibosomal2.fa'
    Ref = main_analysis_fastq.Reference(ref)
    n = 0
    for exp in glob.glob('/mnt/data/nanopore_cas9_bc/*/'):
        if sample_name in exp:
            for f5 in glob.glob(exp + '**/*.fast5', recursive=True):
                r = Read(f5, ref, mod_loc=0)
                if r.is_this_end_to_end(gRNA='18S') and r._is_this_healthy_rDNA():
                    if r.large_scale_by_mapping():
                        r.plot_structure('abnormal_reads/{}'.format(sample_name))
                    n += 1
                    fq = main_analysis_fastq.Fastqread(r.fastq_str, Ref)
                    if fq.get_leaps():
                        r.plot_structure('abnormal_reads/{}'.format(sample_name))


def find_short_deletions_from_sample(sample_name):
    ref = 'rDNA_index/humRibosomal2.fa'
    for d in glob.glob('/mnt/data/nanopore_cas9_bc/**/'):
        if sample_name in d:
            for f5 in glob.glob(d + '**/*.fast5', recursive=True):
                r = Read(f5, ref, mod_loc=0)
                del_res = r.find_short_deletions()
                if del_res and r.ave_quality > 15:
                    r.plot_structure('potential_del/{}'.format(sample_name))


def average_pattern(directory, ref, mod_loc=1):
    scores = []
    nm_scores = []
    m_scores = []
    ref_seq = ''
    with open(ref) as f:
        for line in f:
            if '>' not in line:
                ref_seq += line.strip()
    rDNA_len = len(ref_seq)
    met_per_cod = []
    hts = []
    coding_mets = []
    for f5 in glob.glob(f'{directory}/**/*.fast5', recursive=True):
        r = Read(f5, ref, mod_loc)
        if not r._is_this_healthy_rDNA() or r.length < 15000:
            continue
        cpg_pos, mod_score = r.cpg_scores()
        cpg_ref_pos = [r.refp_from_readp(i) for i in cpg_pos]
        head, tail = r.get_coding_methylations()
        if head != -1:
            coding_mets.append(head)
        if tail != -1:
            coding_mets.append(tail)
        if r.length > 30000:
            hts.append((r.read_id, head, tail))
        for p, s in zip(cpg_ref_pos, mod_score):
            if p == -1:
                continue
            scores.append((p, s))
            if 0 < head/255 < 0.1 and 0 < tail/255 < 0.1:
                nm_scores.append((p, s))
            elif head/255 > 0.3 and tail/255 > 0.3:
                m_scores.append((p, s))
    """
    plt.hist(coding_mets, bins=20)
    plt.show()
    with open('head_tail_met.txt', 'w') as fw:
        for rid, h, t in sorted(hts):
            fw.write('{}\t{}\t{}\n'.format(rid, h, t))
    quit()
    """
    met_scores = collections.defaultdict(list)
    nm_met_scores = collections.defaultdict(list)
    met_met_scores = collections.defaultdict(list)
    split_bin = 300
    for item in scores:
        met_scores[(item[0]%rDNA_len)//split_bin].append(item[1])
    for item in nm_scores:
        nm_met_scores[(item[0]%rDNA_len)//split_bin].append(item[1])
    for item in m_scores:
        met_met_scores[(item[0]%rDNA_len)//split_bin].append(item[1])

    return (met_scores, nm_met_scores, met_met_scores, split_bin)
    """
    for co_bin, mets in nm_met_scores.items():
        plt.bar(co_bin * split_bin, np.mean(mets)/255, width=split_bin, color='blue')
    plt.show()
    plt.close()
    for co_bin, mets in met_met_scores.items():
        #mets = [i for i in mets if i != 0]
        plt.bar(co_bin * split_bin, np.mean(mets)/255, width=split_bin, color='blue')
    plt.show()
    plt.close()
    met_binned = []
    with open('average_m6a_scores.txt') as f:
        for line in f:
            met_binned.append(line.strip().split())
    for co, met in met_binned:
        print((-1 * int(co) + 6900) % 9100)
        plt.bar((-1 * int(co) + 6900) % 9100, float(met), width=100, color='blue')
    for co_bin, mets in met_scores.items():
        print(co_bin * 100, np.mean(mets))
        plt.bar(co_bin * 100, np.mean(mets), width=100, color='blue')
    plt.show()
    """


def neighboring_mets():
    ref = 'rDNA_index/musRibosomal.fa'
    n = 0
    mets = []
    x = []
    y = []
    for n, fast5 in enumerate(glob.glob('/mnt/data/nanopore_cas9_bc/210225_B6BM_single/**/*.fast5', recursive=True)):
        r = Read(fast5, ref, mod_loc=1)
        if not r._is_this_healthy_rDNA():
            continue
        if r.length > 60000:
            continue
        n += 1
        codings = r.get_coding_coordinates()
        if codings:
            met_stats = r.get_coding_met_stats()
            if len(met_stats) == 2:
                x.append(met_stats[0])
                y.append(met_stats[1])
    plt.scatter(x, y)
    plt.show()


def clai_site():
    ref = 'rDNA_index/musRibosomal.fa'
    clai_site = 'CGCGTCGTTGCTCACTCTTAGATCGATGTGGTGCTCCGGAGTTCTCTTCGGGCCAGGGCCAAGCCGCGCCAGGCGAGGGACGGACATTCATGGCGAATGGCGGCCGCTCTTCTCGTTCTGCCAGCGGGCCCTCGTCTCTCCACCCCATCCGTCTGCCGGTGGTGTGTGGAAGGCAGGGGTGCGGCTCTCCGGCCCGACGCTGCCCCGCGCGCACTTTTCTCAGTGGTTCGCGTGGTCCTTGTGGATGTGTGAGGCGCCCGGTTGTGCCCTCACGTGTTTCACTTTGGTCGTGTCTCGCTTGACCATGTTCCCAGAGTCGGTGGATGTGGCCGGTGGCGTTGCATACCCTTCCCGTCTGGTGTGTGCACGCGCTGTTTCTTGTAAGCGTCGAGGTGCTCCTGGAGCGTTCCAGGTTTGTCTCCTAGGTGCCTGCT'
    temp_fq = 'temp_files/temp_fq_mouse.fq'
    with open(temp_fq, 'w') as fw:
        fw.write(f'@ClaI_site\n{clai_site}\n+\n{"J"*len(clai_site)}')
    clai_met_aves = []
    for n, fast5 in enumerate(glob.glob('/mnt/data/nanopore_cas9_bc/210225_B6BM_single/**/*.fast5', recursive=True)):
        r = Read(fast5, ref, mod_loc=1)
        if not r._is_this_healthy_rDNA():
            continue
        seq = r.seq
        temp_fa = f'temp_files/temp_{r.read_id}.fa'
        temp_out = f'temp_files/temp_{r.read_id}.sam'
        with open(temp_fa, 'w') as fw:
            fw.write(f'>temp\n{seq}')
        subprocess.run(f'minimap2 -ax map-ont {temp_fa} {temp_fq} > {temp_out}', shell=True)
        with pysam.AlignmentFile(temp_out) as af:
            for read in af:
                start = read.reference_start
                if start == -1:
                    continue
                end = start+500
                clai_met_aves.append(r.get_coordinate_methylations(start, end))
        os.remove(temp_fa)
        os.remove(temp_out)
    plt.hist(clai_met_aves, bins=30)
    plt.show()


if __name__ == '__main__':
    for f5 in glob.glob("/mnt//data//nanopore_cas9_bc/221124_Y8585fN_CR_bc/workspace/**/*.fast5", recursive=True):
        r = Read(f5, config.RDNA_REF_HUMAN)
        if r._is_this_healthy_rDNA() and r.length > 10000:
            print(f5)
            r.plot_structure("temp_figs")
    #average_pattern()
    #basecalling('B6BM_rDNA_f5s/', 'B6BM_rDNA_newbc')
    #basecalling('/mnt/data/nanopore_cas9_bc/210115_WS/', '/mnt/data/nanopore_cas9_bc/210115_WS_bc', ver=5)
    #met_level_per_coding()
    #clai_site()
    quit()
    #just_visualize('/mnt/data/Minion_data/240702_RP11_Cas9/', mouse=False, bc=True)
    #met_pickle()
    #basecalling('/media/owner/809f403f-5b66-4d70-be53-a585528402c5/Minion_220607/210215_B7iPS', 'B7iPS_bc', ver=5)
    find_short_deletions_from_sample('240702_RP11_Cas9')
    quit()
    for d in glob.glob('/mnt/data/nanopore_cas9_bc/**_CR_bc/'):
        just_visualize(d, fig_dir='brain_visualized', cr=True, bc=True, mouse=False, already_bc=True)
    ref = 'rDNA_index/humRibosomal2.fa'
    quit()
    #basecalling('220419')
    truncate_head(reads, '+', '293_EcoGII')
    truncate_head(reads, '-', '293_EcoGII')
    quit()
    #separate_by_cpg_align_and_squiggle(reads, 'RPA43_plus_new', guppy_v=5, only_head='yes')
    """
    CpG plotting requires a list or tuple of coordinates.
    dam plotting requires a scalar coordinate.
    """
