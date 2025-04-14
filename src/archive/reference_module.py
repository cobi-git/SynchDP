import sys
import os
import shutil
import glob
import pickle
import numpy as np
import pandas as pd
from dataclasses import dataclass
from itertools import combinations

from SynchDP_organize import synchDP

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


@dataclass
class Params:
    penalty: float = 0.1
    pcut: float = 0.65
    window: int = 6
    ecut: float = 3
    
    
class ReferenceAnalysis:
    def __init__(self, reference_seq: str, querry_set: str, date_info: str = None, omics_info: str = None, outdir: str = None, params: Params = Params()):
        self.reference_seq = reference_seq
        self.querry_set = querry_set
        self.date_info = date_info
        self.omics_info = omics_info
        self.outdir = outdir
        self.params = params

        
        # Loading input files
        print("[1/3] Loading files...")
        with open(self.querry_set, 'rb') as fr:
            # self.querry_dict = pickle.load(fr)
            self.querry_dict = pd.read_pickle(fr)

        # creating output dirs
        if os.path.exists(self.outdir):
            shutil.rmtree(self.outdir)  # Remove existing directory
        os.makedirs(self.outdir)        # Create new directory
            
        print("[2/3] Loading the reference sequence...")
        with open(self.reference_seq, 'rb') as fr:
            self.reference_seq_set = pickle.load(fr)
        self.Rs = self.reference_seq_set['ref_sev']
        self.Rm = self.reference_seq_set['ref_nor']
        
        
        print("[3/3] Synchronize input to the reference sequence...")
        self.results = self.synchronize_reference_seq(self.querry_dict)
        
        self.synch_plot(self.results, self.date_info, self.omics_info)
    
    def synchronize_reference_seq(self, querry_dict):
        self.total_result = {}
        for idx, (key, item) in enumerate(querry_dict.items()):
            window = int(len(item.iloc[0]) * 0.3)
            if window < 6:
                window = 6
                
            synch_sev = synchDP(item, self.Rs, self.params.penalty, self.params.pcut, window, self.params.ecut)
            for sev_results in synch_sev.results:
                sev_results.append('sev')
            synch_nor = synchDP(item, self.Rm, self.params.penalty, self.params.pcut, window, self.params.ecut)
            for nor_results in synch_nor.results:
                nor_results.append('nor')
            self.results = synch_sev.results + synch_nor.results

            filtered_result = []
            if len(self.results) == 0:
                continue
            for result in self.results:
                score1, score2, path, label = result
                diffs = []
                for x, y in path:
                    if label == 'sev':
                        diff = abs(item.iloc[0][x-1] - self.Rs.iloc[0][y-1])
                    elif label == 'nor':
                        diff = abs(item.iloc[0][x-1] - self.Rm.iloc[0][y-1])
                    diffs.append(diff)

                avg_diff = np.mean(diffs)

                if avg_diff <= self.params.ecut:
                    result_with_diff = result + [avg_diff]
                    filtered_result.append(result_with_diff)
            if len(filtered_result) == 0:
                continue
            filtered_result.sort(key=lambda x: (round(-x[0], 3), round(-x[1], 3)))
            if len(filtered_result) > 1 and filtered_result[0][-1] > filtered_result[1][-1]:
                filtered_result[0], filtered_result[1] = filtered_result[1], filtered_result[0]

            self.results = filtered_result[:10]
            self.total_result[key] = self.results
            
        return self.total_result

    
    def synch_plot(self, total_result, date_info = None, omics_info = None):
        rows, cols = 4, 3
        plots_per_page = rows * cols
        mypdf = PdfPages('%s/alignment_results_refguided.pdf'%(self.outdir))
        
        fig, axes = plt.subplots(rows, cols, figsize=(24, 16))
        axes = axes.flatten()
        plot_idx = 0
        
        for key, item in total_result.items():
            if len(item) == 0:
                continue

            result = item[0]
            score1, score2, path, label, avg_diff = result
            path = [[x-1, y-1] for x, y in path]

            y_1 = self.querry_dict[key].iloc[0].to_list()
            y_2 = self.Rs.iloc[0].to_list()
            y_3 = self.Rm.iloc[0].to_list()
            x_1 = range(len(y_1))
            x_2 = range(len(y_2))
            x_3 = range(len(y_3))

            mat_1_start, mat_1_end = path[0][0], path[-1][0]
            mat_2_start, mat_2_end = path[0][1], path[-1][1]
            mat_1_start_s, mat_1_end_s = path[0][1], mat_1_end + path[0][1] - mat_1_start
            bat = mat_2_start - mat_1_start

            ax = axes[plot_idx % plots_per_page]
            ax.set_xticks(x_2)
            ax.set_xticklabels(self.Rs.columns.to_list(), rotation=90)

            ax.plot(x_2, y_2, color='blue', label='Rs', linewidth=3, zorder = 0.7)
            ax.plot(x_3, y_3, color='green', label='Rm', linewidth=3, zorder = 0.7)


            if label == 'sev':
                ax.plot(x_2[mat_2_start:mat_2_end+1], y_2[mat_2_start:mat_2_end+1], marker='o', color='orange', label='aligned region', linewidth=3, zorder = 0.9)
            elif label == 'nor':
                ax.plot(x_3[mat_2_start:mat_2_end+1], y_3[mat_2_start:mat_2_end+1], marker='o', color='orange', label='aligned region', linewidth=3, zorder = 0.9)

            ax.plot(x_2[mat_1_start_s:mat_1_end_s+1], y_1[mat_1_start:mat_1_end+1], marker='o', color='red', label='aligned input', linewidth=3, zorder = 1)
            
            ax.plot([i + bat for i in x_1], y_1, color='black', label='input sequence', linewidth=1, zorder = 0.5)
            
            for x, y in path:
                x_point = [x + bat, y]
                y_point = [y_1[x], y_2[y] if label == 'sev' else y_3[y]]
                ax.plot(x_point, y_point, linestyle='--', linewidth=1, color='grey')

            if date_info != None and omics_info != None:
                clinical = pd.read_csv(omics_info,index_col=0).loc[self.querry_dict.keys(),:]
                omics_times = clinical.loc[:,'1time':'7time']
                
                with open(date_info,'rb') as f:
                    dates = pickle.load(f)
                data_day = {key: list(dates[key].T.columns) for key in self.querry_dict}
                
                omics_time = omics_times.loc[key,:].to_list()
                
                show_label = False
                for omics_idx, line in enumerate(omics_times.loc[key]):
                    if line in data_day[key]:
                        line_x = data_day[key].index(line)
                        if mat_1_start <= line_x <= mat_1_end:
                            if omics_idx == 0 and show_label == False:
                                ax.scatter(x=line_x + bat, y=y_1[line_x], color='dodgerblue', marker='p', label='Sampled omics', s=100, zorder=2)
                                show_label = True
                            elif omics_idx != 0 and show_label == False:
                                ax.scatter(x=line_x + bat, y=y_1[line_x], color='dodgerblue', marker='p',label='Sampled omics', s=100, zorder=2)
                                show_label = True
                            else:
                                ax.scatter(x=line_x + bat, y=y_1[line_x], color='dodgerblue', marker='p', s=100, zorder=2)
            
            ax.set_title(key, fontsize=12)
            ax.legend()
            ax.set_xlabel('Reference time')
            ax.set_ylabel('Score')

            plot_idx += 1

            if plot_idx % plots_per_page == 0:
                plt.tight_layout()
                mypdf.savefig(fig)
                plt.close(fig)
                fig, axes = plt.subplots(rows, cols, figsize=(24, 16))
                axes = axes.flatten()

        if plot_idx % plots_per_page != 0:
            for empty_ax in axes[plot_idx % plots_per_page:]:
                empty_ax.axis("off")
            plt.tight_layout()
            mypdf.savefig(fig)
            plt.close(fig)

        mypdf.close()