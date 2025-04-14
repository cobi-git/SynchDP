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

import networkx as nx
from community import community_louvain
import matplotlib
matplotlib.use('Agg')  # Before importing pyplot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



@dataclass
class Params:
    penalty: float = 0.01
    pcut: float = 0.7
    window: int = 6
    ecut: float = 3.5

    
class DenovoAnalysis:
    def __init__(self, querry_set: str, cover_ratio: int = 95, date_info: str = None, omics_info: str = None, outdir: str = None, params: Params = Params()):
        self.querry_set = querry_set
        self.cover_ratio = cover_ratio
        self.date_info = date_info
        self.omics_info = omics_info
        self.outdir = outdir
        self.params = params

        
        print("[1/5] Loading files..")
        # Load files
        with open(self.querry_set, 'rb') as fr:
            # self.querry_dict = pickle.load(fr)
            self.querry_dict = pd.read_pickle(fr)

        # creating output dirs
        if os.path.exists(self.outdir):
            shutil.rmtree(self.outdir)  # Remove existing directory
        os.makedirs(self.outdir)        # Create new directory

        
        print("[2/5] Performing all pairwise alignment...")
        # Generate pairwise matrix using SynchDP
        self.pair_matrix = self.calculate_matrix(self.querry_dict)        
        
        
        print("[3/5] Constructing alignment network...")
        # Constructing the alignment network
        self.H = self.construct_network()
        
        # computing network centrality
        self.centrality_df_ = self.calculate_centrality(self.H)
        self.degree_centrality = nx.degree_centrality(self.H)
        
        
        print("[4/5] Generating the reference sequence de novo...")
        # Generate hub sequence
        self.hub_nodes = self.calculate_hub_node(self.H)
        self.reference_sequence_ = self.generate_reference_sequence()
        self.network_plot()
        
        print("[5/5] Synchronize intput to reference sequence...")        
        self.results = self.synchronize_reference_seq(self.querry_dict, self.reference_sequence_)
        
        self.synch_plot(self.results, self.date_info, self.omics_info)
        
        
    def calculate_matrix(self, querry_dict):
        key_list = list(querry_dict.keys())
        result_df = pd.DataFrame(index=key_list, columns=key_list)

        for key1, key2 in combinations(key_list, 2):
            list1 = querry_dict[key1]
            list2 = querry_dict[key2]
            window = max(int(min(len(list1.iloc[0]), len(list2.iloc[0])) * 0.3), 5)
            synch = synchDP(list1, list2, self.params.penalty, self.params.pcut, window, self.params.ecut)
            if len(synch.results):
                result_df.at[key1, key2] = synch.results[0][1]
                result_df.at[key2, key1] = synch.results[0][1]

        result_df=result_df.fillna(0).infer_objects(copy=False)
        return result_df

    
    def construct_network(self):
        G = nx.Graph()
        for i, source in enumerate(self.pair_matrix.index):
            for j, target in enumerate(self.pair_matrix.columns):
                if i < j:
                    weight = self.pair_matrix.loc[source, target]
                    if weight > 0:
                        G.add_edge(source, target, weight=weight)
                        
        # Louvain clustering               
        self.partition = community_louvain.best_partition(G, weight='weight')
        
        # Threshold filtering (top 40%)
        weights = [d['weight'] for _, _, d in G.edges(data=True)]
        threshold = pd.Series(weights).quantile(0.4)
        
        # Limited graph
        H = nx.Graph()
        for u, v, d in G.edges(data=True):
            if d['weight'] >= threshold:
                H.add_edge(u, v, weight=d['weight'])
                
        # Node layout by cluster
        self.pos = {}
        clusters = {}
        for node, cluster_id in self.partition.items():
            clusters.setdefault(cluster_id, []).append(node)

        cluster_offset = 0
        for cluster_id, nodes in clusters.items():
            subgraph = H.subgraph(nodes)
            cluster_pos = nx.spring_layout(subgraph, seed=42)
            for node, coords in cluster_pos.items():
                self.pos[node] = [coords[0] + cluster_offset, coords[1]]
            cluster_offset += 2.2
            
        return H
            
        
    def calculate_centrality(self, H):
        degree = nx.degree_centrality(H)
        betweenness = nx.betweenness_centrality(H, weight='weight')
        closeness = nx.closeness_centrality(H)
        eigenvector = nx.eigenvector_centrality(H, max_iter=1000, tol=1e-6)

        df = pd.DataFrame({
            'Node': list(H.nodes),
            'Degree': [degree[node] for node in H.nodes],
            'Betweenness': [betweenness[node] for node in H.nodes],
            'Closeness': [closeness[node] for node in H.nodes],
            'Eigenvector': [eigenvector[node] for node in H.nodes],
            'Degree_Neighbors': [len(list(H.neighbors(node))) for node in H.nodes]
        }).sort_values(by='Degree', ascending=False).reset_index(drop=True)

        return df

    
    def calculate_hub_node(self, H):
        degree = nx.degree_centrality(H)
        sorted_nodes = sorted(degree.items(), key=lambda x: x[1], reverse=True)

        self.covered = set()
        required = []
        all_nodes = set(H.nodes)

        for node, _ in sorted_nodes:
            required.append(node)
            self.covered.update([node] + list(H.neighbors(node)))

            coverage = (len(self.covered) / len(all_nodes)) * 100
            if coverage >= self.cover_ratio:
                break

        return required
    
    
    def generate_reference_sequence(self):
        top_nodes = self.hub_nodes
        max_pairs = len(self.hub_nodes) - 1
        self.filtered_edges = pd.DataFrame([
            {"Node1": u, "Node2": v, "Weight": d["weight"]}
            for u, v, d in self.H.edges(data=True)
            if u in top_nodes and v in top_nodes
        ])
        self.filtered_edges = self.filtered_edges.sort_values(by="Weight", ascending=False).reset_index(drop=True)

        sequence = []
        used_nodes = set()
        self.results_dict = {}
        self.shift_dict = {}

        while len(sequence) < max_pairs:
            for _, row in self.filtered_edges.iterrows():
                u, v = row['Node1'], row['Node2']

                if not sequence:
                    sequence.append((u, v))
                    used_nodes.update([u, v])
                    break
                elif (u in used_nodes or v in used_nodes) and not (u in used_nodes and v in used_nodes):
                    sequence.append((u, v))
                    used_nodes.update([u, v])
                    break
            else:
                break

        centrality_ranking = {node: i+1 for i, node in enumerate(top_nodes)}
        for node1, node2 in sequence:
            df_1 = self.querry_dict[node1]
            df_2 = self.querry_dict[node2]
            synch = synchDP(df_1, df_2, self.params.penalty, self.params.pcut, self.params.window, self.params.ecut)
            coordinates = synch.results[0][2]

            x_list = [x-1 for x, y in coordinates]
            y_list = [y-1 for x, y in coordinates]

            mean_value_list = []
            for x, y in coordinates:
                adjusted_x = x-1
                adjusted_y = y-1
                if adjusted_x not in df_1.columns or adjusted_y not in df_2.columns:
                    continue
                mean_value = (df_1.iloc[0, adjusted_x] + df_2.iloc[0, adjusted_y]) / 2
                mean_value_list.append(mean_value)

            key = f'{node1} and {node2}'

            if centrality_ranking[node1] <= centrality_ranking[node2]:
                selected_time = x_list
                shift_val = min(x_list) - min(y_list)
                self.shift_dict[key] = {
                    node1: 0,
                    node2: shift_val
                }
            else:
                selected_time = y_list
                shift_val = min(y_list) - min(x_list)
                self.shift_dict[key] = {
                    node1: shift_val,
                    node2: 0
                }

            save_df = pd.DataFrame(columns=selected_time)
            save_df.loc[len(save_df)] = mean_value_list
            self.results_dict[key] = save_df

        
        reference_sequence = None
        centrality_rank = {node: i for i, node in enumerate(self.hub_nodes)}
        shift_map = {}

        for idx, key in enumerate(self.results_dict.keys()):
            df = self.results_dict[key].copy()
            nodes = list(self.shift_dict[key].keys())

            if idx == 0:
                reference_sequence = df
                for sample, shift_num in self.shift_dict[key].items():
                    shift_map[sample] = shift_num
                continue

            if nodes[0] in shift_map.keys() and nodes[1] not in shift_map.keys():
                base_node = nodes[0]
            elif nodes[1] in shift_map.keys() and nodes[0] not in shift_map.keys():
                base_node = nodes[1]
            elif nodes[0] in shift_map.keys() and nodes[1] in shift_map.keys():
                base_node = nodes[0] if centrality_rank[nodes[0]] < centrality_rank[nodes[1]] else nodes[1]
            else:
                continue

            shift_diff = shift_map[base_node] - self.shift_dict[key][base_node]

            df.columns = [c + shift_diff for c in df.columns]
            reference_sequence = pd.concat([reference_sequence, df], axis=1)
            # reference_sequence = reference_sequence.groupby(reference_sequence.columns, axis=1).mean()
            reference_sequence = reference_sequence.T.groupby(reference_sequence.columns).mean().T

            for sample, shift_num in self.shift_dict[key].items():
                if sample not in shift_map:
                    shift_map[sample] = shift_num + shift_diff

        reference_sequence.columns = [c - reference_sequence.columns[0] for c in reference_sequence.columns]
        reference_sequence.index = ['ref_seq']

        # plot reference sequence
        self.reference_seq_plot(reference_sequence)
        
        return reference_sequence
    
    
    def synchronize_reference_seq(self, querry_dict, reference_sequence):
        self.total_results = {}
        for idx, (key, item) in enumerate(querry_dict.items()):
            if key not in self.covered:
                continue
                
            window = int(len(item.iloc[0]) * 0.3)
            if window < 6:
                window = 6
                
            synch = synchDP(item, self.reference_sequence_, self.params.penalty, self.params.pcut, window, self.params.ecut)
            
            if len(synch.results) == 0:
                continue
                
            #synch.results.sort(key=lambda x: (round(-x[0], 3), round(-x[1], 3)))
            self.total_results[key] = synch.results[:10]
        
        return self.total_results
    
    
    def network_plot(self):
        unique_clusters = list(set(self.partition.values()))
        color_map = plt.cm.get_cmap('Set3', len(unique_clusters))
        cluster_colors = {c: color_map(i) for i, c in enumerate(unique_clusters)}
        custom_colors = {0: 'white', 1: 'lightgray'}

        node_colors = [custom_colors.get(self.partition[n], cluster_colors[self.partition[n]]) for n in self.H.nodes()]
        node_border_colors = ['red' if n in self.hub_nodes else 'black' for n in self.H.nodes()]
        node_border_widths = [2.5 if n in self.hub_nodes else 0.0001 for n in self.H.nodes()]
        node_sizes = [self.degree_centrality[n] * 2000 for n in self.H.nodes()]
        
        # 엣지 분리
        hub_edges = [(u, v) for u, v in self.H.edges() if u in self.hub_nodes and v in self.hub_nodes]
        normal_edges = [(u, v) for u, v in self.H.edges() if (u, v) not in hub_edges]

        # 시각화 시작
        plt.figure(figsize=(15, 15))

        normal_weights = [self.H[u][v]['weight'] * 0.0001 for u, v in normal_edges]
        normal_colors = [plt.cm.Greys(self.H[u][v]['weight']) for u, v in normal_edges]
        normal_lc = nx.draw_networkx_edges(self.H, self.pos,edgelist=normal_edges,edge_color=normal_colors,width=normal_weights,style='solid',alpha=0.3)
        if normal_lc:
            normal_lc.set_zorder(1)

        hub_lc = nx.draw_networkx_edges(self.H, self.pos,edgelist=hub_edges,edge_color='red',width=2.0,style='solid',alpha=1.0)
        if hub_lc:
            hub_lc.set_zorder(2)

        nx.draw_networkx_nodes(self.H, self.pos,node_color=node_colors,node_size=node_sizes,
                               edgecolors=node_border_colors,linewidths=node_border_widths)

        nx.draw_networkx_labels(self.H, self.pos, font_size=8)

        plt.axis("off")
        plt.savefig('%s/alignment_network.pdf'%(self.outdir), dpi=300)
        # plt.show()
        

    def synch_plot(self, total_result, date_info = None, omics_info = None):
        rows, cols = 4, 3
        plots_per_page = rows * cols
        mypdf = PdfPages('%s/alignment_results_denovo.pdf'%(self.outdir))
        
        fig, axes = plt.subplots(rows, cols, figsize=(24, 16))
        axes = axes.flatten()
        plot_idx = 0
        
        for key, item in total_result.items():
            if len(item) == 0:
                continue

            result = item[0]
            score1, score2, path = result
            path = [[x-1, y-1] for x, y in path]

            y_1 = self.querry_dict[key].iloc[0].to_list()
            y_2 = self.reference_sequence_.iloc[0].to_list()
            
            x_1 = range(len(y_1))
            x_2 = range(len(y_2))


            mat_1_start, mat_1_end = path[0][0], path[-1][0]
            mat_2_start, mat_2_end = path[0][1], path[-1][1]
            mat_1_start_s, mat_1_end_s = path[0][1], mat_1_end + path[0][1] - mat_1_start
            bat = mat_2_start - mat_1_start

            ax = axes[plot_idx % plots_per_page]
            ax.set_xticks(x_2)
            ax.set_xticklabels(self.reference_sequence_.columns.to_list(), rotation=90)

            ax.plot(x_2, y_2, color='blue', label='Reference Sequence', linewidth=3, zorder = 0.7)

            ax.plot(x_2[mat_2_start:mat_2_end+1], y_2[mat_2_start:mat_2_end+1], marker='o', color='orange', label='aligned region', linewidth=3, zorder = 0.9)
            

            ax.plot(x_2[mat_1_start_s:mat_1_end_s+1], y_1[mat_1_start:mat_1_end+1], marker='o', color='red', label='aligned input', linewidth=3, zorder = 1)
            
            ax.plot([i + bat for i in x_1], y_1, color='black', label='intput sequence', linewidth=1, zorder = 0.5)
            
            for x, y in path:
                x_point = [x + bat, y]
                y_point = [y_1[x], y_2[y]]
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
        
        
        
    def reference_seq_plot(self, reference_seq):
        plt.figure(figsize=(15, 6))
        plt.plot(reference_seq.columns, reference_seq.iloc[0], label='Reference sequence', linestyle='-', marker='o', color = 'red', alpha = 0.7)
        plt.xlabel('sequence time')
        plt.ylabel('score')
        plt.legend()

        plt.savefig('%s/denovo_reference_sequence.pdf'%(self.outdir))

        
    def matrix(self):
        return self.pair_matrix

    def graph(self):
        return self.H
    
    def centrality_df(self):
        return self.centrality_df_

    def hub_node(self):
        return self.hub_nodes
    
    def reference_sequence(self):
        return self.reference_sequence_