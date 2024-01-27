import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from enum import IntEnum
import glob
import seaborn as sns

class Trace(IntEnum):
    STOP = 0
    LEFT = 1
    UP = 2
    DIAGONAL = 3

def max_xy(ndarr):
    argamx = np.argmax(ndarr)
    return (argamx // ndarr.shape[1], argamx % ndarr.shape[1])

def plot_data(data1,data2):
    fig,axs = plt.subplots(nrows=2,ncols=1,figsize=(20,10))
    y_1,y_2 = data1.iloc[0].to_list(),data2.iloc[0].to_list()
    x_1,x_2 = range(1,len(y_1)+1), range(1,len(y_2)+1)
    axs[0].plot(x_1,y_1,marker='o',color='blue')
    axs[0].set_xticks(x_1,data1.columns.to_list())
    axs[0].set_title(data1.name)
    
    axs[1].plot(x_2,y_2,marker='o',color='blue')
    axs[1].set_xticks(x_2, data2.columns.to_list())
    axs[1].set_title(data2.name)
    
    
def plot_data_warp(data1,data2,path,window):
    if(len(data1.iloc[0].to_list()) > len(data2.iloc[0].to_list())):
        data1, data2 = data2, data1
    fig,axs = plt.subplots(nrows=2,ncols=1,figsize=(20,10))
    y_1,y_2 = data1.iloc[0].to_list(),data2.iloc[0].to_list()
    x_1,x_2 = range(1,len(y_1)+1), range(1,len(y_2)+1)
    x_1_start,x_1_end = path[0][0], path[-1][0]
    x_2_start,x_2_end = path[0][1], path[-1][1]
    
    
    axs[0].plot(x_1,y_1,marker='o',color='blue')
    axs[0].plot(x_1[x_1_start:x_1_end], y_1[x_1_start:x_1_end],marker='o', color='red')
    axs[0].set_xticks(x_1,data1.columns.to_list())
    #axs[0].set_title(data1.name)
    
    axs[1].plot(x_2,y_2,marker='o',color='blue')
    axs[1].plot(x_2[x_2_start:x_2_end], y_2[x_2_start:x_2_end],marker='o', color='red')
    axs[1].set_xticks(x_2, data2.columns.to_list())
    #axs[1].set_title(data2.name)
    
    
def plot_data_warp_with_gap(data1,data2,path,backtrace_matrix,window):
    if(len(data1.iloc[0].to_list()) > len(data2.iloc[0].to_list())):
        data1,data2 = data2, data1
        
    y_1,y_2 = data1.iloc[0].to_list(),data2.iloc[0].to_list()
    x_1,x_2 = range(1,len(y_1)+1), range(1,len(y_2)+1)
    x_1_start,x_1_end = path[0][0] - window, path[-1][0]
    x_2_start,x_2_end = path[0][1] - window, path[-1][1]
    
    x1_gap, x2_gap = [],[]
    #Backtrace
    for (i,j) in path:
        if backtrace_matrix[i,j] == Trace.UP:
            x1_gap.append(i)
        
        elif backtrace_matrix[i,j] == Trace.LEFT:
            x2_gap.append(j)
            
    y1_gap = [y_1[i] for i in x1_gap]
    y2_gap = [y_2[i] for i in x2_gap]

            
    fig,axs = plt.subplots(nrows=2, ncols=1, figsize=(20,10))
    axs[0].plot(x_1,y_1,marker='o',color='blue')
    axs[0].plot(x_1[x_1_start:x_1_end], y_1[x_1_start:x_1_end],marker='o', color='red')
    axs[0].plot(x1_gap, y1_gap, marker = 'o', mfc='g')
    axs[0].set_xticks(x_1,data1.columns.to_list())
    axs[0].set_title(data1.name)
    
    axs[1].plot(x_2,y_2,marker='o',color='blue')
    axs[1].plot(x_2[x_2_start:x_2_end], y_2[x_2_start:x_2_end],marker='o', color='red')
    axs[1].plot(x2_gap,y2_gap,marker = 'o', mfc = 'g')
    axs[1].set_xticks(x_2, data2.columns.to_list())
    axs[1].set_title(data2.name)           
    
    
class NEWS_data:
    def __init__(self,data_dir):
        datas = glob.glob(f'{data_dir}/*.csv')
        datas.sort()
        df = dict()
        days = dict()
        for data in datas:
            cols = pd.read_csv(data,nrows=1).columns
            _ = pd.read_csv(data,usecols=cols)
            _.dropna(inplace=True)
            if data[52:63] == 'COV-CCO-342':
                _.at[13,'DATE'] = '2021-09-12'
            
            days[data[52:63]] = _['DATE'].to_list()
            _['DATE'] = pd.to_datetime(_['DATE'])
            first_day = _['DATE'].iloc[0]
            _['DATE_DIFF'] = (_['DATE'] - first_day).dt.days
            _ = _.loc[:,['NEWS', 'DATE_DIFF']]
            # _ = _.loc[:,['WHO_SCORE', 'DATE_DIFF']]
            _.set_index('DATE_DIFF', inplace=True)
        
            _ = _.T
            _.name = data[52:63]
            df[data[52:63]] = _
        
        
        if 'COV-CCO-240' in df.keys():
            df['COV-CCO-240'].rename(columns={294:1}, inplace=True)   
        self.df = df
        self.patients = list(self.df.keys())
        self.days = days 
        
def plot_align_omics(data1,data2, window,path,news_time,omics_time):
    fig,axs = plt.subplots(nrows=2,ncols=1,figsize=(20,10))
    plt.subplots_adjust(hspace=0.5)

    y_1,y_2 = data1.iloc[0].to_list(),data2.iloc[0].to_list()
    x_1,x_2 = range(len(y_1)), range(len(y_2))
    x_1_start,x_1_end = path[0][0] - window, path[-1][0]
    x_2_start,x_2_end = path[0][1] - window, path[-1][1]


    axs[0].plot(x_1,y_1,marker='o',color='blue')
    axs[0].plot(x_1[x_1_start:x_1_end], y_1[x_1_start:x_1_end],marker='o', color='red', label = 'matching location')
    axs[0].set_xticks(x_1,news_time,rotation=45)
    axs[0].set_title(data1.name)

    axs[1].plot(x_2,y_2,marker='o',color='blue')
    axs[1].plot(x_2[x_2_start:x_2_end], y_2[x_2_start:x_2_end],marker='o', color='red')
    axs[1].set_xticks(x_2, data2.columns.to_list())
    axs[1].set_title(data2.name)

    omic_point = list()
    not_warp_point = list()
    for point in news_time[x_1_start:x_1_end]:
        if point in omics_time: omic_point.append(point)
        
    for point in news_time:
        if point in omics_time and point not in omic_point:
            not_warp_point.append(point)
    
    
    cnt = 0
    points = omic_point + not_warp_point
    sorted_points = sorted(points)

    for point in sorted_points:
        if point in omic_point:
            axs[0].axvline(x=news_time.index(point), color='g', linestyle='--', zorder=0.5, label=f'{cnt} omics time')
        else:
            axs[0].axvline(x=news_time.index(point), color='orange', linestyle='--', zorder=0.5, label=f'{cnt} omics time')
        cnt += 1
    axs[0].legend(loc=[1.01,0])
    
class NKNK_data:
    def __init__(self,dir):
        datas = glob.glob(f'{dir}/*.csv')
        datas.sort()
        
        df = dict()
        for csv in datas:
            _ = pd.read_csv(csv,index_col=0)
            days = [int(day[1:]) for day in list(_.columns)]
            for idx in range(len(_)):
                name = _.iloc[idx].name
                df[name] = pd.DataFrame([_.iloc[idx].to_list()],columns = days, index = [name])
                df[name].name = name
                
        self.df = df

    