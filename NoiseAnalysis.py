#put the columns two at a time in a dataframe
# dataframe and visualization tools
import pandas as pd
import numpy as np
import matplotlib as mlp
import time
from matplotlib import pyplot as plt

import wx
import os

import numpy.polynomial.polynomial as poly
import statistics as stats
from statistics import mode
from scipy.fft import *

import warnings
warnings.filterwarnings("ignore")

#style and formating
pd.options.display.float_format = '{:.15f}'.format
mlp.style.use('tableau-colorblind10')
mlp.rcParams['figure.dpi']= 300
mlp.rcParams['font.family'] = 'Arial'
mlp.rcParams['figure.figsize'] = [14, 10]
mlp.rcParams['figure.facecolor'] = 'white'
mlp.rcParams['axes.edgecolor'] = 'grey'
mlp.rcParams['axes.spines.top'] = False
mlp.rcParams['axes.spines.right'] = False
mlp.rcParams['axes.xmargin'] = 0.15
mlp.rcParams['axes.ymargin'] = 0.15

class NoiseAnalysis():
    def __init__(self):
        self.samples_list=[]
        self.noise_list=[]
        self.LoD_list=[]
        self.LoQ_list=[]
        self.true = ['1', 't', 'tr', 'tru', 'true', 'truee', 'y', 'ye', 'yes', 'yess', 'yeah', 'yu', 'yup', 'yupp', 'sure', 'certainly', 'yay']
        self.ChangeDefaults = 'False'
        self.SetSampleSize = 'False'
        self.SampleSize = 20000
        self.SelectRange = 'False'
        self.Start = -100000
        self.End = 100000
        self.DetectSignal ='True'
        self.Threshold = 1.0
        self.PolyFit = 'True'
        self.Degree = 4
        self.RemoveOutliers = 'False'
        self.nStandardDeviations = 0.0
        self.FourierApproximation = 'True'
        self.nHarmonics = 10
        self.RMS_noise_summary = pd.DataFrame()

    #open a windows file explorer and select file path; save file path
    def get_paths(self):
        app = wx.App(None)
        style = wx.FD_MULTIPLE
        dialog = wx.FileDialog(None, 'Select File', wildcard='*.csv;*.arw', style=style)
        if dialog.ShowModal() == wx.ID_OK:
            paths = dialog.GetPaths()
        else:
            paths = None   
        dialog.Destroy()
        return paths
    
    #read file and save data to a Dataframe
    def read_files(self, paths):
        df = pd.DataFrame()
        for path in paths:
            if path[-4:] == '.arw':
                temp_df = pd.read_csv(path, delimiter="\t", header=None)
                temp_df = temp_df[pd.to_numeric(temp_df[0], errors='coerce').notnull()].reset_index(drop=True)
                df = pd.concat([df, temp_df], axis=1)
            elif path[-4:] == '.csv':    
                temp_df = pd.read_csv(path)
                df = pd.concat([df, temp_df], axis=1)
            else:
                pass
        return df.astype('float')
    
    #interactive dialog with user
    def user_input(self):
        print(f'The program\'s default settings are:')
        print(f'''
              SelectRange: {self.SelectRange},
              DetectSignal: {self.DetectSignal}, Threshold: {self.Threshold},
              PolyFit: {self.PolyFit}, Degree: {self.Degree}, RemoveOutliers: {self.RemoveOutliers}, nStandardDeviations: {self.nStandardDeviations},
              FourierApproximation: {self.FourierApproximation}, nHarmonics: {self.nHarmonics}''')
        print('')
        self.ChangeDefaults = input('Would you like to make any changes? ')
        
        if self.ChangeDefaults.lower() in self.true:
            self.SelectRange = input(f'Would you like to enter a specific range? ') or self.SelectRange
            
            if self.SelectRange.lower() in self.true:
                self.Start = input(f'Start: ') or self.Start
                self.End = input(f'End: ') or self.End

            self.DetectSignal = input(f'Detect signals? ') or self.DetectSignal
            if self.DetectSignal.lower() in self.true:
                self.Threshold = input(f'Signal detection threshold: ') or self.Threshold

            self.PolyFit = input(f'Polynomial fit? ') or self.PolyFit
            if self.PolyFit.lower() in self.true:
                self.Degree = input(f'Polynomial fit degree: ') or self.Degree
                self.RemoveOutliers = input(f'Remove Outliers? ') or self.RemoveOutliers
                if self.RemoveOutliers.lower() in self.true:
                    self.nStandardDeviations = input(f'Number of standard deviation: ') or self.nStandardDeviations

            self.FourierApproximation = input(f'Fourier approximation? ') or self.FourierApproximation
            if self.FourierApproximation.lower() in self.true:
                self.nHarmonics = input(f'Number of harmonics to use: ') or self.nHarmonics
                
            print('')        
            print(f'Your settings are:')
            print(f'''
                  SelectRange: {self.SelectRange},
                  DetectSignal: {self.DetectSignal}, Threshold: {self.Threshold}, 
                  PolyFit: {self.PolyFit}, Degree: {self.Degree}, RemoveOutliers: {self.RemoveOutliers}, nStandardDeviations: {self.nStandardDeviations},
                  FourierApproximation: {self.FourierApproximation}, nHarmonics: {self.nHarmonics}''')
            print('')
        return None
    
    #option to control sample size
    def set_sample_size(self, x, y, sample_size):
        x_new = np.linspace(min(x), max(x), sample_size)
        # Where you want to interpolate    
        y_new = np.interp(x_new, x, y) 
        return x_new, y_new
    
    #option to select a specific range to operate on
    def select_range(self, x, y, Start, End):
        keep = np.zeros(len(x))
        for i in range(len(x)):
            if x[i] > Start and x[i] < End:
                keep[i] = 1
        return x[keep==1], y[keep==1]

    #classify each data point as either baseline (0) or signal (1) 
    def signal_baseline_classifier(self, y, signal_threshold, lag_fraction=0.03, draw_baseline=True):
        #use a SMA as a lagging baseline to determine signal
        lag = int(len(y)*lag_fraction)
        len_data = len(y)
        threshold = signal_threshold*y.std() #min(y.std()/10 , y[:lag].std())
        signal = np.zeros(len_data)

        for i in range(lag, len_data):
            SMA_i = np.sum(y[i-lag:i+lag])/len(y[i-lag:i+lag])

            if abs(y[i]-SMA_i) >= threshold:
                signal[i] = 1

        #correct any false negatives points by conforming to nearest n neighboors
        n_neighbors = max(1, int(lag/5))
        s = signal.copy()
        for i in range(n_neighbors, len_data):
            if s[i] == 0 and mode(s[i-n_neighbors:i+n_neighbors]) == 1:
                signal[i-n_neighbors:i+n_neighbors] = mode(s[i-n_neighbors:i+n_neighbors])

        #characterize baseline points around signals as signals to reduce false negatives
        s = signal.copy()
        for i in range(n_neighbors,len_data):
            if s[i] == 1:
                signal[i-3*n_neighbors:i+3*n_neighbors] = 1

        #recreate baseline as a copy of y without signals  
        if draw_baseline:          
            baseline = pd.Series(y.copy())
            baseline[signal==1] = np.nan
            baseline = baseline.interpolate()
            for i in range(len_data):
                baseline[i] = min(y[i], baseline[i])
        else: 
            baseline = 'N/a'
        return signal
        
    #creat a tunnel-like polynomial fit of the baseline; this is can be used to flatten and/or remove outliers
    def polynomial_tunnel_fit(self, x, y, deg, n_std, n_iters=1, remove_outliers=True, flatten=True):
        #Runge's phenomenon is a known issue with this method
        i = 0
        outlier = np.zeros(len(y))
        while i < n_iters:  
            coefs = poly.polyfit(x, y, deg)
            ffit = poly.polyval(x, coefs)

            top = ffit + n_std*y.std()
            base = ffit - n_std*y.std()
            
            if remove_outliers:
                toutlier = y > top
                boutlier = y < base
                outlier = toutlier | boutlier
                x = x[~outlier]
                y = y[~outlier]
                
                top = top[~outlier] 
                base = base[~outlier]
                
            if flatten:
                y = y-base 
                y = y-y.mean()
                
            if i == 0:
                int_top = top
                int_base = base
            
            i += 1     
        return x, y, int_top, int_base
    
    def fourier_transformation_approximation(self, y, nHarmonics):
    
        # Number of samples in our input signal
        N = len(y)
        #This function computes a fourier series representation of our input signal; 
        #using the 1-D discrete Fourier Transform (DFT) of a real-valued array by means of an efficient algorithm called the Fast Fourier Transform (FFT).
        fourier_series = rfft(y)
        #reconstruct signal from the inverse of the real valued components of the fourier series
        #only use the first 'n' number pf preodic components from the fourier series to reconstruct the signal
        y_approx = irfft(fourier_series[:nHarmonics], N)
        return y_approx
    
    #produce short report 
    def short_report_grapher(self, x2, y2, LoB, LoD, LoQ):

        #plot cleaned baselines + LOD/LOQ thresholds
        fig, ax = plt.subplots(figsize=(12, 6))
        plt.gca().ticklabel_format(axis='both', style='plain', useOffset=False)
        plt.suptitle(f'{df.columns[1]}', fontsize=12, y = 0.94, fontweight='bold')

        ax.scatter(x2, y2, s=0.5)
        LoB_=ax.hlines(LoB, xmin=min(x), xmax=max(x), linestyle= '-',  alpha=0.4, linewidth=0.6)
        LoD_=ax.hlines(LoB+LoD, xmin=min(x2), xmax=max(x2), linestyle= ':', alpha=0.9, linewidth=1.2)
        LoQ_=ax.hlines(LoB+LoQ, xmin=min(x2), xmax=max(x2), linestyle= '-',  alpha=0.9, linewidth=1.2)

        ax.set_xlabel(f'{df.columns[0]}', fontsize=11)
        ax.set_ylabel('signal', fontsize=11)
        ax.legend([LoQ_, LoD_], ['LoQ', 'LoD'], frameon=False, bbox_to_anchor=(1.05, 1), loc='upper right', handlelength=0.5)
        return fig.savefig(f'{os.path.dirname(path)}\\{fig._suptitle.get_text()}_rms_noise.png', facecolor=fig.get_facecolor(), dpi=fig.dpi)
    
    
na = NoiseAnalysis()

paths = na.get_paths()
input_data = na.read_files(paths=paths)

na.user_input()

print('')
print(f'working...')

numcols = len(input_data.columns)
for i in range(numcols):
    #i = 0,2,4,6,8...etc.
    if i%2 == 0:
        #generate temporary dataframe and define temp x,y 
        df = pd.DataFrame(input_data.iloc[:,i:i+2]).dropna()
        N = int(len(df)*0.015)
        x = df.iloc[N:-N,0].values
        y = df.iloc[N:-N,1].values
        
        if na.SetSampleSize.lower() in na.true:
            x, y = na.set_sample_size(x, y, sample_size= na.SampleSize)
            
        x2 = x.copy()
        y2 = y.copy()
        signal=np.zeros(len(x2))
        
        fig, axs = plt.subplots(2, 2)
        plt.gca().ticklabel_format(axis='both', style='plain', useOffset=False)
        plt.suptitle(f'{df.columns[1]}', fontsize=12, y = 0.94, fontweight='bold')
        
        if na.SelectRange.lower() in na.true:
            x, y = na.select_range(x2, y2, int(na.Start), int(na.End))
            x2 = x.copy()
            y2 = y.copy()
            signal=np.zeros(len(x2))
        
        if na.DetectSignal.lower() in na.true:
            signal = na.signal_baseline_classifier(y=y, signal_threshold= float(na.Threshold), lag_fraction=0.03, draw_baseline=True)
            
            b=axs[0, 0].scatter(x[signal==0], y[signal==0], s=0.2)
            s=axs[0, 0].scatter(x[signal==1], y[signal==1], s=0.2)
            axs[0, 0].set_title(f'Signal Detection (threshold={na.Threshold})', fontsize=10)
            axs[0, 0].set_xlabel(f'{df.columns[0]}', fontsize=9)
            axs[0, 0].set_ylabel(f'signal', fontsize=9)
            axs[0, 0].legend([b, s],[ 'Baseline', 'Signal'], fontsize=8, frameon=False, bbox_to_anchor=(1.05, 1), loc='upper right', handlelength=1)
            
            x2 = x2[signal==0]
            y2 = y2[signal==0]
            
        if na.PolyFit.lower() in na.true:
            x2, y2, topline, bottomline = na.polynomial_tunnel_fit(x2, y2, deg=int(na.Degree), n_std=float(na.nStandardDeviations), n_iters=1,
                                                                 remove_outliers= na.RemoveOutliers.lower() in na.true, flatten=True)
        
            b=axs[0, 1].scatter(x[signal==0], y[signal==0], s=0.2)
            t=axs[0, 1].scatter(x2, topline, s=0.2, color='#ff7f0e', alpha=0.4, linewidth=0.6)
            b2=axs[0, 1].scatter(x2, bottomline, s=0.2, color='#ff7f0e', alpha=0.4, linewidth=0.6)
            axs[0, 1].set_title(f'Polynomial Fit (degree={na.Degree})', fontsize=10)
            axs[0, 1].set_xlabel(f'{df.columns[0]}', fontsize=9)
            axs[0, 1].set_ylabel(f'signal', fontsize=9)
            axs[0, 1].legend([b, t],['Baseline', 'Polynomial Fit'], fontsize=8, frameon=False, bbox_to_anchor=(1.05, 1), loc='upper right', handlelength=1)
        
        if na.FourierApproximation.lower() in na.true:
            y_approx = na.fourier_transformation_approximation(y=y2, nHarmonics=int(na.nHarmonics))
            
            b=axs[1, 0].scatter(x2, y2, s=0.2)
            a=axs[1, 0].scatter(x2, y_approx, s=0.2, color='#ff7f0e', alpha=0.4, linewidth=0.4)
            axs[1, 0].set_title(f'Fourier Approximation Using The First {na.nHarmonics} Harmonics', fontsize=10)
            axs[1, 0].set_xlabel(f'{df.columns[0]}', fontsize=9)
            axs[1, 0].set_ylabel(f'signal', fontsize=9)
            axs[1, 0].legend([b, a],[ 'Baseline', 'Fourier Approximation'], fontsize=8, frameon=False, bbox_to_anchor=(1.05, 1), loc='upper right', handlelength=1)
            
            y2=y2-y_approx

        #calculate LOD/LOQ
        y2 = y2 - y2.mean()
        noise = y2.std()
        LoD = 3*noise
        LoQ = 10*noise
        
        #graph 4th quadrant with final baseline and LoD/LoQ horizontal lines
        axs[1, 1].scatter(x2, y2, s=0.2)
        axs[1, 1].hlines(0, xmin=min(x2), xmax=max(x2), linestyle= '-', color='#000000', alpha=0.4, linewidth=0.6)
        LoD_line=axs[1, 1].hlines(LoD, xmin=min(x2), xmax=max(x2), linestyle= ':', color='#000000', alpha=0.8, linewidth=1.0)
        LoQ_line=axs[1, 1].hlines(LoQ, xmin=min(x2), xmax=max(x2), linestyle= '-', color='#000000', alpha=0.8, linewidth=1.0)
        axs[1, 1].set_title('Baseline Noise Evaluation', fontsize=10)
        axs[1, 1].set_xlabel(f'{df.columns[0]}', fontsize=9)
        axs[1, 1].set_ylabel('signal', fontsize=9)
        axs[1, 1].legend([LoQ_line, LoD_line], ['LoQ', 'LoD'], fontsize=8, frameon=False, bbox_to_anchor=(1.05, 1), loc='upper right', handlelength=0.5)

        fig.savefig(f'{os.path.dirname(paths[0])}\\_{fig._suptitle.get_text()}_rms_noise.png', facecolor=fig.get_facecolor(), dpi=fig.dpi)
    
        #collect LOD/LOQ data in lists
        na.samples_list.append(f'{df.columns[1]}')
        na.noise_list.append(noise)
        na.LoD_list.append(LoD)
        na.LoQ_list.append(LoQ)
        
        print('')
        print(f'({int((i/2)+1)}/{int(numcols/2)})')
                
#build summary dataframe
idx = int(numcols/2)
if idx>1:
    #Final average
    na.noise_list.append(stats.mean(na.noise_list[:idx]))
    na.LoD_list.append(stats.mean(na.LoD_list[:idx]))
    na.LoQ_list.append(stats.mean(na.LoQ_list[:idx]))
    na.samples_list.append(f'Average')

    #Final standard deviation`
    na.noise_list.append(stats.stdev(na.noise_list[:idx]))
    na.LoD_list.append(stats.stdev(na.LoD_list[:idx]))
    na.LoQ_list.append(stats.stdev(na.LoQ_list[:idx]))
    na.samples_list.append(f'Standard Deviation')

# [Samples,LoD,LoQ]
na.RMS_noise_summary = pd.DataFrame({'Summary': na.samples_list, 'noise (Standard deviation of baseline)': na.noise_list, 'LoD (3*noise)': na.LoD_list, 'LoQ (10*noise)': na.LoQ_list})

#save summary dataframe to .csv
if len(paths)>1:
    na.RMS_noise_summary.to_csv(f'{os.path.dirname(paths[0])}\\_{len(paths)}_rms_noise.csv', index = False)
else:
    na.RMS_noise_summary.to_csv(f'{os.path.dirname(paths[0])}\\_{os.path.splitext(os.path.basename(paths[0]))[0]}_rms_noise.csv', index = False)
print('')
print('Your analysis is complete.')
print('')
time.sleep(5)