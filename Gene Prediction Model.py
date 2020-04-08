#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 20:23:48 2020

@author: jaoming
"""

import os
os.chdir('file directory')

import pandas as pd 
import numpy as np
from scipy.stats import norm
import seaborn as sns

# importing all the datasets
## bam holds the different reads; 'guesses' as to where the genes are
bam_data = pd.read_csv('bam_data.csv')
bam_data = bam_data.drop('Unnamed: 0', axis = 1)

## txdb holds the 'real' data
txdb_data = pd.read_csv('txdb_data.csv')
txdb_data = txdb_data.drop('Unnamed: 0', axis = 1)

# Data for trials
trial_data = bam_data[bam_data['seqnames'] == 'chr11']
trial_data = trial_data[trial_data['strand'] == '-']

class GeneLocater():
       """
       This class is meant to look into the data of a BAM file to average out
       the location of all the genes
       
       This class is proprietary to the human genome.
       """
       def __init__(self, bam_data):
              self.data = bam_data
              
       def Locate_by_Strand(self, data, threshold, variance):
              """
              Function:     Locates the potential start and end sites of a gene
                            by the strand within a chromosome
                            
              Input:        Data of only that strand of that chromosome (n x 5)
                            Threshold = for p-value
                            Variance = how close to you think the reads should be to be counted as one gene
              
              Returns:      Start and End position of all genes
                            From that strand of that chromosome
              """       
              # Making sure that the start guesses are aligned
              data = data.sort_values(['start', 'end'])
              data = data.reset_index()
              data = data.drop('index', axis = 1)
              
              # Prepping the variables to create the resultant dataframe
              seqname = data['seqnames'][0]
              strand = data['strand'][0]
              result_start = []
              result_end = []
              result_width = []
              
              temp_start = []
              temp_end = []
              
              # For Tracking
              percentiles = [int(x) for x in np.linspace(0, data.shape[0], 11, endpoint = False)]
              print('Working on the', strand + 've', 'strand of Seqname', seqname)
              print('[                              ]', '0' + '% complete')
              
              # Making the guess as to where the genes are located\
              for i in data.index:
                     if len(temp_start) <= 1:
                            temp_start.append(data.at[i, 'start'])
                            temp_end.append(data.at[i, 'end'])
                            continue
                     
                     next_start, next_end = data.at[i, 'start'], data.at[i, 'end']
                            
                     if norm.sf(next_start, loc = np.mean(temp_start), scale = variance) > threshold:
                            temp_start.append(next_start)
                            temp_end.append(next_end)

                     else:
                            result_start.append(np.mean(temp_start, dtype = int))
                            result_end.append(np.mean(temp_end, dtype = int))
                            result_width.append(int(result_end[-1] - result_start[-1]))
                                   
                            temp_start.clear()
                            temp_start.append(next_start)
                                   
                            temp_end.clear()
                            temp_end.append(next_end)
                     
                     # For Tracking
                     if i in percentiles:
                            progress = percentiles.index(i)
                            print('[' + '###'*progress + '   '*(10 - progress) + ']', str(progress*10) + '% complete')
              
              # More variables for the resultant dataframe
              result_seqnames = [seqname]*len(result_start)
              result_strand = [strand]*len(result_start)
              
              # Creating the resultant dataframe
              result = {'seqnames': result_seqnames,
                        'strand': result_strand,
                        'start': result_start,
                        'end': result_end,
                        'width': result_width}
              
              result = pd.DataFrame(result)
              return result
                
       def Locate_by_Chrom(self, data, threshold, variance):
              """
              Function:     Locates the potential start and end sites of a gene 
                            by Chromosome
                            
              Input:        Data of only that particular Chromosome (n x 5)
              
              Returns:      Start and End position of all the genes
                            From only that stated chromosome
              """
              pos_strand, neg_strand = data[data['strand'] == '+'], data[data['strand'] == '-']
              
              if pos_strand.shape[0] != 0:
                     result_pos = self.Locate_by_Strand(pos_strand, threshold, variance)
              if neg_strand.shape[0] != 0:
                     result_neg = self.Locate_by_Strand(neg_strand, threshold, variance)
              
              result = result_pos
              if neg_strand.shape[0] != 0:
                     result = result.append(result_neg)
              
              return result
       
       def Locate(self, threshold, variance = 15000):
              """
              Function:     Locates all the potential start and end sites of a gene
              
              Input:        None. Uses the BAM file from the initiation
              
              Returns:      Start and End position of all the genes
                            From all chromosomes and from either side of the strand
              """
              result = 0
              
              for chrom in self.data['seqnames'].unique():
                     print(chrom)
                     if type(result) == int:
                            result = self.Locate_by_Chrom(self.data[self.data['seqnames'] == chrom], threshold, variance)
                     else:
                            chrom_data = self.Locate_by_Chrom(self.data[self.data['seqnames'] == chrom], threshold, variance)
                            result = result.append(chrom_data)
              
              result = result.reset_index()
              result = result.drop('index', axis = 1)
              
              return result
              
gl = GeneLocater(bam_data)
gene_predict = gl.Locate(0.1)

# Evaluation Method
def gene_locater_eval(predicted, actual, threshold = 20692): # 20692 is the median width
       """
       
       Parameters
       ----------
       predicted : DataFrame (n x 5)
              Data of the predicted indexes of the genes in the genome
       actual : DataFrame (n x 5)
              Data of the actual indexes of the genes in the genome

       Returns
       -------
       A scalar metric. The lower the value the better

       """
       # Resetting all the Indexes of the Dataset for Uniformity
       predicted = predicted.reset_index(drop = True)
       actual = actual.reset_index(drop = True)
       
       # Setting up variables for the evaluation 
       strands = ['+', '-']
       chroms = actual['seqnames'].unique()
       
       predicted_obs = predicted.shape[0]
       actual_obs = actual.shape[0]
       predicted_obs_TRUE = 0 # to count how many of the predicted observations were True
       
       MSE_of_Predicted = [] # collates the error in predicted vs. actual for both start and end values
       
       best_predicted = []
       actual_predicted = []
       
       for chrom in chroms:
              print('\nStarting on seqname:', chrom)
              for strand in strands:
                     print('Starting on strand:', strand)
                     # Filtering out all the data
                     p_test = predicted[predicted['seqnames'] == chrom]
                     p_test = p_test[p_test['strand'] == strand]
                     p_test.reset_index(inplace = True, drop = True)
                     
                     a_test = actual[actual['seqnames'] == chrom]
                     a_test = a_test[a_test['strand'] == strand]
                     a_test.reset_index(inplace = True, drop = True)
                     
                     # Actual Evaluation
                     for i in a_test.index:
                            actual_start, actual_end = a_test.at[i, 'start'], a_test.at[i, 'end']
                            
                            start_thres = np.array(p_test['start'] <= actual_start + threshold) * np.array(p_test['start'] >= actual_start - threshold)
                            end_thres = np.array(p_test['end'] <= actual_end + threshold) * np.array(p_test['end'] >= actual_end - threshold)
                            if p_test[pd.Series(start_thres)].shape[0] != 0:
                                   r_guesses = p_test[pd.Series(start_thres * end_thres)]
                                   if r_guesses.shape[0] != 0:
                                          r_guesses.reset_index(inplace = True, drop = True)
                                          all_mse = [sum(x) for x in list(zip(r_guesses['start'], r_guesses['end']))]
                                          smallest_mse = min(all_mse)
                                          best_guess = r_guesses.iloc[[i for i, x in enumerate(all_mse) if x == smallest_mse][0], :].tolist()
                                          
                                          best_predicted.append(best_guess)
                                          actual_predicted.append(a_test.iloc[i, :].tolist())
                                          predicted_obs_TRUE += 1
                                          MSE_of_Predicted.append(smallest_mse**2)
       
       MSE_of_Predicted = np.round(sum(MSE_of_Predicted)**0.5, 2) # round to 2 decimal places
       print('\nMSE of Predicted Values:', MSE_of_Predicted)
       print('\nProportion of Correct Predictions over Actual:', str(np.round((predicted_obs_TRUE/actual_obs)*100, 1)) + '%')
       print('\nProportion of Correct Predictons over Predicted:', str(np.round((predicted_obs_TRUE/predicted_obs)*100, 1)) + '%')
       
       result = pd.DataFrame(best_predicted, columns = ['seqnames', 'strand', 'start', 'end', 'width'])
       result_actual = pd.DataFrame(actual_predicted, columns = ['seqnames', 'strand', 'start', 'end', 'width'])
       print('Actual Gene that was Predicted:\n', result_actual)
       print('\n\nCompleted!')
       return result

gene_locater_eval(gene_predict, txdb_data)

# Exploratory Data Analysis (EDA)
## Looking into the distribution of the widths to better estimate the threshold
w = txdb_data['width'].tolist()
w.sort()
sns.distplot(w)
sns.distplot(w[:20000])
sns.distplot(w[:10000])
sns.distplot(w[:1000])

def main():
       gl = GeneLocater(bam_data)
       gene_predict = gl.Locate(0.1)
       gene_locater_eval(gene_predict, txdb_data)
