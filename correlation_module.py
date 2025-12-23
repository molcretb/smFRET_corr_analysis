# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 09:55:08 2025

@author: molcre0000
"""

import numpy as np




def Correlation_1D(X1_traj, X0_traj):
    """
    Calculates a linear correlation function between two trajectories of equal length.

    Parameters
    ----------
    X1_traj : Numpy array
        Delayed trajectory
    X0_traj : Numoy arrray
        First trajectory

    Returns
    -------
    Correlation : Numpy array
        1D array of correlation values as a function of delays in Tau
    Tau :  Numpy array
        1D array of delays corresponding to Correlation
    DTau : Numpy array
        1D array of time-bin widths corresponding to Tau
    Weight : Numpy array
        1D array of the number of original time points within each time bin

    """
    MinLength = 8 # Below this number of time points, the calculation stops
    CurrentTimeSpacing = 1 # Current time resolution
    CurrentLength = len(X0_traj) # Current length trajectory
    
    Correlation = [] # np.zeros(CurrentLength)
    DTau = [] # np.zeros(CurrentLength)
    Tau = [] # np.zeros(CurrentLength)
    Weight = [] # np.zeros(CurrentLength)
    
    X1_traj = np.array(X1_traj)
    X0_traj = np.array(X0_traj)
    
    for k in range(4):
        
        # Calculate the correlation function
        Correlation.append(np.sum(X1_traj[k:CurrentLength-1] * X0_traj[0:CurrentLength-k-1]) / (CurrentLength - k))
        
        # Calculate the other outputs
        DTau.append(CurrentTimeSpacing)
        Tau.append(k * CurrentTimeSpacing)
        Weight.append((CurrentLength-k) * CurrentTimeSpacing)
        
    
    while CurrentLength >= MinLength:
        # Calculate four delay points at the current resolution
        for k in range(4,8):
            
            # Calculate the correlation function
            Correlation.append(np.sum(X1_traj[k:CurrentLength-1] * X0_traj[0:CurrentLength-k-1]) / (CurrentLength - k))
            
            # Calculate the other outputs
            DTau.append(CurrentTimeSpacing)
            Tau.append(k * CurrentTimeSpacing)
            Weight.append((CurrentLength-k) * CurrentTimeSpacing)
        
        # Bin trajectories by 2
        X0_traj = (X0_traj[range(0,CurrentLength-1,2)] + X0_traj[range(1,CurrentLength,2)]) / 2
        X1_traj = (X1_traj[range(0,CurrentLength-1,2)] + X1_traj[range(1,CurrentLength,2)]) / 2
        
        CurrentLength = len(X0_traj)
        
        CurrentTimeSpacing  = 2 * CurrentTimeSpacing 
        
    return np.array(Correlation), np.array(Tau), np.array(DTau), np.array(Weight)

def Primed1DCorrelation(X1_traj, m, X0_traj, n):
    """
    Calculates a primed moment-correlation between two individual trajectories of equal length

    Parameters
    ----------
    X1_traj : Numpy array
        Delayed trajectory
    m : Int
        Order of primed power of X1_traj
    X0_traj : Numoy arrray
        First trajectory
    n : Int
        Order of primed power of X0_traj

    Returns
    -------
    Correlation : Numpy array
        1D array of correlation values as a function of delays in Tau
    Tau :  Numpy array
        1D array of delays corresponding to Correlation
    DTau : Numpy array
        1D array of time-bin widths corresponding to Tau
    Weight : Numpy array
        1D array of the number of original time points within each time bin

    """
    X1_traj = np.array(X1_traj)
    X0_traj = np.array(X0_traj)
    
    # Calculate X0nPrime
    if n == 0:
        X0nPrime = np.ones(len(X0_traj))
    else:
        X0nPrime = X0_traj
    
    if n >= 2:
        for i in range(2, n + 1):
            # Truncate the first point of X0nPrime and the last point of X0
            X0nPrime = X0nPrime[1:]
            X0_traj = X0_traj[:-1]
            
            X0nPrime = X0nPrime * X0_traj
    
    #  Calculate X1mPrime = X1 to the mʹ power
    
    if m == 0:
        X1mPrime = np.ones(len(X1_traj))
    else:
        X1mPrime = X1_traj
        
    # Truncate the first n points of X1mPrime and the last n points of X1
    X1mPrime = X1mPrime[n:]
    X1_traj = X1_traj[:-n]
    
    if m >= 2:
        for i in range(2, m + 1):
            X1mPrime = X1mPrime[1:]
            X1_traj = X1_traj[:-1]
            X1mPrime = X1mPrime * X1_traj
    
    return Correlation_1D(X1mPrime, X0nPrime[:len(X1mPrime)])


def  EnAvePrimedCorrelation(X1Ensemble, m, X0Ensemble, n):
    """
    Calculates a primed moment between two ensembles of trajectories. 
    Corresponding trajectories must have equal lengths, but different members of the ensemble 
    can have different lengths. 

    Parameters
    ----------
    X1Ensemble : List of Numpy arrays
        Structure of delayed trajectories
    m : Int
        order of primed power of X1
    X0Ensemble : List of Numpy arrays
        Structure of first trajectories
    n : Int
        order of primed power of X0

    Returns
    -------
    Correlation : Numpy array
        1D array of correlation values as a function of delays in Tau
    Tau :  Numpy array
        1D array of delays corresponding to Correlation
    DTau : Numpy array
        1D array of time-bin widths corresponding to Tau
    Weight : Numpy array
        1D array of the number of original time points within each time bin

    """
    SMStructure = [] # list containing the ensemble of single-molecule correlation structures
    SMWeight = [] # structure containing the ensemble of single-molecule correlation weights
    SMCorrelation = [] # structure containing the ensemble of single-molecule correlation functions
    
    
    MaxLength = 0 # length of longuest sm correlation
    MaxIndex = 0 # index of longuest sm correlation
    
    EnsembleSize = len(X0Ensemble) # number of molecules in the dataset
    
    # Create SMCorrelations and find the length and index of the longest single-molecule correlation 
    for i in range(EnsembleSize):
        SMStructure.append([Primed1DCorrelation(X1Ensemble[i], m, X0Ensemble[i], n)])
        CurrentLength = len(SMStructure[i][0][0])
        if MaxLength < CurrentLength:
            MaxLength = CurrentLength
            MaxIndex = i
    
    # Form rectangular arrays for correlations and weights by padding to the maximum length
    for i in range(EnsembleSize):
        PaddingLength = MaxLength - len(SMStructure[i][0][0]) 
        SMWeight.append(list(SMStructure[i][0][3]) + [0 for k in range(PaddingLength)])
        SMCorrelation.append(list(SMStructure[i][0][0]) + [0 for k in range(PaddingLength)])  # Padding value is unimportant, because the weights are zero
    
    SMWeight = np.array(SMWeight)
    SMCorrelation = np.array(SMCorrelation)
    
    # Calculate the weighted average of the correlation functions and sum of weights 
    # Correlation = np.sum(SMWeight * SMCorrelation, axis = 0) / EnsembleSize # original code from pseudocode, but I suspect an error with the weights
    Weight = np.sum(SMWeight, axis = 0)
    
    Correlation = np.sum(SMWeight * SMCorrelation, axis = 0) / Weight # my version of the above line, with correct weights; need to be confirmed
    
    # For Tau and DTau, take the longest examples
    Tau = SMStructure[MaxIndex][0][1]
    DTau = SMStructure[MaxIndex][0][2]
    
    # Assemble the results into a structure
    return [Correlation, Tau, DTau, Weight]
    
    