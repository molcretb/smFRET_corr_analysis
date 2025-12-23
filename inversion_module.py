# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 16:10:58 2025

@author: molcre0000
"""

from correlation_module import *
import numpy as np
import cvxpy as cp



def CreateDerivative(ProbPoints):
    """
    Calculate the first derivative matrix D defined in Eq.(58) of https://doi.org/10.1063/5.0284658 and used in InvertMoments

    Parameters
    ----------
    ProbPoints : Numpy array
        1D array of points where the probability distribution will be calculated

    Returns
    -------
    D :  Numpy array
        2D array representing the derivative matrix 

    """
    
    R = len(ProbPoints)
    D = np.zeros((R-1, R))
    for i in range(R-1):
        for j in range(R):
            if i == j:
                D[i, j] = 1 # Diagonal elements
            elif i+1 == j:
                D[i, j] = -1  #  Off-diagonal points
    return D

def CreateVandermonde(ProbPoints, InversionOrder):
    """
    Calculate the Vandermonde matrix V defined in Eq. (54) of https://doi.org/10.1063/5.0284658 and used in InvertMoments

    Parameters
    ----------
    ProbPoints : Numpy array
        1D array of points where the probability distribution will be calculated
    InversionOrder : Int
        Maximum moment-order to be calculated

    Returns
    -------
    V : Numpy array
        2D array representing the Vandermonde matrix

    """
    R = len(ProbPoints)
    
    V = np.zeros((InversionOrder, R))
    
    for i in range(0, InversionOrder):
        V[i, :] = ProbPoints**(i+1)
    return V
    

def InvertMoments(Moments, Weights, ProbPoints, V, D, Beta, InitialProb):
    """
    Inverts a sequence of univariate moments to its probability distribution for a given smoothing parameter.
    See Eq. (55) of https://doi.org/10.1063/5.0284658

    Parameters
    ----------
    Moments : Numpy array
        1D array of moments, starting at first order
    Weights : Numpy array
        1D array of estimated weights for each moment 
    ProbPoints : Numpy array
        1D array of the points where the probability will be calculated
    V : Numpy array
        2D array containing the Vandermonde matrix
    D : Numpy array
        2D array containing the first derivative matrix
    Beta : float
        scalar value of the smoothing parameter
    InitialProb : Numpy array
        1D array of initial probabilities

    Returns
    -------
    Probabilities : Numpy array
        1D array of probabilities corresponding to ProbPoints
    Chi : float
        scalar measuring the fit to the moments
    Smoothness : float
        scalar measuring the smoothness of the probability distribution

    """
    # Combine the two terms in the minimization with Weights (original from pseudocode)
    #C = np.concat((Weights.reshape(-1, 1)*V, Beta * D)) # Concatenate coefficients
    #d = np.concat((Weights * Moments, np.zeros((len(ProbPoints)-1))))
    
    # BM corrected
    new_V = V/Moments.reshape(-1, 1)
    new_V[~np.isfinite(new_V)] = 0
    
    C = np.concat((new_V, Beta * D)) # Concatenate coefficients
    d = np.concat((np.ones(len(Moments)), np.zeros((len(ProbPoints)-1))))
    
    # Constraints on the minimization
    Aeq = np.ones(len(ProbPoints)) # Linear constraint. Coefficient in sum of Probabilities
    beq = 1 # Linear constraint constant = Total probability 
    # LowerBound = np.zeros((len(ProbPoints), 1)) # Boundary constraint = Positive probabilities
    
    # Quadratic minimization, see equation (55), using CVXPY library
    x = cp.Variable(len(ProbPoints))
    objective = cp.Minimize(cp.sum_squares(C @ x - d))
    constraints = [0 <= x, 1 >= x, np.matmul(Aeq,x) == beq]
    
    
    
    prob = cp.Problem(objective, constraints)
    prob.solve(verbose= True, max_iter=100000, solver=cp.CLARABEL)
    Probabilities = x.value
    
    
    # mean-squared fractional deviation of the solution moments from the input moments, see equation (56)
    # Chi = (1/np.sqrt(len(Moments))) * np.linalg.norm(Weights.reshape(-1, 1)*(np.matmul(V, Probabilities) - Moments.reshape(-1, 1))) # from pseudocode
    Chi = (1/np.sqrt(len(Moments))) * np.linalg.norm(np.matmul(V, Probabilities)[1:]/Moments[1:] - 1)

    # scalar measuring the smoothness of the probability distribution, see equation (57)
    Smoothness = (1/np.sqrt(len(ProbPoints)-1)) * np.linalg.norm(np.matmul(D, Probabilities))
    
    return Probabilities, Chi, Smoothness
    
    