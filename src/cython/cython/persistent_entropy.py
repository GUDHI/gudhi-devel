#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    
    Persistent entropy is a function from barcodes to the reals. It can be used
    to perform statistical tests on barcodes. It was defined in [1], although
    a precursor definition appears in [2].
    
    [1] M. Rucco, F. Castiglione, E. Merelli, M. Pettini, Characterisation of 
    the idiotypic immune network through persistent entropy, in: Proc. Complex,
    2015.
    [2] H. Chintakunta, T. Gentimis, R. Gonzalez-Diaz, M.-J. Jimenez,
    H. Krim, An entropy-based persistence barcode, Pattern Recognition
    48 (2) (2015) 391–401.
        
    Implementation of persistent entropy
    Author: Manuel Soriano-Trigueros (cimagroup, University of Seville)
    contact: msoriano4@us.es
    
"""

import numpy as np

def persistent_entropy(
    dgms, val_dim=None, val_inf=None, normalize=False, basis=2
):
    """
    Calculate persistent entropy values of a family of persistence barcodes
    or diagrams. If no dimension is specified, it will consider all bars in the
    barcode. If the infinity bars have any meaning in your experiment and you 
    want to use them for the calculations, remember to give the value you 
    desire to val_Inf.
    Parameters
    -----------
    dgms: list of diagrams, each of them formed by a list of tuples 
        (dimension, (birth, death)).
    val_dim: int, default None.
        calculate persistent entropy only for bars with a fixed 
        dimension. If this value is None, calculate persistent entropy for all  
        the bars.
    val_inf: float, default None
        value infinity bars will be given. If it is None, infinity bars are 
        ignored.
    normalize: bool, default False
        if False, the persistent entropy values are not normalized.
        if True, the persistent entropy values are normalized (not recommended
        for statistical tests).
          
    Returns
    --------
    p_ent: ndarray (n_pairs,)
        array of persistent entropy values corresponding to each barcode.
    """
    
    # Step 0: If the input is only one barcode, put it in a list.
    if isinstance(dgms[0], list) == False:
        dgms = [dgms]
        
    # Step 1: extract the barcode of dimension we are working with as an array.
    if val_dim == None:
        dgms = [np.array([[element[1][0], element[1][1]]
                   for element in dgm])
                    for dgm in dgms]
    else:
        dgms = [np.array([[element[1][0], element[1][1]] 
                   for element in filter(lambda x: x[0]==val_dim, dgm)])
                    for dgm in dgms]
        
    # Step2: remove infinity bars or change them by val_inf.
    if val_inf==None:
        dgms = [np.stack(filter(lambda x: x[1]!=np.infty, dgm)) for dgm in dgms]
    else:
        dgms = [np.where(dgm == np.inf, val_inf, dgm) for dgm in dgms]
        
    # Step3: calculate persistent entropy
    p_ent = []
    for dgm in dgms:
        l = dgm[:, 1] - dgm[:, 0]
        if all(l > 0):
            L = np.sum(l)
            p = l / L
            E = -np.sum(p * (np.log2(p)/np.log2(basis)))
            if normalize == True:
                E = E / (np.log2(len(l))/np.log2(basis))
            p_ent.append(E)
        else:
            raise Exception("A bar has negative length")

    return np.array(p_ent)
    