"""
EQ Shelving Filtering using numpy

Author: rasmushyldgaard

Date: 20-05-2023
"""

import numpy as np


def shelving_filter(gain: int,
                    fc: int,
                    fs: int,
                    mode: str = 'low' or 'high'):
    """
    DAFX: Digital Audio Effects by Udo Zolzer (Chapter 2.3.)
    
    :param gain              Gain in dB
    :param fc                Cutoff frequency in Hz
    :param fs                Sampling frequency in Hz
    :param mode              Shelving to apply: 'low' shelf
                             or 'high' shelf
    
    :return                  Filter coefficients
    """
    if mode != 'low' and mode != 'high':
        print("Shelving mode is not recognized! Please use either 'low' or 'high'...")
        exit(1)
    
    # filter coefficients
    a = np.empty((1, 3))
    b = np.empty((1, 3))
    
    if gain != 0:
        K = np.tan((np.pi * fc) / fs)
        V0 = 10**(gain / 20)
        den1 = 1 + np.sqrt(2) * K + K**2
        
        if mode == 'low':
            b, a = low_shelf(gain=gain,
                             K=K,
                             V0=V0,
                             den1=den1)
        else:
            b, a = high_shelf(gain=gain,
                              K=K,
                              V0=V0,
                              den1=den1)
            
    else:
        b = np.array([1, 1, 1])
        a = np.array([1, 1, 1])
        
    return b, a
    
    
def low_shelf(gain, K, V0, den1):
    if gain > 0:
        b0 = (1 + (np.sqrt(2*V0) * K) + V0 * K**2) / (den1)
        b1 = (2 * ((V0 * K**2)-1)) / (den1)
        b2 = (1 - (np.sqrt(2*V0) * K) + V0 * K**2) / (den1)
        a1 = (2 * (K**2 - 1)) / (den1)
        a2 = (1 - np.sqrt(2) * K + K**2) / (den1)
    
    else:
        den2 = V0 + np.sqrt(2*V0) * K + K**2
        b0 = (V0 * (1 + np.sqrt(2) * K + K**2)) / (den2)
        b1 = (2 * V0 * (K**2 - 1)) / (den2)
        b2 = (V0 * (1 - np.sqrt(2) * K + K**2)) / (den2)
        a1 = (2 * (K**2 - V0)) / (den2)
        a2 = (V0 - np.sqrt(2*V0) * K + K**2) / (den2)
    
    b = np.array([b0, b1, b2])
    a = np.array([1, a1, a2])
    return b, a


def high_shelf(gain, K, V0, den1):
    if gain > 0:
        b0 = (V0 + np.sqrt(2*V0) * K + K**2) / (den1)
        b1 = (2 * (K**2 - V0)) / (den1)
        b2 = (V0 - np.sqrt(2*V0) * K + K**2) / (den1)
        a1 = (2 * (K**2 - 1)) / (den1)
        a2 = (1 - np.sqrt(2) * K + K**2) / (den1)
        
    else:
        den2 = 1 + np.sqrt(2*V0) * K + V0 * K**2
        b0 = (V0 * (1 + np.sqrt(2) * K + K**2)) / (den2)
        b1 = (2 * V0 * (K**2 - 1)) / (den2)
        b2 = (V0 * (1 - np.sqrt(2) * K + K**2)) / (den2)
        a1 = (2 * (V0 * K**2 - 1)) / (den2)
        a2 = (1 - np.sqrt(2*V0) * K + V0 * K**2) / (den2)
        
    b = np.array([b0, b1, b2])
    a = np.array([1, a1, a2])
    return b, a