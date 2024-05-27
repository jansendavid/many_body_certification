#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 14:44:43 2024

@author: david
"""

import numpy as np
sym=["z", "x", "y"]
sx=np.array([[0,1],[1,0]],dtype=complex)
sz=np.array([[1,0],[0,-1]],dtype=complex)
sy=np.array([[0,-1j],[1j,0]],dtype=complex)
I=np.eye(2, dtype=complex)