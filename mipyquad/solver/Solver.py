
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 15:39:25 2017

@author: Dr. Marco Pingaro 
@ref   : La Spienza University of Rome

Python Version 2.7
"""

import numpy as np
import time
from scipy.sparse import csc_matrix, linalg as sla

# Using Trilinos (Iterative solver)
#from PyTrilinos import Epetra
#from PyTrilinos import AztecOO


class Solver:
    
    @staticmethod 
    def LinearSolver(options,msh,analysis,stiffK,rhs):
        
        if options.silent==False:
            print("Solving linear system")
            start_time = time.time()

        # Local system for the nodes
        index_lc_pt = [] #
        for pt in range(len(msh.Ang)):
            if msh.Ang[pt][1] != 0.0:
                
                index_lc_pt.append(pt)
                
                node  = int( msh.Ang[pt][0] )
                theta = np.deg2rad( msh.Ang[pt][1] )
                
                gdl = 2*node
                
                # Create rotational matrix relative to the i-th node
                T = np.zeros((2,2))
                T[0][0] = np.cos(theta)
                T[0][1] = -np.sin(theta)
                T[1][0] = np.sin(theta)
                T[1][1] = np.cos(theta)  
                # Rotation of stiffness matrix  
                stiffK[[gdl,gdl+1],:] = np.dot(np.transpose(T), stiffK[[gdl,gdl+1],:])
                stiffK[:,[gdl,gdl+1]] = np.dot(stiffK[:,[gdl,gdl+1]], T)
                # Rotation of load vector
                rhs[gdl:gdl+2] = np.dot( np.transpose(T), rhs[gdl:gdl+2] )
                        
        # Impose Dirichlet Boundary Conditions
        ngdl  = len(stiffK)
        n_BCd = len(msh.BCd) 
        
        for n in range(n_BCd):
            stiffK[int(msh.BCd[n][0])][:] = np.zeros(ngdl)
            stiffK[int(msh.BCd[n][0])][int(msh.BCd[n][0])] = 1. 
            rhs[int(msh.BCd[n][0])] = msh.BCd[n][1]
        
        ##
        # Solve linear system 
        ##
        
        # Use Trilinos AztecOO
#        comm = Epetra.PyComm() # Comunications
#        numRows = ngdl
#        stdMap = Epetra.Map(numRows, 0, comm)        
#        K = K = Epetra.SerialDenseMatrix(stiffK)
#        
#        U   = Epetra.Vector(stdMap)
#        rhs = Epetra.Vector(rhs)
#        
#        ## Iteration solver
#        linearProblem = Epetra.LinearProblem(K, U, rhs)
#        solver = AztecOO.AztecOO(linearProblem)
#        solver.Iterate(10000,1.e-5)


        # Use superLU decompositions    
        K = csc_matrix(stiffK) 

        # SuperLU solver
        if analysis.TypeSolver == 'SuperLU':
            lu = sla.splu(K)
            U = lu.solve(rhs)
        
        # Krylov subspace method
        if analysis.TypeSolver == 'Krylov':
            U = sla.minres(K,rhs, x0=np.zeros(ngdl), shift=0.0, 
                           tol=1e-05, maxiter=None, xtype=None, 
                           M=None, callback=None, show=False, check=False)
            U = U[0]
        
        # Numpy linalg
        if analysis.TypeSolver == 'Numpy':
            U = np.linalg.solve(stiffK, rhs) # np linal


        # Recovery global system
        for i in range(len(index_lc_pt)):
            gdl = 2* int( msh.Ang[index_lc_pt[i]][0] )
            theta = np.deg2rad( msh.Ang[index_lc_pt[i]][1] )
            
            p = U[gdl]
            q = U[gdl+1]
            
            t11 = np.cos(theta)
            t12 = -np.sin(theta)
            t21 = np.sin(theta)
            t22 = np.cos(theta)
            
            U[gdl]   = p*t11 + q*t12
            U[gdl+1] = p*t21 + q*t22
        
        if options.silent==False:
            end_time = time.time()
            secs = end_time - start_time
            print("  Elapsed time .................................: %.2e Seconds" % secs)
   
        return U
