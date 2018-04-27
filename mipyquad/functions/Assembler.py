
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 21:37:23 2017

@author: Dr. Marco Pingaro 
@ref   : La Spienza University of Rome

Python Version 2.7
"""
import numpy as np
from LinearElasticVEM import LinearElasticVEM
from Element import Element
from Load import Load
import time

class Assembler:
    
    @staticmethod
    def local_to_globalK(localK,fem,analysis,Ngdl_p,stiffK):
        
        dim_phy = analysis.dim_phy  # Physic dimension of domain
        
        if fem.DegEl==1:
            mc_lg = np.zeros(dim_phy*fem.NPointsEl)
        elif fem.DegEl==2:
            mc_lg = np.zeros(dim_phy*(fem.NPointsEl+1))
        elif fem.DegEl==3:
            mc_lg = np.zeros(dim_phy*(fem.NPointsEl+3))
        
        for k in range(fem.NPointsEl):
            
            if dim_phy==2:
                mc_lg[dim_phy*k]   = fem.ElementEl[k]*dim_phy 
                mc_lg[dim_phy*k+1] = fem.ElementEl[k]*dim_phy + 1             
            
            if dim_phy==3:
                mc_lg[dim_phy*k]   = fem.ElementEl[k]*dim_phy 
                mc_lg[dim_phy*k+1] = fem.ElementEl[k]*dim_phy + 1
                mc_lg[dim_phy*k+2] = fem.ElementEl[k]*dim_phy + 2
        # Add internal degrees of freedom (Degree > 1)
        if fem.DegEl>1:
            if dim_phy==2:
                if fem.DegEl==2:
                    mc_lg[fem.Ndofs-2] = Ngdl_p + fem.IdEl*dim_phy
                    mc_lg[fem.Ndofs-1] = Ngdl_p + fem.IdEl*dim_phy + 1
                elif fem.DegEl==3:
                    mc_lg[fem.Ndofs-6] = Ngdl_p + fem.IdEl*dim_phy
                    mc_lg[fem.Ndofs-5] = Ngdl_p + fem.IdEl*dim_phy +1
                    mc_lg[fem.Ndofs-4] = Ngdl_p + fem.IdEl*dim_phy +2
                    mc_lg[fem.Ndofs-3] = Ngdl_p + fem.IdEl*dim_phy +3
                    mc_lg[fem.Ndofs-2] = Ngdl_p + fem.IdEl*dim_phy +4
                    mc_lg[fem.Ndofs-1] = Ngdl_p + fem.IdEl*dim_phy +5
        
        # TODO: do not use the symmetry (change this part)
        for row in range(fem.Ndofs):
            for col in range(fem.Ndofs):
                stiffK[int(mc_lg[row])][int(mc_lg[col])] += localK[row][col]
                
        return stiffK
    
    @staticmethod 
    def local_to_global_rhs(local_rhs,fem,analysis,rhs):
        
        # TODO: try to implement      
      
        
        return rhs
    
    @staticmethod
    def StandardFEM(options,msh,analysis):
        
        if options.silent==False:
            print("Assemblying matrix")
            start_time = time.time()

        # Inizialize stiffK
        Ngdl_p = analysis.dim_phy*msh.NNodes
        t_ngdl = Ngdl_p
        for k in range(msh.NElements):
            if msh.DegElements[k]==2:
                t_ngdl += analysis.dim_phy
            elif msh.DegElements[k]==3:
                t_ngdl += analysis.dim_phy*3
        
        stiffK = np.zeros((t_ngdl,t_ngdl))
        
        # Store Pi nabla star of all elements
        PiN_s = []

        # Cycle on all Elements
        for nel in range(msh.NElements):
            
            # Compute structure fem
            fem = Element.LocalElement(msh,analysis,nel)
            
            # Global Stiffness matrix        
            localK, PiNs = LinearElasticVEM.StiffnessMatrix2D(fem,analysis)
            PiN_s.append(PiNs)
            
            stiffK = Assembler.local_to_globalK(localK,fem,analysis,Ngdl_p,stiffK)
            
            # Compute rhs body load
            # TODO: try to implement
        
        if options.silent==False:
            end_time = time.time()
            secs = end_time - start_time
            print("  Elapsed time .................................: %.2e Seconds" % secs)
        
        # Compute rhs external forces
        if options.silent==False:
            print("Assemblying load vector")
            start_time = time.time()
        
        rhs = Load.NeumannConditions(msh,analysis,t_ngdl)
        
        if options.silent==False:
            end_time = time.time()
            secs = end_time - start_time
            print("  Elapsed time .................................: %.2e Seconds" % secs) 
                     
        return stiffK, rhs, PiN_s
    
