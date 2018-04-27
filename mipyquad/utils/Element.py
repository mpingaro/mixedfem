
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 14:20:23 2017

@author: Dr. Marco Pingaro 
@ref   : La Spienza University of Rome

Python Version 2.7
"""
import numpy as np
from sys import exit
from VEMClasses import femEl
from Functions import Functions

class Element:
    
    @staticmethod
    def LocalElement(msh,analysis,nel):
        fem = femEl()
        dim_phy = analysis.dim_phy
        #
        fem.IdEl = msh.IdEl[nel]
        fem.ElementEl = msh.Elements[nel]
        fem.DegEl = msh.DegElements[nel]
        fem.NPointsEl = len(fem.ElementEl)
        
        if fem.DegEl==1:
            fem.NVertexEl = len(fem.ElementEl)
            fem.Ndofs = dim_phy*fem.NPointsEl
        elif fem.DegEl==2:
            fem.NVertexEl = len(fem.ElementEl)/2
            fem.Ndofs = dim_phy*(fem.NPointsEl + 1)
        elif fem.DegEl==3:
            fem.NVertexEl = len(fem.ElementEl)/3
            fem.Ndofs = dim_phy*(fem.NPointsEl + 3)
        else:
            exit ("Warning: Degree element %d do not valid" % (nel))
        #                                 
        fem.VertexEl = np.zeros((fem.NVertexEl,dim_phy))
        for i in range(fem.NVertexEl):
            for j in range(dim_phy):
                fem.VertexEl[i][j] = msh.Nodes[int(fem.ElementEl[i])][j]
        #
        fem.PointsEl = np.zeros((fem.NPointsEl,dim_phy))
        for i in range(fem.NPointsEl):
            for j in range(dim_phy):
                fem.PointsEl[i][j] = msh.Nodes[int(fem.ElementEl[i])][j]        

        fem.NEdgesEl   = fem.NVertexEl
        if fem.DegEl==1:
            fem.EdgesEl = np.zeros((fem.NEdgesEl,2))
            for i in range(fem.NEdgesEl):
                
                if i==fem.NEdgesEl-1:
                    fem.EdgesEl[i][0] = fem.ElementEl[i]
                    fem.EdgesEl[i][1] = fem.ElementEl[0]
                else:
                    fem.EdgesEl[i][0] = fem.ElementEl[i]
                    fem.EdgesEl[i][1] = fem.ElementEl[i+1]
        
        elif fem.DegEl==2:
            fem.EdgesEl = np.zeros((fem.NEdgesEl,3))
            for i in range(fem.NEdgesEl):
                
                if i==fem.NEdgesEl-1:
                    fem.EdgesEl[i][0] = fem.ElementEl[i]
                    fem.EdgesEl[i][1] = fem.ElementEl[0]
                    fem.EdgesEl[i][2] = fem.ElementEl[fem.NEdgesEl+i]
                else:
                    fem.EdgesEl[i][0] = fem.ElementEl[i]
                    fem.EdgesEl[i][1] = fem.ElementEl[i+1]
                    fem.EdgesEl[i][2] = fem.ElementEl[fem.NEdgesEl+i]

        elif fem.DegEl==3:
            fem.EdgesEl = np.zeros((fem.NEdgesEl,4))
            for i in range(fem.NEdgesEl):
                
                if i==fem.NEdgesEl-1:
                    fem.EdgesEl[i][0] = fem.ElementEl[i]
                    fem.EdgesEl[i][1] = fem.ElementEl[0]
                    fem.EdgesEl[i][2] = fem.ElementEl[fem.NEdgesEl+2*i]
                    fem.EdgesEl[i][3] = fem.ElementEl[fem.NEdgesEl+2*i+1]
                else:
                    fem.EdgesEl[i][0] = fem.ElementEl[i]
                    fem.EdgesEl[i][1] = fem.ElementEl[i+1]
                    fem.EdgesEl[i][2] = fem.ElementEl[fem.NEdgesEl+2*i]        
                    fem.EdgesEl[i][3] = fem.ElementEl[fem.NEdgesEl+2*i+1]        
        else:
            exit ("Warning: Degree element %d do not valid" % (nel))


        fem.NNormalEl  = fem.NVertexEl
        fem.NormalEl   = np.zeros((fem.NNormalEl,dim_phy))
        fem.LEdgesEl   = np.zeros(fem.NNormalEl)
        
        for nedg in range(fem.NNormalEl):
            ptA = msh.Nodes[int(fem.EdgesEl[nedg][1])]
            ptB = msh.Nodes[int(fem.EdgesEl[nedg][0])]
            
            fem.LEdgesEl[nedg] = np.sqrt( (ptA[0]-ptB[0])**2 + (ptA[1]-ptB[1])**2 )
            # Tangent(i,[1 2]) = (ptA-ptB)./len_edges(i);
            fem.NormalEl[nedg][0] = (ptA[1]-ptB[1])/fem.LEdgesEl[nedg]
            fem.NormalEl[nedg][1] = (ptB[0]-ptA[0])/fem.LEdgesEl[nedg]   
                 
        # TODO: this part not implemented (Only for 3-D)
        fem.FacesEl    = []
        fem.NFacesEl   = []

        # 
        fem.Area       = msh.Areas[nel]
        fem.B          = msh.Centroid[nel]
        # Compute fem.He 
        fem            = Functions.DiameterPolyElements(fem)
        
        # Compute Integral of monomials (Only for degree >= to 2)
        if fem.DegEl>=2:
            fem = Functions.IntegralMonomials(msh,fem)
        
        # Material part (Elastic moduli) 
        mat = int(msh.MAT[nel])

        if analysis.Materials[mat][0]=='Cauchy':
            fem.YoungEl    = analysis.Materials[mat][1]
            fem.PoissonEl  = analysis.Materials[mat][2]       
        elif analysis.TypeEl=='Orthotropic':              
            fem.YoungXEl      = analysis.Materials[mat][1]
            fem.YoungYEl      = analysis.Materials[mat][2]
            fem.PoissonXEl    = analysis.Materials[mat][3]
            fem.PoissonYEl    = analysis.Materials[mat][4]        
        else:
            exit ("Warning: type element do not exist")        
            
        return fem
