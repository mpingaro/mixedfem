
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 12:14:58 2017

@author: Dr. Marco Pingaro 
@ref   : La Spienza University of Rome

Python Version 2.7
"""
# Classes used in PyVEM


# Contains values of input data
class Analysis():
    
    def __init__(self):
        self.Title        = None
        self.TypeAnalysis = None
        self.TypeSolver   = None
        self.dim_phy      = None
        #
        self.Materials    = []
        self.TypeEl       = []
        #
        self.fname_inp    = None
        self.fname_res    = None
        
# Contains reading of the mesh data        
class GeometryMesh():
    
    def __init__(self):
        self.Nodes       = None
        self.NNodes      = None
        self.IdEl        = None
        self.Elements    = None
        self.NElements   = None
        self.DegElements = None
        self.MAT         = None
        self.Centroid    = None
        self.Areas       = None        
        self.BCd         = None
        self.Ang         = None
        self.BCn         = None
        
        
# Contains local geometry and material moduli        
class femEl():
    
    def __init__(self):
        self.VertexEl     = None
        self.NVertexEl    = None
        self.PointsEl     = None
        self.NPointsEl    = None
        self.IdEl         = None
        self.DegEl        = None
        self.Ndofs        = None
        self.NormalEl     = None
        self.NNormalEl    = None
        self.ElementEl    = None
        self.EdgesEl      = None
        self.LEdgesEl     = None
        self.NEdgesEl     = None
        self.FacesEl      = None
        self.NFacesEl     = None
        self.B            = None
        self.Area         = None
        self.He           = None
        self.IntMonomials = None
        self.YoungEl      = None
        self.PoissonEl    = None
        self.YoungXEl     = None
        self.PoissonXEl   = None
        self.YoungYEl     = None
        self.PoissonYEl   = None

# Contains the results structure         
class Results():
    
    def __init__(self):
        self.displacement  = None
        self.deformation   = None
        self.stresses      = None
        self.von_mises     = None
        self.strains       = None
        self.values        = None
        self.energyComp    = None
        self.energyEl      = None
        self.energyTot     = None
        self.energyTotComp = None