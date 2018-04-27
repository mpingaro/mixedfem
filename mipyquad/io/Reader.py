
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:13:01 2017

@author: Dr. Marco Pingaro 
@ref   : La Spienza University of Rome

Python Version 2.7
"""

import xml.etree.ElementTree as ET
import numpy as np
from sys import exit
import os

from VEMClasses import Analysis
from VEMClasses import GeometryMesh

class VEMReader:
        
    # Read input file xml
    @staticmethod
    def InputFile( name_f ):
        analysis = Analysis()

        if (os.path.exists(name_f)==False):
            exit ("Warning: File %s do not find" % (name_f)) 
            
        tree = ET.parse( name_f )
        root = tree.getroot()

        analysis.Title        = str( root[0][0].text )
        analysis.TypeAnalysis = str( root[0][1].text )
        analysis.TypeSolver   = str( root[0][2].text )
        analysis.dim_phy      = int( root[0][3].text )
        #
        #for mt in root.iter('Materials'):
        #    mat = int(mt.attrib['N'])

        for pt in root.iter('Material'): 
            analysis.TypeEl.append(str(pt.attrib['Type']))
            
            if str(pt.attrib['Type'])=='Cauchy':
                analysis.Materials.append([str(pt.attrib['Type']),
                    float(pt.attrib['Young']),float(pt.attrib['Poisson'])])
                
            elif str(pt.attrib['Type'])=='Orthotropic':
                analysis.Materials.append([str(pt.attrib['Type']),
                    float(pt.attrib['YoungX']),float(pt.attrib['YoungY']),
                    float(pt.attrib['PoissonX']),float(pt.attrib['PoissonY'])])
            else:
                exit ("Warning: type element do not exist")  
 
        analysis.fname_inp = str( root[2][0].text )
        analysis.fname_geo = str( root[2][1].text )
        analysis.fname_BCd = str( root[2][2].text )
        analysis.fname_res = str( root[2][3].text )
      
        return analysis

    # Read Geometry file xml 2/3-D problems
    @staticmethod
    def GeometryMesh(name_f,dim_phy,TypeEl):

        if (os.path.exists(name_f)==False):
            exit ("Warning: File %s do not find" % (name_f))
            
        tree = ET.parse( name_f )
        root = tree.getroot()
        
        msh = GeometryMesh()
        
        # Number of Nodes 
        for points in root.iter('Nodes'):
            msh.NNodes =(int(points.attrib['N']))
            
        msh.Nodes = np.zeros((msh.NNodes,dim_phy))
        i=0
        for point in root.iter('Node'):
            if dim_phy == 2:
                msh.Nodes[i][0] = float(point.attrib['X'])
                msh.Nodes[i][1] = float(point.attrib['Y'])
            elif dim_phy == 3:
                msh.Nodes[i][0] = float(point.attrib['X'])
                msh.Nodes[i][1] = float(point.attrib['Y'])
                msh.Nodes[i][3] = float(point.attrib['Z'])            
            i+=1
            
        # Number of Elements 
        for el in root.iter('Elements'):
            msh.NElements =(int(el.attrib['N']))
            
        msh.Elements    = []
        msh.IdEl        = np.zeros(msh.NElements)
        msh.DegElements = np.zeros(msh.NElements)
        msh.MAT         = np.zeros(msh.NElements)

        i=0
        for els in root.iter('Element'):

            msh.IdEl[i] = int(els.attrib['ElementId'])
            msh.MAT[i] = int(els.attrib['MAT'])

            # Degree
            msh.DegElements[i] = int(els.attrib['Degree'])
            
            # Element
            el = (els.attrib['Nodes'])
            el = el.split(' ')
            
            elm = np.zeros(len(el))
            for j in range(len(el)):
                elm[j] = int(el[j])
            msh.Elements += [elm]              
            
            i+=1
                    
        return msh
        
    # Read Boundary conditions file xml 2/3-D problems
    @staticmethod    
    def BoundaryConditions(name_f,dim_phy,msh):
        
        if (os.path.exists(name_f)==False):
            exit ("Warning: File %s do not find" % (name_f))
        
        tree = ET.parse( name_f )
        root = tree.getroot()
                
        msh.BCd = np.zeros((10000,2))
        msh.BCn = np.zeros((10000,2))
        msh.Ang = []
        # Dirichlet boundary conditions
        i=0
        for node in root.iter('Displacement'):
    
            Node = int(node.attrib['NodeId'])
            Comp = list(node.attrib['Components'])
            msh.Ang.append( (Node, float(node.attrib['Angle'])) )
            val  = node.text  
            val  = val.split(' ')
    
            for j in range(len(Comp)):
                
                if dim_phy==2:
                
                    if int(Comp[j])==0:
                        msh.BCd[i][0] = 2*Node
                    elif int(Comp[j])==1:
                        msh.BCd[i][0] = 2*Node+1
                    
                elif dim_phy==3:
                
                    if int(Comp[j])==0:
                        msh.BCd[i][0] = 3*Node
                    elif int(Comp[j])==1:
                        msh.BCd[i][0] = 3*Node+1
                    elif int(Comp[j])==2:
                        msh.BCd[i][0] = 3*Node+2
        
                msh.BCd[i][1] = float(val[j])
                i+=1
    
        msh.BCd = msh.BCd[0:i][:] 
        
        # Neumann boundary conditions
        i=0
        for force in root.iter('Force'):
    
            Node = int(force.attrib['NodeId'])
            Comp = list(force.attrib['Components'])
            val  = force.text  
            val  = val.split(' ')
    
            for j in range(len(Comp)):

                if dim_phy==2:
                
                    if int(Comp[j])==0:
                        msh.BCn[i][0] = 2*Node
                    elif int(Comp[j])==1:
                        msh.BCn[i][0] = 2*Node+1
        
                elif dim_phy==3:
                
                    if int(Comp[j])==0:
                        msh.BCn[i][0] = 3*Node
                    elif int(Comp[j])==1:
                        msh.BCn[i][0] = 3*Node+1
                    elif int(Comp[j])==2:
                        msh.BCn[i][0] = 3*Node+2
        
                msh.BCn[i][1] = float(val[j])
                i+=1
    
        msh.BCn = msh.BCn[0:i][:] 
        
        return msh
