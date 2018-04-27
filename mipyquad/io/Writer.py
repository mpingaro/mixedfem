
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:18:50 2017

@author: Dr. Marco Pingaro 
@ref   : La Spienza University of Rome

Python Version 2.7
"""
#from tabulate import tabulate

class VEMWriter:
    
    @staticmethod
    def GeometryVTK(msh,path,analysis,results):
        name_f = path.OUT + '/' + analysis.fname_res + '.0.000.vtk'  
        f = open(name_f, 'w')
                
        f.write("# vtk DataFile Version 2.0\n")
        f.write(" Marco Pingaro data ...\n") 
        f.write("ASCII\n") 
        f.write("DATASET POLYDATA\n") 
        f.write("POINTS    %d float\n" % msh.NNodes)
        for i in range(msh.NNodes):
            f.write("   %.8f   %.8f   %.8f\n" % (msh.Nodes[i][0], msh.Nodes[i][1],0.0))
        
        n = 0
        for i in range(msh.NElements):
            n += len(msh.Elements[i]) 
        n += msh.NElements
        
        f.write("POLYGONS %d %d\n" % (msh.NElements,n)) 
        for nl in range(msh.NElements):
            nl_e = len(msh.Elements[nl])
            f.write("%d " % nl_e)
            for nd in range(nl_e):
                El = msh.Elements[nl][nd]
                f.write("%d " % El)
            f.write("\n")
        f.write("\n")
        f.close()
 
    # TODO: change this part to plot in the same file geometry and deformation
    #       and the stresses --> for 3-D plotting
      
    @staticmethod
    def ResultsVTK(msh,path,analysis,results):
        
        name_f = path.OUT + '/' + analysis.fname_res + '.vtk'
        f = open(name_f, 'w')
                
        f.write("# vtk DataFile Version 2.0\n")
        f.write(" Marco Pingaro data ...\n") 
        f.write("ASCII\n") 
        f.write("DATASET POLYDATA\n") 
        f.write("POINTS    %d float\n" % msh.NNodes)
        for i in range(msh.NNodes):
            f.write("   %.8f   %.8f   %.8f\n" % (msh.Nodes[i][0], msh.Nodes[i][1],0.0))
        n = 0
        for i in range(msh.NElements):
            n += len(msh.Elements[i]) 
        n += msh.NElements
        
        f.write("POLYGONS %d %d\n" % (msh.NElements,n)) 
        for nl in range(msh.NElements):
            nl_e = len(msh.Elements[nl])
            f.write("%d " % nl_e)
            
            if msh.DegElements[nl]==1:
                for nd in range(nl_e):
                    El = msh.Elements[nl][nd]
                    f.write("%d " % El)
            
            elif msh.DegElements[nl]==2:
                for nd in range(nl_e/2):
                    El_vx = msh.Elements[nl][2*nd]
                    El_pm = msh.Elements[nl][2*nd+1]
                    f.write("%d %d" % (El_vx,El_pm))
            
            elif msh.DegElements[nl]==3:
                for nd in range(nl_e/3):
                    El_vx = msh.Elements[nl][3*nd]
                    El_pm1 = msh.Elements[nl][3*nd+1]
                    El_pm2 = msh.Elements[nl][3*nd+2]
                    f.write("%d %d %d" % (El_vx,El_pm1,El_pm2))   
                    
            f.write("\n")
        f.write("\n")
        
        # Cells 
        f.write("CELL_DATA %d\n" % (msh.NElements))
        f.write("SCALARS cell_scalars int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for i in range(msh.NElements):
            f.write("%d\n" % i)
        f.write("\n")    
        
        # Displacement X,Y,Z
        f.write("POINT_DATA %d\n" % msh.NNodes)
        f.write("VECTORS Displacement float\n")
        #f.write("LOOKUP_TABLE default\n")
        for i in range(msh.NNodes):        
            f.write("%0.8f %0.8f %0.8f\n" % (results.displacement[i][0], results.displacement[i][1], 0.0))        
        
        # Von Mises Stress
        f.write("CELL_DATA %d\n" % (msh.NElements))
        f.write("SCALARS Von-Mises float\n")
        f.write("LOOKUP_TABLE default\n")
        for i in range(msh.NElements):        
            f.write("%0.8f\n" % (results.von_mises[i][0])) # Change this!!!
        
        # Stresses
        f.write("TENSORS Cauchy-Stress float\n")  
        for i in range(msh.NElements):
            f.write("%.8f %.8f %.8f\n" % (results.stresses[i][0], results.stresses[i][2], 0.0) )
            f.write("%.8f %.8f %.8f\n" % (results.stresses[i][2], results.stresses[i][1], 0.0) )
            f.write("%.8f %.8f %.8f\n" % (0.0, 0.0, 0.0) )
            f.write("\n")
          
        # Strains
        f.write("TENSORS Strain float\n")
        for i in range(msh.NElements):
            f.write("%.8f %.8f %.8f\n" % (results.strains[i][0], results.strains[i][2], 0.0) )
            f.write("%.8f %.8f %.8f\n" % (results.strains[i][2], results.strains[i][1], 0.0) )
            f.write("%.8f %.8f %.8f\n" % (0.0, 0.0, 0.0) )
            f.write("\n")
        
        f.close()
        
        
    @staticmethod
    def IntegralResultsTXT(msh,path,analysis,results):
        name_f = path.OUT + '/' + analysis.fname_res + '.txt'
        f = open(name_f, 'w')
        
        f.write("Compute the integrals of the compliance for all elements\n")
        f.write("\n")
        f.write("Elements Id \t   Energy Comp_xx \t   Energy Comp_yy \t   Energy Comp_xy \t   Energy Comp_yx \t   Total Energy\n")
        f.write("\n")
        for i in range(msh.NElements):
            f.write("\t %d \t\t\t\t %.4e \t\t %.4e \t\t %.4e \t\t %.4e \t\t %.4e \n" 
                    % (i,results.energyComp[i][0],
                       results.energyComp[i][1],
                       results.energyComp[i][2],
                       results.energyComp[i][3], results.energyEl[i]))
        f.write("\n")
        f.write("Total Energy Components: %.4e \t\t %.4e \t\t %.4e \t\t %.4e \n" 
                % (results.energyTotComp[0],results.energyTotComp[1],results.energyTotComp[2],results.energyTotComp[3]))
        f.write("Total Energy Problem: %.4e\n" % results.energyTot)
        
#        for i in range(msh.NElements):
#            f.write("Elements Id %d \t %.2f\n" % (i,results.stresses[i][0]))
            
        f.close()
        
    @staticmethod
    def ResultsTXT(msh,path,analysis,results):
        name_f = path.OUT + '/' + analysis.fname_res + '.txt'
        f = open(name_f, 'w')
        
        f.write("Compute Area, Stress and Strain components\n")
        f.write("\n")
        f.write("Elements Id \t Areas \t S_xx \t S_yy \t S_xy \t S_yx \t E_xx \t E_yy \t E_xy \t E_yx \n")
        f.write("\n")
        for i in range(msh.NElements):
            f.write("\t %d \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \n" 
                    % (i,msh.Areas[i],results.stresses[i][0],results.stresses[i][1],
                       results.stresses[i][2],results.stresses[i][2],
                       results.strains[i][0],results.strains[i][1],
                       results.strains[i][2],results.strains[i][2]))

        f.close()        
        