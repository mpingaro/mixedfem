
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 16:46:30 2017

@author: Dr. Marco Pingaro 
@ref   : La Spienza University of Rome

Python Version 2.7
"""
import numpy as np
from scipy.linalg import block_diag

# Linear elastic VEM Element 

class LinearElasticVEM:
    
    #  2-D Case (Plane Strain)   
    @staticmethod
    def StiffnessMatrix2D(fem,analysis, alf=0.5):
        
        # Dimension of Physics domain
        #dim_phy = analysis.dim_phy

        # Linear Elastic
        T = fem.YoungEl/(2*(1+fem.PoissonEl))
        C =  [[T*2*(1-fem.PoissonEl)/(1-2*fem.PoissonEl), T*2*fem.PoissonEl/(1-2*fem.PoissonEl), 0],
               [T*2*fem.PoissonEl/(1-2*fem.PoissonEl), T*2*(1-fem.PoissonEl)/(1-2*fem.PoissonEl), 0], 
               [0, 0, T]]

        if fem.DegEl==1:        
            I = np.identity(3)
            G    = np.dot(fem.Area, I)
            invG = np.dot((1./fem.Area), I)
            B = np.zeros((3,fem.Ndofs))
            D = [] # np.zeros((fem.Ndofs,6))
            for i in range(fem.NVertexEl):
                # Compute B matrix
                if i==0: 
                    n = 0 
                    m = fem.NVertexEl-1
                else: 
                    n = i-1 
                    m = i
                
                bound_nx = 0.5*(fem.NormalEl[n][0]*fem.LEdgesEl[n] + fem.NormalEl[m][0]*fem.LEdgesEl[m])
                bound_ny = 0.5*(fem.NormalEl[n][1]*fem.LEdgesEl[n] + fem.NormalEl[m][1]*fem.LEdgesEl[m])
                #
                B[0][2*i]   = bound_nx 
                #
                B[1][2*i+1] = bound_ny;
                #
                B[2][2*i]   = bound_ny;
                B[2][2*i+1] = bound_nx;
                # Compute D matrix
                xi  = ( fem.VertexEl[i][0] - fem.B[0] )/fem.He
                eta = ( fem.VertexEl[i][1] - fem.B[1] )/fem.He             
                D.append([1, 0, xi, 0, eta, 0])
                D.append([0, 1, 0, xi, 0, eta])               
        
            # Compute Pi operator
            PN = np.dot( invG, B)
            # Compute consistent part of stiffness matrix
            M = np.dot( np.dot( np.transpose(PN),C), np.dot(G,PN) )
            # Compute stabilization part of stiffness matrix
            St = np.dot( np.dot(alf,np.trace(M)), np.add(np.identity(fem.Ndofs) ,
                            -np.dot(np.dot(D, np.linalg.inv(np.dot(np.transpose(D), D))), np.transpose(D))))
            # Stiffness matrix 
            localK = np.zeros(fem.Ndofs)        
            localK = np.add(M, St) 

            # Vectorize Pi nabla star
            PiN_s = np.reshape(PN,(3*fem.Ndofs))
        
        elif fem.DegEl==2:
            G    = [[fem.IntMonomials[0], 0, 0, fem.IntMonomials[1], 0, 0, fem.IntMonomials[2], 0, 0],
                    [0, fem.IntMonomials[0], 0, 0, fem.IntMonomials[1], 0, 0, fem.IntMonomials[2], 0],
                    [0, 0, fem.IntMonomials[0], 0, 0, fem.IntMonomials[1], 0, 0, fem.IntMonomials[2]],
                    [fem.IntMonomials[1], 0, 0, fem.IntMonomials[4], 0, 0, fem.IntMonomials[3], 0, 0],
                    [0, fem.IntMonomials[1], 0, 0, fem.IntMonomials[4], 0, 0, fem.IntMonomials[3], 0],
                    [0, 0, fem.IntMonomials[1], 0, 0, fem.IntMonomials[4], 0, 0, fem.IntMonomials[3]],
                    [fem.IntMonomials[2], 0, 0, fem.IntMonomials[3], 0, 0, fem.IntMonomials[5], 0, 0],
                    [0, fem.IntMonomials[2], 0, 0, fem.IntMonomials[3], 0, 0, fem.IntMonomials[5], 0],
                    [0, 0, fem.IntMonomials[2], 0, 0, fem.IntMonomials[3], 0, 0, fem.IntMonomials[5]]]
            invG = np.linalg.inv(G)
            B = np.zeros((9,fem.Ndofs))
            D = [] #np.zeros((fem.Ndofs,12))               
        
            for i in range(fem.NPointsEl): #
                # Compute D matrix
                xi  = ( fem.PointsEl[i][0] - fem.B[0] )/fem.He
                eta = ( fem.PointsEl[i][1] - fem.B[1] )/fem.He
                # Nodal points / Middle nodes
                D.append([1, 0, xi, 0, eta, 0, xi**2, 0, xi*eta, 0, eta**2, 0])
                D.append([0, 1, 0, xi, 0, eta, 0, xi**2, 0, xi*eta, 0, eta**2])
            
                # Compute B matrix
                if i<fem.NEdgesEl:
                    if i==0: 
                        n = 0 ; 
                        m = fem.NEdgesEl-1
                    else:
                        n = i-1
                        m = i

                    nx_n = fem.NormalEl[n][0] 
                    ny_n = fem.NormalEl[n][1]
                    nx_m = fem.NormalEl[m][0]
                    ny_m = fem.NormalEl[m][1]
             
                    Pa = fem.VertexEl[i]
                    #                   
                    B[0][2*i] = (fem.LEdgesEl[n]*nx_n + fem.LEdgesEl[m]*nx_m)/6
                    B[2][2*i] = (fem.LEdgesEl[n]*ny_n + fem.LEdgesEl[m]*ny_m)/6
                    B[3][2*i] = (fem.LEdgesEl[n]*nx_n + fem.LEdgesEl[m]*nx_m)*( Pa[0] - fem.B[0] )/(6*fem.He)
                    B[5][2*i] = (fem.LEdgesEl[n]*ny_n + fem.LEdgesEl[m]*ny_m)*( Pa[0] - fem.B[0] )/(6*fem.He)
                    B[6][2*i] = (fem.LEdgesEl[n]*nx_n + fem.LEdgesEl[m]*nx_m)*( Pa[1] - fem.B[1] )/(6*fem.He)
                    B[8][2*i] = (fem.LEdgesEl[n]*ny_n + fem.LEdgesEl[m]*ny_m)*( Pa[1] - fem.B[1] )/(6*fem.He)                
                    #
                    B[1][2*i+1] = (fem.LEdgesEl[n]*ny_n + fem.LEdgesEl[m]*ny_m)/6     
                    B[2][2*i+1] = (fem.LEdgesEl[n]*nx_n + fem.LEdgesEl[m]*nx_m)/6
                    B[4][2*i+1] = (fem.LEdgesEl[n]*ny_n + fem.LEdgesEl[m]*ny_m)*( Pa[0] - fem.B[0] )/(6*fem.He)    
                    B[5][2*i+1] = (fem.LEdgesEl[n]*nx_n + fem.LEdgesEl[m]*nx_m)*( Pa[0] - fem.B[0] )/(6*fem.He)
                    B[7][2*i+1] = (fem.LEdgesEl[n]*ny_n + fem.LEdgesEl[m]*ny_m)*( Pa[1] - fem.B[1] )/(6*fem.He)
                    B[8][2*i+1] = (fem.LEdgesEl[n]*nx_n + fem.LEdgesEl[m]*nx_m)*( Pa[1] - fem.B[1] )/(6*fem.He)
                
                else:
                    z = i-fem.NEdgesEl
                    nx = fem.NormalEl[z][0]
                    ny = fem.NormalEl[z][1]
                    if z==fem.NEdgesEl-1:
                        Pa = fem.VertexEl[z] 
                        Pb = fem.VertexEl[0]
                    else:
                        Pa = fem.VertexEl[z]
                        Pb = fem.VertexEl[z+1]
                                    
                    #
                    B[0][2*i] = fem.LEdgesEl[z]*nx/3
                    B[2][2*i] = fem.LEdgesEl[z]*ny/3
                    B[3][2*i] = fem.LEdgesEl[z]*nx*( Pa[0] + Pb[0] - 2*fem.B[0] )/(3*fem.He)
                    B[5][2*i] = fem.LEdgesEl[z]*ny*( Pa[0] + Pb[0] - 2*fem.B[0] )/(3*fem.He)
                    B[6][2*i] = fem.LEdgesEl[z]*nx*( Pa[1] + Pb[1] - 2*fem.B[1] )/(3*fem.He)
                    B[8][2*i] = fem.LEdgesEl[z]*ny*( Pa[1] + Pb[1] - 2*fem.B[1] )/(3*fem.He)
                    #
                    B[1][2*i+1] = fem.LEdgesEl[z]*ny/3        
                    B[2][2*i+1] = fem.LEdgesEl[z]*nx/3   
                    B[4][2*i+1] = fem.LEdgesEl[z]*ny*( Pa[0] + Pb[0] - 2*fem.B[0] )/(3*fem.He)    
                    B[5][2*i+1] = fem.LEdgesEl[z]*nx*( Pa[0] + Pb[0] - 2*fem.B[0] )/(3*fem.He)  
                    B[7][2*i+1] = fem.LEdgesEl[z]*ny*( Pa[1] + Pb[1] - 2*fem.B[1] )/(3*fem.He)
                    B[8][2*i+1] = fem.LEdgesEl[z]*nx*( Pa[1] + Pb[1] - 2*fem.B[1] )/(3*fem.He)         
                
        
            # Internal node (integral)
            D.append([1, 0, fem.IntMonomials[1]/fem.Area, 0, fem.IntMonomials[2]/fem.Area, 0, 
                      fem.IntMonomials[4]/fem.Area, 0, fem.IntMonomials[3]/fem.Area, 0, 
                      fem.IntMonomials[5]/fem.Area, 0])
            D.append([0, 1, 0, fem.IntMonomials[1]/fem.Area, 0, fem.IntMonomials[2]/fem.Area, 
                      0, fem.IntMonomials[4]/fem.Area, 0, fem.IntMonomials[3]/fem.Area, 
                      0, fem.IntMonomials[5]/fem.Area])
            # B matrix add divN
            B[3][fem.Ndofs-2] = -fem.Area/fem.He
            B[8][fem.Ndofs-2] = -fem.Area/fem.He   
            B[5][fem.Ndofs-1] = -fem.Area/fem.He
            B[7][fem.Ndofs-1] = -fem.Area/fem.He
        
            # Compute Pi operator
            PN = np.dot( invG, B)
            # Compute consistent part of stiffness matrix
            Mst = np.dot(block_diag(C, C, C), G)
            M = np.dot( np.dot(np.transpose(PN),Mst), PN )
            # Compute stabilization part of stiffness matrix
            St = np.dot( np.dot(alf,np.trace(M)), np.add(np.identity(fem.Ndofs) ,
                            -np.dot(np.dot(D, np.linalg.inv(np.dot(np.transpose(D), D))), np.transpose(D))))
            # Stiffness matrix 
            localK = np.zeros(fem.Ndofs)        
            localK = np.add(M, St) 

            # Vectorize Pi nabla star
            PiN_s = np.reshape(PN,(9*fem.Ndofs))
            
        # TODO: Insert here the degree 3
        #elif fem.DegEl==3:
        #        
        
        return localK, PiN_s
    
    @staticmethod
    def MassMatrix2D(fem,analysis):
        
        # number of Dofs
        dim_phy = analysis.dim_phy
        ngdl = dim_phy*fem.NVertexEl
        localM = np.zeros((ngdl,ngdl))
        
        return localM
    
    @staticmethod
    def BodyLoad2D(fem,analysis):
        
        # number of Dofs
        dim_phy = analysis.dim_phy
        ngdl = dim_phy*fem.NVertexEl
        body_load = np.zeros(ngdl)
    
        return body_load
    
    # 3-D Case
    @staticmethod
    def StiffnessMatrix3D(fem,analysis):
        
        # number of Dofs
        dim_phy = analysis.dim_phy
        ngdl = dim_phy*fem.NVertexEl
        localK = np.zeros((ngdl,ngdl))
        
        return localK
    
    @staticmethod
    def MassMatrix3D(fem,analysis):
        
        # number of Dofs
        dim_phy = analysis.dim_phy
        ngdl = dim_phy*fem.NVertexEl
        localM = np.zeros((ngdl,ngdl))
        
        return localM


'''     OLD FUNCTION
        lm = (fem.PoissonEl*fem.YoungEl)/((1+fem.PoissonEl)*(1-2*fem.PoissonEl))
        mu = (fem.YoungEl)/(2*(1+fem.PoissonEl))

        # MATRIX D
        D = np.zeros((fem.Ndofs,6));
                  
        for i in range(fem.NVertexEl):
            # baric-scaled coord of i-vertex
            x = fem.VertexEl[i][0] - fem.B[0]  
            y = fem.VertexEl[i][1] - fem.B[1]
            
            # first col
            D[2*i][0]   = 1.; D[2*i+1][0] = 0.
            # second col
            D[2*i][1]   = 0.; D[2*i+1][1] = 1.
            # third col
            D[2*i][2]   = y;  D[2*i+1][2] = -x
            # fourth col
            D[2*i][3]   = y;  D[2*i+1][3] = x
            # fiveten col 
            D[2*i][4]   = x;  D[2*i+1][4] = 0.
            # sixth col
            D[2*i][5]   = 0.; D[2*i+1][5] = y
            
             
        # Everage of vertex position
        vb = np.zeros(dim_phy)
        for n in range(fem.NVertexEl):
            vb[0] += fem.VertexEl[n][0]/fem.NVertexEl
            vb[1] += fem.VertexEl[n][1]/fem.NVertexEl
            
        # MATRIX B (transpose)             
        B = np.zeros((6,fem.Ndofs))       
           
        for j in range(fem.NVertexEl):
            # j=0
            B[0][2*j]   = 1./fem.NVertexEl
            B[0][2*j+1] = 0.
        
            #j=1
            B[1][2*j]   = 0.;
            B[1][2*j+1] = 1./fem.NVertexEl
        
            # j=2
            B[2][2*j]   = (1./fem.NVertexEl)*(fem.VertexEl[j][1] - vb[1])
            B[2][2*j+1] = (1./fem.NVertexEl)*(-fem.VertexEl[j][0] + vb[0])

            # j=3  
            B[3][2*j]   = mu*( -fem.VertexEl[np.mod(j+1,fem.NVertexEl)][0] 
                            + fem.VertexEl[np.mod(j-1,fem.NVertexEl)][0] )
            B[3][2*j+1] = mu*(  fem.VertexEl[np.mod(j+1,fem.NVertexEl)][1] 
                            - fem.VertexEl[np.mod(j-1,fem.NVertexEl)][1] )

            # j=4
            B[4][2*j]   = (2*mu+lm)/2.*( fem.VertexEl[np.mod(j+1,fem.NVertexEl)][1] 
                            - fem.VertexEl[np.mod(j-1,fem.NVertexEl)][1] );
            B[4][2*j+1] = (lm/2.)*( -fem.VertexEl[np.mod(j+1,fem.NVertexEl)][0] 
                            + fem.VertexEl[np.mod(j-1,fem.NVertexEl)][0] );

            # j=5
            B[5][2*j]   = (lm/2.)*( fem.VertexEl[np.mod(j+1,fem.NVertexEl)][1] 
                            - fem.VertexEl[np.mod(j-1,fem.NVertexEl)][1] ); 
            B[5][2*j+1] = (2*mu+lm)/2.*( -fem.VertexEl[np.mod(j+1,fem.NVertexEl)][0] 
                            + fem.VertexEl[np.mod(j-1,fem.NVertexEl)][0] )

        # Build projectors
        G = np.dot(B, D)           # matrix G from Hitchhikers paper
        invG = np.linalg.inv(G)    # inverse of matrix G         
        PNs  = np.dot(invG, B)     # PiNabla star projector (polynomial basis)
        PN   = np.dot(D, PNs)      # PiNabla projector (Vh basis) 
        
        #  localK
        Gt = np.zeros((6,6))
        Gt[0:3][:] = np.zeros((3,6))   # matrix G tilde (null on the kernel)
        Gt[3:6][:] = G[3:6][:]
                
        # Consistent part
        M  = np.dot( np.dot(np.transpose(PNs),Gt),PNs)  
        # Stabilization
        St = np.dot(alf/2.,np.trace(M))*np.dot(np.transpose( np.add( np.identity(fem.Ndofs),-PN) ), np.add(np.identity(fem.Ndofs),-PN ))
                          
        # Local K          
        localK = np.zeros(fem.Ndofs)
        localK = np.add(M,St)

        # Vectorize Pi nabla star
        PiN_s = np.reshape(PNs,(6*fem.Ndofs))
        
        # alternative option for stabilization
        # localK = M + alf*(trace(M)/2)*(eye(2*m)-D*inv((D')*D)*(D'));
        # KERNEL PROJECTION
        #
        # Let the vertex based scalar product 
        # <v,w> = (1/N) sum_{i=1}^m v(vertex_i)\cdot w(vertex_i) .
        # Then, the three projectors on the kernel are
        # P_1(v) := <[1;0],v>
        # P_2(v) := <[0;1],v>
        # P_3(v) := <[(y-yvb);-(x-xvb)],v> with vb the "vertex baricenter"

'''
       