""" quadrilateral element four node element
 
    Author: Dr. Marco Pingaro
"""

import numpy as np
#from sys import exit

class Functions:

    @staticmethod
    def gauss_quadrature(npt_gauss):
        """ Gauss-Quadrature four points """
    
        # Points of quadrature 4 points
        if npt_gauss==4:
            pt = np.zeros((4,2))
            pt[0][0] = -1.0/np.sqrt(3)
            pt[0][1] = -1.0/np.sqrt(3)
            pt[1][0] =  1.0/np.sqrt(3)
            pt[1][1] = -1.0/np.sqrt(3)
            pt[2][0] =  1.0/np.sqrt(3)
            pt[2][1] =  1.0/np.sqrt(3)
            pt[3][0] = -1.0/np.sqrt(3)
            pt[3][1] =  1.0/np.sqrt(3)

            # Weigh of quadrature
            wg = np.zeros(4)
            wg[0] = 1.0
            wg[1] = 1.0
            wg[2] = 1.0
            wg[3] = 1.0
            
        return pt,wg

    @staticmethod
    def eval_grd_bilinear_functions(npt_gauss):
        """ Evaluetion of the gradient fo the shape functions in the parametric domain """

        pt,wg = Functions.gauss_quadrature(npt_gauss) # call gauss quadrature
        grd_phi = np.zeros((4,2,npt_gauss))
    
        # Cycle on gauss points
        for ipt in range(npt_gauss):
            grd_phi[0][0][ipt] = -(1-pt[ipt][1])/4.
            grd_phi[0][1][ipt] = -(1-pt[ipt][0])/4.
    
            grd_phi[1][0][ipt] =  (1-pt[ipt][1])/4.
            grd_phi[1][1][ipt] = -(1+pt[ipt][0])/4.

            grd_phi[2][0][ipt] =  (1+pt[ipt][1])/4.
            grd_phi[2][1][ipt] =  (1+pt[ipt][0])/4.
            
            grd_phi[3][0][ipt] = -(1+pt[ipt][1])/4.
            grd_phi[3][1][ipt] =  (1-pt[ipt][0])/4.
            # end cycle on gauss points

        return grd_phi
        
    @staticmethod
    def eval_bilinear_functions(npt_gauss):
        """ Evaluetion of the shape functions in the parametric domain """

        pt,wg = Functions.gauss_quadrature(npt_gauss) # call gauss quadrature
        fun_phi = np.zeros((4,npt_gauss))
    
        # Cycle on gauss points
        for ipt in range(npt_gauss):
            fun_phi[0][ipt] = (1-pt[ipt][0])*(1-pt[ipt][1])/4.
            fun_phi[1][ipt] = (1+pt[ipt][0])*(1-pt[ipt][1])/4.
            fun_phi[2][ipt] = (1+pt[ipt][0])*(1+pt[ipt][1])/4.
            fun_phi[3][ipt] = (1-pt[ipt][0])*(1+pt[ipt][1])/4.
            # end cycle on gauss points

        return fun_phi    
    
    @staticmethod
    def jacobian(npt_gauss):
        """Compute the Jacobian matrix"""

        pt,wg = Functions.gauss_quadrature(npt_gauss) # call gauss quadrature
        J = np.zeros((2,2,npt_gauss))
        # call the gradient of shape functions
        grd_phi = Functions.eval_grd_bilinear_functions(npt_gauss)
    
        # Cycle on gauss points
        for ipt in range(npt_gauss):
            # Cycle on shape functions
            for ifun in range(4):
                J[0][0][ipt] += grd_phi[ifun][0][ipt]*pt[ifun][0]
                J[0][1][ipt] += grd_phi[ifun][0][ipt]*pt[ifun][1]
                J[1][0][ipt] += grd_phi[ifun][1][ipt]*pt[ifun][0]
                J[1][1][ipt] += grd_phi[ifun][1][ipt]*pt[ifun][1]
                # end cycle on shape functions
            # end cycle on gauss points 

        return J
    
    @staticmethod
    def det_jacobian(J):
        """ Compute the determinant of the Jacobian matrix """
        
        npt = len(J)
        dJ = np.zeros(npt)    

        for ipt in range(npt):
            dJ[ipt] = J[0][0][ipt]*J[1][1][ipt] - J[0][1][ipt]*J[1][0][ipt]

        return dJ
    
    @staticmethod
    def inv_jacobian(npt_gauss):
        """ Compute the inverse transpose of the jacobian matrix """
    
        J  = Functions.jacobian(npt_gauss)
        dJ = Functions.det_jacobian(J)
    
        invJ = np.zeros((2,2,npt_gauss))

        for ipt in range(npt_gauss):
            invJ[0][0][ipt] =  J[1][1][ipt]/dJ[ipt]
            invJ[0][1][ipt] = -J[1][0][ipt]/dJ[ipt]
            invJ[1][0][ipt] = -J[0][1][ipt]/dJ[ipt]
            invJ[1][1][ipt] =  J[0][0][ipt]/dJ[ipt]

        return invJ

    @staticmethod
    def eval_phi_grd_bilinear_function(pt):
        """ Compute the physical gradient of the shape functions """
    
        grd_phi = Functions.eval_grd_bilinear_functions()
        invJ    = Functions.inv_jacobian(pt)

        grd = np.zeros((4,2,4))
        # Cycle on gauss points
        for ipt in range(0,4):
            # Cycle on shape functions
            for ifun in range(0,4):
                grd[ifun][0][ipt] = invJ[0][0][ipt]*grd_phi[ifun][0][ipt] + \
                    invJ[0][1][ipt]*grd_phi[ifun][1][ipt]
                grd[ifun][1][ipt] = invJ[1][0][ipt]*grd_phi[ifun][0][ipt] + \
                    invJ[1][1][ipt]*grd_phi[ifun][1][ipt]
                # end cycle on shape functions
            # end cycle on gauss points
        return grd
    
    @staticmethod
    def eval_stress_functions(npt_gauss,TypeEl):
        """ Evaluation of the stresses shape functions in the physical domain """

        pt,wg = Functions.gauss_quadrature(npt_gauss) # call gauss quadrature
        J    = Functions.jacobian(npt_gauss)          # compute Jacobian Matrix
        invJ = Functions.inv_jacobian(npt_gauss)      # compute inverse of Jacobian
        detJ = Functions.det_jacobian(J)              # compute determinant of Jacobian
                
        if TypeEl=='PEERSQ2B': # PEERS QUAD WITH 2 BUBBLES
            stress_phi = np.zeros((6,2,npt_gauss))
            # Cycle on gauss points
            for ipt in range(npt_gauss):
                # RT0 - 1
                stress_phi[0][0][ipt] = J[0][1][ipt]*(-0.5+0.5*pt[ipt][1])/detJ[ipt]
                stress_phi[0][1][ipt] = J[1][1][ipt]*(-0.5+0.5*pt[ipt][1])/detJ[ipt]
                # RT0 - 2
                stress_phi[1][0][ipt] = J[0][0][ipt]*(0.5+0.5*pt[ipt][0])/detJ[ipt]
                stress_phi[1][1][ipt] = J[1][0][ipt]*(0.5+0.5*pt[ipt][0])/detJ[ipt]
                # RT0 - 3
                stress_phi[2][0][ipt] = J[0][1][ipt]*(0.5+0.5*pt[ipt][1])/detJ[ipt]
                stress_phi[2][1][ipt] = J[1][1][ipt]*(0.5+0.5*pt[ipt][1])/detJ[ipt]
                # RT0 - 4
                stress_phi[3][0][ipt] = J[0][0][ipt]*(-0.5+0.5*pt[ipt][0])/detJ[ipt]
                stress_phi[3][1][ipt] = J[1][0][ipt]*(-0.5+0.5*pt[ipt][0])/detJ[ipt]
                # BUBBLE - 1
                stress_phi[4][0][ipt] = invJ[0][0][ipt]*(2*pt[ipt][0]*(pt[ipt][1]**2-1))+\
                    invJ[0][1][ipt]*(2*pt[ipt][1]*(pt[ipt][0]**2-1))
                stress_phi[4][1][ipt] = invJ[1][0][ipt]*(2*pt[ipt][0]*(pt[ipt][1]**2-1))+\
                    invJ[1][1][ipt]*(2*pt[ipt][1]*(pt[ipt][0]**2-1))
                # BUBBLE - 2 11
                stress_phi[5][0][ipt] = invJ[0][0][ipt]*((1-3*pt[ipt][0]**2-2*pt[ipt][0]*pt[ipt][1])*(1-pt[ipt][1]**2))+\
                    invJ[0][1][ipt]*((1-3*pt[ipt][1]**2-2*pt[ipt][0]*pt[ipt][1])*(1-pt[ipt][0]**2))
                stress_phi[5][1][ipt] = invJ[1][0][ipt]*((1-3*pt[ipt][0]**2-2*pt[ipt][0]*pt[ipt][1])*(1-pt[ipt][1]**2))+\
                    invJ[1][1][ipt]*((1-3*pt[ipt][1]**2-2*pt[ipt][0]*pt[ipt][1])*(1-pt[ipt][0]**2))
               
            # end cycle on gauss points

        return stress_phi    