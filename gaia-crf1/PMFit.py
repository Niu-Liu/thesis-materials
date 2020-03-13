# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 14:15:05 2016

@author: Neo

proper motion fitting.

observed total proper motion = peculiar proper motion of star + 
    effect of peculiar proper motion of Sun (Barycenter of Solar system) + 
    Rotation motion of the Galaxy(rigid part and differential part)
    
However the peculiar proper motions of stars are assumed to be isotropic or random 
and then the corresponding mean is zero. 

As a result, we can obtain(proper motions are expressed in Galactic coordinate just for convenience):

k*(mu_l*, mu_b)^T = [[-sin(l)          cos(l)             0],
                     [-sib(b)*cos(l)] -sin(b)*sin(l) cos(b)]]*(-X -Y -Z)^T/r +
                    [[-sib(b)*cos(l)] -sin(b)*sin(l) cos(b)],
                     [-sin(l)          cos(l)             0]]*( 0  0  B)^T*k   +
                    A*(cos(2l)*cos(b)  -sin(2l)*sib(b)*cos(b))^T*k

NB:                     
k-> kapa = 4.74047
r: heliocentric distance
^T: transpose
A, B: Oort constants.
l*: l*cos(b)

Oct 7: fix some bugs by Niu.
"""

import numpy as np
sin = np.sin
cos = np.cos
k =4.74047

def PMfit(pml, pmb, pmerrl, pmerrb, l, b, r):
    '''
    Actually here 'pml' = pml*cos(b), which is usually labelled as \mu_{l*}
    
    A*x = b
    x = (X Y Z A B)^T, A: 5x5 matrix
    '''
    kpml = k*pml
    kpmb = k*pmb
    
    A11 = np.sum( sin(l)**2/r**2/pmerrl**2 ) + \
          np.sum( sin(b)**2*cos(l)**2/r**2/pmerrb**2 )
    A12 = np.sum(-sin(l)*cos(l)/r**2/pmerrl**2 ) + \
          np.sum( sin(b)**2*cos(l)*sin(l)/r**2/pmerrb**2 )
    A13 = 0.0 + \
          np.sum(-sin(b)*cos(l)*cos(b)/r**2/pmerrb**2 )
    A14 = np.sum( sin(l)/r*cos(2*l)*cos(b)/pmerrl**2 ) + \
          np.sum(-sin(b)*cos(l)/r*sin(2*l)*cos(b)*sin(b)/pmerrb**2 )
    A15 = np.sum( sin(l)/r*cos(b)/pmerrl**2 ) + \
          0.0
        
    A22 = np.sum( cos(l)**2/r**2/pmerrl**2 ) + \
          np.sum( sin(b)**2*sin(l)**2/r**2/pmerrb**2 )
    A23 = 0.0 + \
          np.sum(-sin(b)*sin(l)*cos(b)/r**2/pmerrb**2 )
    A24 = np.sum(-cos(l)/r*cos(2*l)*cos(b)/pmerrl**2 ) + \
          np.sum( sin(b)*sin(l)/r*-sin(2*l)*cos(b)*sin(b)/pmerrb**2 )
    A25 = np.sum(-cos(l)/r*cos(b)/pmerrl**2 ) + \
          0.0
        
    A33 = 0.0 + \
          np.sum( cos(b)**2/r**2/pmerrb**2 )
    A34 = 0.0 + \
          np.sum(-cos(b)/r*-sin(2*l)*cos(b)*sin(b)/pmerrb**2 )
    A35 = 0.0 + \
          0.0
    
    A44 = np.sum( cos(2*l)**2* cos(b)**2/pmerrl**2 ) + \
          np.sum( sin(2*l)**2*cos(b)**2*sin(b)**2/pmerrb**2 )
    A45 = np.sum( cos(2*l)*cos(b)**2/pmerrl**2 ) + \
          0.0
    
    A55 = np.sum( cos(b)**2/pmerrl**2 ) + \
          0.0
        
    b1 = np.sum( sin(l)/r*kpml/pmerrl**2 ) + \
         np.sum( sin(b)*cos(l)/r*kpmb/pmerrb**2 )
    b2 = np.sum( -cos(l)/r*kpml/pmerrl**2 ) + \
         np.sum( sin(b)*sin(l)/r*kpmb/pmerrb**2 )
    b3 = 0.0 + \
         np.sum( -cos(b)/r*kpmb/pmerrb**2 )
    b4 = np.sum( cos(2*l)*cos(b)*kpml/pmerrl**2 ) + \
         np.sum( -sin(2*l)*cos(b)*sin(b)*kpmb/pmerrb**2 )
    b5 = np.sum( cos(b)*kpml/pmerrl**2 ) + \
        0.0
        
    A = np.mat([[A11, A12, A13, A14, A15],
                [A12, A22, A23, A24, A25],
                [A13, A23, A33, A34, A35],
                [A14, A24, A34, A44, A45],
                [A15, A25, A35, A45, A55]])
                
    B = np.array([b1, b2, b3, b4, b5])
    
    x = np.linalg.solve(A, B)
    
    cov = np.linalg.inv(A)  
    sig = np.sqrt(cov.diagonal())
    sig = np.array(sig).reshape(len(x),)
    
    corrcoef = np.array([ cov[i,j]/sig[i]/sig[j] \
                for j in range(len(x)) for i in range(len(x))])
    corrcoef.resize((len(x), len(x)))
                 
    return x, sig, corrcoef

def PMfit1(pml, pmb, pmerrl, pmerrb, l, b, r):
    '''
    Here we consider 9 parameters x = (S1, S2, S3, D_32, D_13, D_21, D+12, D+13, D+32)^T
    D_21, D+12 similar to Oort constant B, A.
    '''
    kpml = k*pml
    kpmb = k*pmb
    
    A11 = np.sum( sin(l)**2/r**2/pmerrl**2 ) + \
          np.sum( cos(l)**2*sin(b)**2/r**2/pmerrb**2 )
    A12 = np.sum(-sin(l)*cos(l)/r**2/pmerrl**2 ) + \
          np.sum( cos(l)*sin(l)*sin(b)**2/r**2/pmerrb**2 )
    A13 = 0.0 + \
          np.sum(-cos(l)*sin(b)*cos(b)/r**2/pmerrb**2 )
    A14 = np.sum(-sin(l)*cos(l)*sin(b)/r/pmerrl**2 ) + \
          np.sum( sin(b)*cos(l)*sin(l)/r/pmerrb**2 )
    A15 = np.sum(-sin(l)**2*sin(b)/r/pmerrl**2 ) + \
          np.sum(-sin(b)*cos(l)**2/r/pmerrb**2 )
    A16 = np.sum( sin(l)*cos(b)/r/pmerrl**2 ) + \
          0.0
    A17 = np.sum( sin(l)*cos(2*l)*cos(b)/r/pmerrl**2 ) + \
          np.sum(-cos(l)*sin(2*l)*cos(b)*sin(b)**2/r/pmerrb**2 )   
    A18 = np.sum(-sin(l)**2*sin(b)/r/pmerrl**2 ) + \
          np.sum( cos(l)**2*sin(b)*cos(2*b)/r/pmerrb**2 )
    A19 = np.sum( sin(l)*cos(l)*sin(b)/r/pmerrl**2 ) + \
          np.sum( cos(l)*sin(l)*sin(b)*cos(2*b)/r/pmerrb**2 )
        
    A22 = np.sum( cos(l)**2/r**2/pmerrl**2 ) + \
          np.sum( sin(b)**2*sin(l)**2/r**2/pmerrb**2 )
    A23 = 0.0 + \
          np.sum(-sin(l)*sin(b)*cos(b)/r**2/pmerrb**2 )
    A24 = np.sum( cos(l)**2*sin(b)/r/pmerrl**2 ) + \
          np.sum( sin(l)**2*sin(b)/r/pmerrb**2 )
## A25 = - A14
    A25 = np.sum( cos(l)*sin(l)*sin(b)/r/pmerrl**2 ) + \
          np.sum(-sin(l)*cos(l)*sin(b)/r/pmerrb**2 )
    A26 = np.sum(-cos(l)*cos(b)/r/pmerrl**2 ) + \
          0.0
    A27 = np.sum(-cos(l)*cos(2*l)*cos(b)/r/pmerrl**2 ) + \
          np.sum(-sin(l)*sin(2*l)*cos(b)*sin(b)**2/r/pmerrb**2 )
## A28 = A19
    A28 = np.sum( cos(l)*sin(l)*sin(b)/r/pmerrl**2 ) + \
          np.sum( sin(l)*cos(l)*sin(b)*cos(2*b)/r/pmerrb**2 )
    A29 = np.sum(-cos(l)**2*sin(b)/r/pmerrl**2 ) + \
          np.sum( sin(l)**2*sin(b)*cos(2*b)/r/pmerrb**2 )
        
    A33 = 0.0 + \
          np.sum( cos(b)**2/r**2/pmerrb**2 )
    A34 = 0.0 + \
          np.sum(-cos(b)*sin(l)/r/pmerrb**2 )
    A35 = 0.0 + \
          np.sum( cos(b)*cos(l)/r/pmerrb**2 )
    A36 = 0.0 + \
          0.0
    A37 = 0.0 + \
          np.sum( sin(2*l)*cos(b)**2*sin(b)/r/pmerrb**2 )
    A38 = 0.0 + \
          np.sum(-cos(l)*cos(b)*cos(2*b)/r/pmerrb**2 )
    A39 = 0.0 + \
          np.sum(-sin(l)*cos(b)*cos(2*b)/r/pmerrb**2 )
        
    A44 = np.sum( sin(b)**2*cos(l)**2/pmerrl**2 ) + \
          np.sum( sin(l)**2/pmerrb**2 )
    A45 = np.sum( sin(b)**2*cos(l)*sin(l)/pmerrl**2 ) + \
          np.sum( -sin(l)*cos(l)/pmerrb**2 )
    A46 = np.sum( -cos(l)*cos(b)*sin(b)/pmerrl**2 ) + \
          0.0
    A47 = np.sum(-cos(l)*cos(2*l)*cos(b)*sin(b)/pmerrl**2 ) + \
          np.sum(-sin(2*l)*sin(l)*cos(b)*sin(b)/pmerrb**2 )  
    A48 = np.sum( sin(b)**2*cos(l)*sin(l)/pmerrl**2 ) + \
          np.sum( sin(l)*cos(l)*cos(2*b)/pmerrb**2 )
    A49 = np.sum(-sin(b)**2*cos(l)**2/pmerrl**2 ) + \
          np.sum( sin(l)**2*cos(2*b)/pmerrb**2 )   
        
    A55 = np.sum( sin(b)**2*sin(l)**2/pmerrl**2 ) + \
          np.sum( cos(l)**2/pmerrb**2 )
    A56 = np.sum(-cos(b)*sin(b)*sin(l)/pmerrl**2 ) + \
          0.0
    A57 = np.sum(-cos(2*l)*sin(l)*cos(b)*sin(b)/pmerrl**2 ) + \
          np.sum( sin(2*l)*cos(l)*cos(b)*sin(b)/pmerrb**2 )  
    A58 = np.sum( sin(b)**2*sin(l)**2/pmerrl**2 ) + \
          np.sum(-cos(l)**2*cos(2*b)/pmerrb**2 )
## A59 = - A48
    A59 = np.sum(-sin(l)*cos(l)*sin(b)**2/pmerrl**2 ) + \
          np.sum(-cos(l)*sin(l)*cos(2*b)/pmerrb**2 )

    A66 = np.sum( cos(b)**2/pmerrl**2 ) + \
          0.0
    A67 = np.sum( cos(2*l)*cos(b)**2/pmerrl**2 ) + \
          0.0
    A68 = np.sum(-cos(b)*sin(b)*sin(l)/pmerrl**2 ) + \
          0.0
    A69 = np.sum( cos(b)*sin(b)*cos(l)/pmerrl**2 ) + \
          0.0
        
    A77 = np.sum( cos(2*l)**2*cos(b)**2/pmerrl**2 ) + \
          np.sum( sin(2*l)**2*sin(b)**2*cos(b)**2/pmerrb**2 )
    A78 = np.sum(-cos(2*l)*sin(l)*cos(b)*sin(b)/pmerrl**2 ) + \
          np.sum(-sin(2*l)*cos(l)*sin(b)*cos(b)*cos(2*b)/pmerrb**2 )
    A79 = np.sum( cos(2*l)*cos(l)*cos(b)*sin(b)/pmerrl**2 ) + \
          np.sum(-sin(2*l)*sin(l)*sin(b)*cos(b)*cos(2*b)/pmerrb**2 )

    A88 = np.sum( sin(b)**2*sin(l)**2/pmerrl**2 ) + \
          np.sum( cos(l)**2*cos(2*b)**2/pmerrb**2 )
    A89 = np.sum(-sin(l)*cos(l)*sin(b)**2/pmerrl**2 ) + \
          np.sum( cos(l)*sin(l)*cos(2*b)**2/pmerrb**2 )
    
    A99 = np.sum( sin(b)**2*cos(l)**2/pmerrl**2 ) + \
          np.sum( sin(l)**2*cos(2*b)**2/pmerrb**2 )
        
    b1 = np.sum( sin(l)/r*kpml/pmerrl**2 ) + \
         np.sum( sin(b)*cos(l)/r*kpmb/pmerrb**2 )
    b2 = np.sum( -cos(l)/r*kpml/pmerrl**2 ) + \
         np.sum( sin(b)*sin(l)/r*kpmb/pmerrb**2 )
    b3 = 0.0 + \
         np.sum( -cos(b)/r*kpmb/pmerrb**2 )
    b4 = np.sum( -sin(b)*cos(l)*kpml/pmerrl**2 ) + \
         np.sum( sin(l)*kpmb/pmerrb**2 )
    b5 = np.sum( -sin(b)*sin(l)*kpml/pmerrl**2 ) + \
         np.sum( -cos(l)*kpmb/pmerrb**2 )
    b6 = np.sum( cos(b)*kpml/pmerrl**2 ) + \
         0.0
    b7 = np.sum( cos(2*l)*cos(b)*kpml/pmerrl**2 ) + \
         np.sum( -sin(2*l)*cos(b)*sin(b)*kpmb/pmerrb**2 )
    b8 = np.sum( -sin(b)*sin(l)*kpml/pmerrl**2 ) + \
         np.sum( cos(l)*cos(2*b)*kpmb/pmerrb**2 )
    b9 = np.sum( cos(l)*sin(b)*kpml/pmerrl**2 ) + \
         np.sum( sin(l)*cos(2*b)*kpmb/pmerrb**2 )
        
    A = np.mat([[A11, A12, A13, A14, A15, A16, A17, A18, A19],
                [A12, A22, A23, A24, A25, A26, A27, A28, A29],
                [A13, A23, A33, A34, A35, A36, A37, A38, A39],
                [A14, A24, A34, A44, A45, A46, A47, A48, A49],
                [A15, A25, A35, A45, A55, A56, A57, A58, A59],
                [A16, A26, A36, A46, A56, A66, A67, A68, A69],
                [A17, A27, A37, A47, A57, A67, A77, A78, A79],
                [A18, A28, A38, A48, A58, A68, A78, A88, A89],
                [A19, A29, A39, A49, A59, A69, A79, A89, A99]])
                
    B = np.array([b1, b2, b3, b4, b5, b6, b7, b8, b9])
    
    x = np.linalg.solve(A, B)
    
    cov = np.linalg.inv(A)  
    sig = np.sqrt(cov.diagonal())
    sig = np.array(sig).reshape(len(x),)
    
    corrcoef = np.array([ cov[i,j]/sig[i]/sig[j] \
                for j in range(len(x)) for i in range(len(x))])
    corrcoef.resize((len(x), len(x)))
                 
    return x, sig, corrcoef

def PMfit2(pml, pmb, pmerrl, pmerrb, l, b, r):
    '''
    Here we consider 7 parameters X = (S1, S2, S3, D_32, D_13, D_21, D+12, D+13, D+32)^T
    D_21, D+12 similar to Oort constant B, A.
    This model is similar to 9-parameter model, but we set D_32, D+32 to zero.
    '''
    kpml = k*pml
    kpmb = k*pmb
    
    A11 = np.sum( sin(l)**2/r**2/pmerrl**2 ) + \
          np.sum( cos(l)**2*sin(b)**2/r**2/pmerrb**2 )
    A12 = np.sum(-sin(l)*cos(l)/r**2/pmerrl**2 ) + \
          np.sum( cos(l)*sin(l)*sin(b)**2/r**2/pmerrb**2 )
    A13 = 0.0 + \
          np.sum(-cos(l)*sin(b)*cos(b)/r**2/pmerrb**2 )

    A15 = np.sum(-sin(l)**2*sin(b)/r/pmerrl**2 ) + \
          np.sum(-sin(b)*cos(l)**2/r/pmerrb**2 )
    A16 = np.sum( sin(l)*cos(b)/r/pmerrl**2 ) + \
          0.0
    A17 = np.sum( sin(l)*cos(2*l)*cos(b)/r/pmerrl**2 ) + \
          np.sum(-cos(l)*sin(2*l)*cos(b)*sin(b)**2/r/pmerrb**2 )   
    A18 = np.sum(-sin(l)**2*sin(b)/r/pmerrl**2 ) + \
          np.sum( cos(l)**2*sin(b)*cos(2*b)/r/pmerrb**2 )

        
    A22 = np.sum( cos(l)**2/r**2/pmerrl**2 ) + \
          np.sum( sin(b)**2*sin(l)**2/r**2/pmerrb**2 )
    A23 = 0.0 + \
          np.sum(-sin(l)*sin(b)*cos(b)/r**2/pmerrb**2 )

## A25 = - A14
    A25 = np.sum( cos(l)*sin(l)*sin(b)/r/pmerrl**2 ) + \
          np.sum(-sin(l)*cos(l)*sin(b)/r/pmerrb**2 )
    A26 = np.sum(-cos(l)*cos(b)/r/pmerrl**2 ) + \
          0.0
    A27 = np.sum(-cos(l)*cos(2*l)*cos(b)/r/pmerrl**2 ) + \
          np.sum(-sin(l)*sin(2*l)*cos(b)*sin(b)**2/r/pmerrb**2 )
## A28 = A19
    A28 = np.sum( cos(l)*sin(l)*sin(b)/r/pmerrl**2 ) + \
          np.sum( sin(l)*cos(l)*sin(b)*cos(2*b)/r/pmerrb**2 )

        
    A33 = 0.0 + \
          np.sum( cos(b)**2/r**2/pmerrb**2 )

    A35 = 0.0 + \
          np.sum( cos(b)*cos(l)/r/pmerrb**2 )
    A36 = 0.0 + \
          0.0
    A37 = 0.0 + \
          np.sum( sin(2*l)*cos(b)**2*sin(b)/r/pmerrb**2 )
    A38 = 0.0 + \
          np.sum(-cos(l)*cos(b)*cos(2*b)/r/pmerrb**2 ) 
        
    A55 = np.sum( sin(b)**2*sin(l)**2/pmerrl**2 ) + \
          np.sum( cos(l)**2/pmerrb**2 )
    A56 = np.sum(-cos(b)*sin(b)*sin(l)/pmerrl**2 ) + \
          0.0
    A57 = np.sum(-cos(2*l)*sin(l)*cos(b)*sin(b)/pmerrl**2 ) + \
          np.sum( sin(2*l)*cos(l)*cos(b)*sin(b)/pmerrb**2 )  
    A58 = np.sum( sin(b)**2*sin(l)**2/pmerrl**2 ) + \
          np.sum(-cos(l)**2*cos(2*b)/pmerrb**2 )

    A66 = np.sum( cos(b)**2/pmerrl**2 ) + \
          0.0
    A67 = np.sum( cos(2*l)*cos(b)**2/pmerrl**2 ) + \
          0.0
    A68 = np.sum(-cos(b)*sin(b)*sin(l)/pmerrl**2 ) + \
          0.0
        
    A77 = np.sum( cos(2*l)**2*cos(b)**2/pmerrl**2 ) + \
          np.sum( sin(2*l)**2*sin(b)**2*cos(b)**2/pmerrb**2 )
    A78 = np.sum(-cos(2*l)*sin(l)*cos(b)*sin(b)/pmerrl**2 ) + \
          np.sum(-sin(2*l)*cos(l)*sin(b)*cos(b)*cos(2*b)/pmerrb**2 )


    A88 = np.sum( sin(b)**2*sin(l)**2/pmerrl**2 ) + \
          np.sum( cos(l)**2*cos(2*b)**2/pmerrb**2 )


        
    b1 = np.sum( sin(l)/r*kpml/pmerrl**2 ) + \
         np.sum( sin(b)*cos(l)/r*kpmb/pmerrb**2 )
    b2 = np.sum( -cos(l)/r*kpml/pmerrl**2 ) + \
         np.sum( sin(b)*sin(l)/r*kpmb/pmerrb**2 )
    b3 = 0.0 + \
         np.sum( -cos(b)/r*kpmb/pmerrb**2 )

    b5 = np.sum( -sin(b)*sin(l)*kpml/pmerrl**2 ) + \
         np.sum( -cos(l)*kpmb/pmerrb**2 )
    b6 = np.sum( cos(b)*kpml/pmerrl**2 ) + \
         0.0
    b7 = np.sum( cos(2*l)*cos(b)*kpml/pmerrl**2 ) + \
         np.sum( -sin(2*l)*cos(b)*sin(b)*kpmb/pmerrb**2 )
    b8 = np.sum( -sin(b)*sin(l)*kpml/pmerrl**2 ) + \
         np.sum( cos(l)*cos(2*b)*kpmb/pmerrb**2 )

        
    A = np.mat([[A11, A12, A13, A15, A16, A17, A18],
                [A12, A22, A23, A25, A26, A27, A28],
                [A13, A23, A33, A35, A36, A37, A38],
                [A15, A25, A35, A55, A56, A57, A58],
                [A16, A26, A36, A56, A66, A67, A68],
                [A17, A27, A37, A57, A67, A77, A78],
                [A18, A28, A38, A58, A68, A78, A88]])
                
    B = np.array([b1, b2, b3, b5, b6, b7, b8])
    
    x = np.linalg.solve(A, B)
    
    cov = np.linalg.inv(A)  
    sig = np.sqrt(cov.diagonal())
    sig = np.array(sig).reshape(len(x),)
    
    corrcoef = np.array([ cov[i,j]/sig[i]/sig[j] \
                for j in range(len(x)) for i in range(len(x))])
    corrcoef.resize((len(x), len(x)))
                 
    return x, sig, corrcoef