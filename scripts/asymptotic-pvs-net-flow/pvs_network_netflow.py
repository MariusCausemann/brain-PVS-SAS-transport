import numpy as np
import numpy.linalg

"""
Data structures for the network and parameters:

*Convention*: Node n is numbered equal to its mother edge/element e.

* indices (list): list of three-tuples representing the network via
  its bifurcations. For each bifurcation/junction, give the mother
  (e0) and daughter edge indices (e1, e2) as (e0, e1, e2).

* paths (list): list of root-to-leaf-paths, one tuple (n0, n1, n2,
  ..., nl) for each leaf node where n0 is the root node, nl is the
  leaf node, and n1, n2, ... are the bifurcation nodes on the path
  between n0 and nl.

The remaining parameters can be given in any dimensions, but sample
units are given here for illustration purposes:

* f (float):             Wave frequency (Hz)
* omega (float):         Angular frequency (Hz)
* lmbda (float):         Wave length (mm)
* varepsilon (float):    Relative wave amplitude (AU)
* r_o (list of floats):  Inner radii (mm) for each element/edge
* r_e (list of floats):  Outer radii (mm) for each element/edge
* L (list of floats):    Length (mm) of each element/edge

Idealized networks are included in the accompanying test script.
"""

def debug(x):
    print(x)

# Bunch of helper functions for ensuring close match between code and
# paper presentation.
def _invR(beta):
    "Helper function for evaluating \mathcal{R}^{-1}(beta), eq. (25)."
    val = np.pi/8*(beta**4 - 1 - (beta**2 - 1)**2/np.log(beta))
    return val

def _R(beta):
    "Helper function for evaluating \mathcal{R}(beta)."
    return 1.0/_invR(beta)

def _delta(beta):
    "Helper function for evaluating Delta(beta), eq. (26)"
    d_dividend = (2 - (beta**2-1)/np.log(beta))**2
    d_divisor = beta**4 - 1 - (beta**2 - 1)**2/np.log(beta)
    delt = d_dividend/d_divisor
    return delt
    
def _beta(r_e, r_o):
    "beta is the ratio of outer-to-inner radius (r_e/r_o) for an annular PVS cross-section."
    return r_e/r_o

def _xi(l):
    "Helper function for expressions in eq. (38)."
    z = 1j
    xi1 = (np.exp(z*l) - 1)/l
    return xi1

def _alpha(l, P, R):
    "Helper function for expressions in eq. (39)."
    z = 1j
    A1 = (0.5 - (1 - np.cos(l))/(l**2))
    A2 = P*(1 - np.exp(z*l))/(2*l**2*R)
    return (A1 + A2.real)

def avg_Q_1_n(P, dP, ell, R, Delta):  
    "Evaluate <Q_1_n>."
    z = 1j
    val = (- dP/(R*ell) + Delta*(1./2 - (1 - np.cos(ell))/(ell**2)) 
           + Delta/(2*ell**2*R)*(P*(1 - np.exp(z*ell))).real)
    return val

def solve_for_P(indices, paths, R, ell, gamma):
    """Main function #1: Solving for the P_i, eq. (38)."""
    
    # Recall: A bifurcating tree has N junctions with n = 2N + 1
    # edges. The number of junctions N plus number of downstream ends
    # (N+1) is also 2N + 1.

    # Form the n x n system of linear equations for determining P: A P
    # = b (complex):
    n = len(R)
    A = np.zeros((n, n), dtype=complex)
    b = np.zeros(n, dtype=complex)
    
    # The complex number i
    z = 1j

    for (i, j, k) in indices:
        I = i # This junction (the index of the mother edge by convention)

        # Set the correct matrix columns for this junction constraint, eq. (38)
        A[I, i] = np.exp(z*ell[i])*gamma[i]/(R[i]*ell[i])
        A[I, j] = - gamma[j]/(R[j]*ell[j])
        A[I, k] = - gamma[k]/(R[k]*ell[k])

        # Set the correct vector row for this junction constraint
        b[I] = gamma[i]*_xi(ell[i]) - gamma[j]*_xi(-ell[j]) \
            - gamma[k]*_xi(-ell[k]) + z*(- gamma[i] + gamma[j] + gamma[k])

    # Apply the additional constraints given in terms of the
    # root-to-leaf paths
    I = len(indices)
    for (k, path) in enumerate(paths):
        x_n = 0.0
        for n in path:
            A[I+k, n] = np.exp(-z*x_n)/gamma[n] # - compared to eq. (40)
            x_n += ell[n]

    # Solve the linear systems for real and imaginary parts of P:
    P = np.linalg.solve(A, b)
    return P

def solve_for_dP(R, ell, Delta, gamma, indices, paths, P):
    "Main function #2: Solving for the dP_i."

    # n x n system of linear (real) equations for determining dP: A dP = b 
    n = len(P)
    A = np.zeros((n, n))
    b = np.zeros(n)
    
    # Define the real linear system A P = b 
    for (i, j, k) in indices:
        I = i # This junction (the index of the mother edge by convention)
        
        # Set right matrix columns for this junction constraint
        A[I, i] = gamma[i]/(R[i]*ell[i])
        A[I, j] = - gamma[j]/(R[j]*ell[j])
        A[I, k] = - gamma[k]/(R[k]*ell[k])

        # Set right vector row for this junction constraint
        b[I] = (gamma[i]*Delta[i]*_alpha(ell[i], P[i], R[i]) 
                - gamma[j]*Delta[j]*_alpha(ell[j], P[j], R[j])
                - gamma[k]*Delta[k]*_alpha(ell[k], P[k], R[k]))

    # Apply additional constraints given in terms of the root-to-leaf paths
    I = len(indices)
    for (k, path) in enumerate(paths):
        for n in path:
            A[I+k, n] = 1.0/gamma[n]
            
    # Solve the linear systems for dP:
    dP = np.linalg.solve(A, b)
    
    return dP

def solve_bifurcating_tree(network):

    # Extract individual parameters from argument list (for readability)
    (indices, paths, r_o, r_e, L, k, omega, varepsilon) = network

    # Compute scaled parameters for solver functions
    beta = [_beta(re, ro) for (re, ro) in zip(r_e, r_o)]
    ell = [k*l for l in L]

    # Computation of dimension-less parameters
    R = [_R(b) for b in beta]
    Delta = [_delta(b) for b in beta]
    gamma = [(r_o_n/r_o[0])**2 for r_o_n in r_o]
    
    debug("Solving for P")
    P = solve_for_P(indices, paths, R, ell, gamma)
    
    debug("Solving for dP")
    dP = solve_for_dP(R, ell, Delta, gamma, indices, paths, P)

    debug("Evaluating the flow rates Q1[i]: ...")
    Q1 = np.zeros(len(dP))
    for (i, _) in enumerate(dP):
        Q1[i] = avg_Q_1_n(P[i], dP[i], ell[i], R[i], Delta[i])

    return (P, dP, Q1)

def Qprime(Q, varepsilon, omega, L, k, Delta, r_o):
    scale = 2*np.pi*varepsilon*omega*r_o**2/k
    val = scale*Q
    return val

