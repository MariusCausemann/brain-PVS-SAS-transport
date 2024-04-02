import numpy as np
import numpy.linalg
from numpy import pi

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

Idealized networks included:
* single_bifurcation_data
* three_junction_data
* murray_tree_data

"""

__DEBUG_MODE=False
def debug(x):
    if __DEBUG_MODE:
        print(x)

# Bunch of helper functions for ensuring close match between code and
# paper presentation.
def _invR(beta):
    "Helper function for evaluating \mathcal{R}^{-1}(beta)."
    val = np.pi/8*(beta**4 - 1 - (beta**2 - 1)**2/np.log(beta))
    return val

def _R(beta):
    "Helper function for evaluating \mathcal{R}(beta)."
    return 1.0/_invR(beta)

def _delta(beta):
    "Helper function for evaluating Delta(beta)."
    delt = ((2 - (beta**2-1)/np.log(beta))**2)/(beta**4 - 1 - (beta**2 - 1)**2/np.log(beta))
    return delt
    
def _beta(r_e, r_o):
    return r_e/r_o

def _alpha(l, P, R):
    "Helper function for matrix/vector expression."
    z = 1j
    A1 = (0.5 - (1 - np.cos(l))/(l**2))
    A2 = P*(1 - np.exp(z*l))/(2*l**2*R)
    return (A1 + A2.real)

def _xi(l):
    "Helper function for matrix/vector expression."
    z = 1j
    xi1 = (np.exp(z*l) - 1)/l
    return xi1

def avg_Q_1_n(P, dP, ell, R, Delta):  
    "Evaluate <Q_1_n>."
    z = 1j
    val = (- dP/(R*ell) + Delta*(1./2 - (1 - np.cos(ell))/(ell**2)) 
           + Delta/(2*ell**2*R)*(P*(1 - np.exp(z*ell))).real)
    return val


def solve_for_P(indices, paths, R, ell):
    """Main function #1: Solving for the P_i."""
    
    # Recall: A bifurcating tree has N junctions with n = 2N + 1
    # edges. The number of junctions N plus number of downstream ends
    # (N+1) is also 2N + 1.

    # Form the n x n system of linear equations for determining P: A P = b (complex):
    n = len(R)
    A = np.zeros((n, n), dtype=complex)
    b = np.zeros(n, dtype=complex)
    
    # The complex number i
    z = 1j

    for (i, j, k) in indices:
        I = i # This junction (the index of the mother edge by convention)

        # Set the correct matrix columns for this junction constraint
        A[I, i] = np.exp(z*ell[i])/(R[i]*ell[i])
        A[I, j] = - 1.0/(R[j]*ell[j])
        A[I, k] = - 1.0/(R[k]*ell[k])

        # Set the correct vector row for this junction constraint
        b[I] = _xi(ell[i]) + z - _xi(-ell[j]) - _xi(-ell[k])

    # Apply the additional constraints given in terms of the root-to-leaf paths
    I = len(indices)
    for (k, path) in enumerate(paths):
        x_n = 0.0
        for n in path:
            A[I+k, n] = np.exp(-z*x_n)
            x_n += ell[n]

    # Solve the linear systems for real and imaginary parts of P:
    P = np.linalg.solve(A, b)
    return P

def solve_for_dP(R, ell, Delta, indices, paths, P):
    "Main function #2: Solving for the dP_i."

    # n x n system of linear (real) equations for determining dP: A dP = b 
    n = len(P)
    A = np.zeros((n, n))
    b = np.zeros(n)
    
    # Define the real linear system A P = b 
    for (i, j, k) in indices:
        I = i # This junction (the index of the mother edge by convention)
        
        # Set right matrix columns for this junction constraint
        A[I, i] = 1.0/(R[i]*ell[i])
        A[I, j] = - 1.0/(R[j]*ell[j])
        A[I, k] = - 1.0/(R[k]*ell[k])

        # Set right vector row for this junction constraint
        b[I] = (Delta[i]*_alpha(ell[i], P[i], R[i]) 
                - Delta[j]*_alpha(ell[j], P[j], R[j])
                - Delta[k]*_alpha(ell[k], P[k], R[k]))

    # Apply additional constraints given in terms of the root-to-leaf paths
    I = len(indices)
    for (k, path) in enumerate(paths):
        for n in path:
            A[I+k, n] = 1.0
            
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
    
    debug("Solving for P")
    P = solve_for_P(indices, paths, R, ell)
    
    debug("Solving for dP")
    dP = solve_for_dP(R, ell, Delta, indices, paths, P)

    debug("Evaluating the flow rates Q1[i]: ...")
    Q1 = np.zeros(len(dP))
    for (i, _) in enumerate(dP):
        Q1[i] = avg_Q_1_n(P[i], dP[i], ell[i], R[i], Delta[i])

    return (P, dP, Q1)

def solve_three_junction(network):

    # Extract individual parameters from argument list (for readability)
    (indices, paths, r_o, r_e, L, k, omega, varepsilon) = network

    # Compute scaled parameters for solver functions
    beta = [_beta(re, ro) for (re, ro) in zip(r_e, r_o)]
    ell = [k*l for l in L]

    # Computation of dimension-less parameters
    R = [_R(b) for b in beta]
    Delta = [_delta(b) for b in beta]
    
    n = 3
    A = np.zeros((n, n), dtype=complex)
    b = np.zeros(n, dtype=complex)
    
    # The complex number i
    z = 1j

    # Define the complex linear system A P = b in terms of the real
    # and imaginary part: Ar Pr = Br and Ai Pi = Bi:

    # Set right matrix columns for this junction constraint
    # [P0, P1, P3]
    A[0, :] = [np.exp(z*ell[0])/(R[0]*ell[0]), - 2*1./(R[1]*ell[1]), 0]
    A[1, :] = [0, np.exp(z*ell[1])/(R[1]*ell[1]), - 2*1./(R[3]*ell[3])]
    A[2, :] = [1, np.exp(-z*ell[0]), np.exp(-z*ell[0] - z*ell[1])]

    # Set right vector row for this junction constraint
    b[0] = _xi(ell[0]) + z - 2*_xi(-ell[1])  
    b[1] = _xi(ell[1]) + z - 2*_xi(-ell[3])  
    b[2] = 0

    # Solve the linear systems for real and imaginary parts of P:
    debug("Solving A P = b")
    P = np.linalg.solve(A, b)
    (P0, P1, P3) = P

    P = (P0, P1, P1, P3, P3, P3, P3)

    debug("Solving for dP")
    dP = solve_for_dP(R, ell, Delta, indices, paths, P)

    debug("Evaluating the flow rates Q1[i]: ...")
    Q1 = np.zeros(len(dP))
    for (i, _) in enumerate(dP):
        Q1[i] = avg_Q_1_n(P[i], dP[i], ell[i], R[i], Delta[i])
    
    return (P, dP, Q1)

# =================================================================================
# Helper functions for setting up idealized networks:
# 
# * single_bifurcation_data
# * three_junction_data
# * murray_tree_data
#
# =================================================================================
def single_bifurcation_data():
    """Data for a single bifurcation test case. Data for network
    consisting of three edges and a single bifurcation, with inner and
    outer radius and lengths of the daughter vessels half those of the
    mother vessel.

    """
    indices = [(0, 1, 2),]
    paths = [(0, 1), (0, 2)]
    
    # Peristaltic wave parameters: wave length lmbda and (angular) wave number k
    f = 1.0                 # frequency (Hz = 1/s)
    omega = 2*np.pi*f       # Angular frequency (Hz)
    lmbda = 2.0             # mm    
    k = 2*pi/lmbda          # wave number (1/mm)
    varepsilon = 0.1        # AU 
    ro0 = 0.1               # Base inner radius (mm)
    re0 = 0.2               # Base outer radius (mm)
        
    r_o = [ro0, ro0/2, ro0/2]  # Inner radii (mm) for each element/edge
    r_e = [re0, re0/2, re0/2]  # Outer radii (mm) for each element/edge
    L = [1.0, 0.5, 0.5]    # Element lengths (mm)
        
    data =  (indices, paths, r_o, r_e, L, k, omega, varepsilon)
    return data
    
def three_junction_data():
    """Data for a three bifurcation test case. Data for network consisting
    of seven edges and three bifurcation. Inner and outer radius of
    the daughter vessels are half that of the mother vessel.
    """
    
    indices = [(0, 1, 2), (1, 3, 4), (2, 5, 6)]
    paths = [(0, 1, 3), (0, 1, 4), (0, 2, 5), (0, 2, 6)]
    
    f = 1.0                 # frequency (Hz = 1/s)
    omega = 2*np.pi*f       # Angular frequency (Hz)
    lmbda = 2.0             # mm    
    k = 2*np.pi/lmbda       # wave number (1/mm)
    varepsilon = 0.1        # AU 
    ro0 = 0.1               # Base inner radius (mm)
    re0 = 0.2               # Base outer radius (mm)
    ro1 = ro0/2; ro2 = ro1/2
    re1 = re0/2; re2 = re1/2
        
    r_o = [ro0, ro1, ro1, ro2, ro2, ro2, ro2]  # Inner radii (mm) for each element/edge
    r_e = [re0, re1, re1, re2, re2, re2, re2]  # Outer radii (mm) for each element/edge
    L = [1.0, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25]    # Element lengths (mm)
        
    data = (indices, paths, r_o, r_e, L, k, omega, varepsilon)
    return data
    
def murray_tree_data(m=1, r=0.1, gamma=1.0, beta=2.0, L0=10):
    """Generate a Murray tree of m (int) generations. 

    Default corresponds to the single bifurcation test data.

    Starting one branch for m = 0. Given r is the radius of the mother
    branch. Radii of daughter branches (j, k) with parent i is
    goverend by gamma, such that

    r_o_j + r_o_k = r_o_i.

    and

    r_o_j/r_o_k = gamma

    i.e 

    r_o_j = gamma r_o_k
    """
    
    # Compute the number of elements/branches
    N = int(sum([2**i for i in range(m+1)]))

    # Create the generations
    generations = [[0,]]
    for i in range(1, m+1):
        _N = int(sum([2**j for j in range(i)]))
        siblings = list(range(_N, _N + 2**i))
        generations.append(siblings)

    # Create the indices (bifurcations)
    indices = []
    for (i, generation) in enumerate(generations[:-1]):
        children = generations[i+1]
        for (j, parent) in enumerate(generation):
            (sister, brother) = (children[2*j], children[2*j+1])
            indices += [(parent, sister, brother)]
    n = len(indices)
    assert N == (2*n+1), "N = 2*n+1"

    # Iterate through generations from children and upwards by
    # reversing the generations list for creating the paths
    generations.reverse()

    # Begin creating the paths
    ends = generations[0]
    paths = [[i,] for i in ends]
    for (i, generation) in enumerate(generations[1:]):
        for (j, path) in enumerate(paths):
            end_node = path[0]
            index = j // (2**(i+1))
            path.insert(0, generations[i+1][index])

    # Reverse it back
    generations.reverse()

    # Create inner radii based on Murray's law with factor gamma
    r_o = numpy.zeros(N)
    r_o[0] = r  
    for (g, generation) in enumerate(generations[1:]):
        for index in generation[::2]:
            ri = r_o[0]/(2**g)
            rk = 1/(1+gamma)*ri
            rj = gamma*rk
            r_o[index] = rk
            r_o[index+1] = rj
            
    L = L0*r_o
    r_e = beta*r_o
    
    # Peristaltic wave parameters: wave length lmbda and (angular) wave number k
    f = 1.0                 # frequency (Hz = 1/s)
    omega = 2*pi*f          # Angular frequency (Hz)
    lmbda = 2.0             # mm    
    k = 2*pi/lmbda          # wave number (1/mm)
    varepsilon = 0.1        # AU 
    
    data = (indices, paths, r_o, r_e, L, k, omega, varepsilon)
    return data
    
def Qprime(Q, varepsilon, omega, L, k, Delta, r_o):
    scale = 2*np.pi*varepsilon*omega*r_o**2/k
    val = scale*Q
    return val

def run_single_bifurcation_test():

    print("Running single bifurcation test case via general solution algorithm")
    data = single_bifurcation_data()
    (indices, paths, r_o, r_e, L, k, omega, varepsilon) = data
    
    (P, dP, avg_Q_1) = solve_bifurcating_tree(data)

    print("P = ", P)
    print("dP = ", dP)
    print("<Q_1> = ", avg_Q_1)
    print("eps*<Q_1_0> = %.3e" % (varepsilon*avg_Q_1[0]))
    
    beta0 = _beta(r_e[0], r_o[0])
    delta0 = _delta(beta0)
    Q10 = varepsilon*avg_Q_1[0]
    Qp = Qprime(Q10, varepsilon, omega, L[0], k, delta0, r_o[0])
    print("eps*<Q_1_0>' (mm^3/s) = %.3e" % Qp)

def run_three_junction_verification():
    """Verification test: Compare explicit method and general method on a
    three junction network and check that the computed flows are the same.
    """
    
    data = three_junction_data()
    (indices, paths, r_o, r_e, L, k, omega, varepsilon) = data

    print("Solving three junction test case via general algorithm")
    (P_g, dP_g, avg_Q_1_g) = solve_bifurcating_tree(data)

    print("Solving three-junction case explicitly")
    (P_e, dP_e, avg_Q_1_e) = solve_three_junction(data)

    debug("eps*<Q_1_0> = %.3e" % (varepsilon*avg_Q_1_e[0]))

    a = numpy.array(avg_Q_1_g)
    b = numpy.array(avg_Q_1_e)
    d = numpy.max(numpy.abs(a-b))
    
    print("Explicit method gives", a)
    print("General method gives", b)
    print("Max difference: ", d)

    assert d < 1.e-10, "Difference is larger than expected!"

def test_murray_data():

    debug("Murray tree data (m=1)")
    data = murray_tree_data(m=1, r=0.1, gamma=1.0, beta=2.0, L0=10)
    debug(data)

    debug("Single bifurcation data")
    data = single_bifurcation_data()
    debug(data)

    debug("Murray tree data (m=2)")
    data = murray_tree_data(m=2, r=0.1, gamma=1.0, beta=2.0, L0=10)
    debug(data)

    debug("Three junction data")
    data = three_junction_data()
    debug(data)

def run_murray_tree():
    print("Solving Murray tree.")

    data = murray_tree_data(m=2, r=0.1, gamma=1.0, beta=2.0, L0=10)
    (indices, paths, r_o, r_e, L, k, omega, varepsilon) = data

    (P, dP, avg_Q_1) = solve_bifurcating_tree(data)
    
    print("P = ", P)
    print("dP = ", dP)
    print("<Q_1> = ", avg_Q_1)
    print("eps*<Q_1_0> = %.3e" % (varepsilon*avg_Q_1[0]))
    
    beta0 = _beta(r_e[0], r_o[0])
    delta0 = _delta(beta0)
    Q10 = varepsilon*avg_Q_1[0]
    Qp = Qprime(Q10, varepsilon, omega, L[0], k, delta0, r_o[0])
    print("eps*<Q_1_0>' (mm^3/s) = %.3e" % Qp)

def test():

    print("")
    run_single_bifurcation_test()

    print("")
    run_three_junction_verification()

    print("")
    test_murray_data()

    print("")
    run_murray_tree()

if __name__ == "__main__":

    # FIXME: Update to include varying radius estimates. WRONG WITHOUT
    # IT. See Gjerde, Sanchez, Rognes, JAP (2023):
    # https://arxiv.org/pdf/2310.02429.pdf Section II and onwards.
    
    __DEBUG_MODE=False
    test()


