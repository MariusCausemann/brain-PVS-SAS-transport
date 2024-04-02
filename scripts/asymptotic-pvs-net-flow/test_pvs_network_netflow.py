from pvs_network_netflow import *
from pvs_network_netflow import _beta, _delta

# =================================================================================
# Helper functions for setting up idealized networks included:
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
    k = 2*np.pi/lmbda          # wave number (1/mm)
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
    the daughter vessels are identical to that of the mother vessel in
    order to compare with analytical expression.

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
        
    r_o = [ro0, ro0, ro0, ro0, ro0, ro0, ro0]  # Inner radii (mm) for each element/edge
    r_e = [re0, re0, re0, re0, re0, re0, re0]  # Outer radii (mm) for each element/edge
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
    omega = 2*np.pi*f          # Angular frequency (Hz)
    lmbda = 2.0             # mm    
    k = 2*np.pi/lmbda          # wave number (1/mm)
    varepsilon = 0.1        # AU 
    
    data = (indices, paths, r_o, r_e, L, k, omega, varepsilon)
    return data

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

    return True
    
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

    success = d < 1.e-10
    assert success, "Difference is larger than expected!"

    return success
    
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

    return True
    
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

    return True

def test():

    print("")
    success = run_single_bifurcation_test()

    print("")
    success = run_three_junction_verification()

    exit()

    print("")
    success = test_murray_data()

    print("")
    success = run_murray_tree()
    
if __name__ == "__main__":

    # FIXME: Update to include varying radius estimates. WRONG WITHOUT
    # IT. See Gjerde, Sanchez, Rognes, JAP (2023):
    # https://arxiv.org/pdf/2310.02429.pdf Section II and onwards.
    
    test()


