from dolfin import *
from solver import read_vtk_network, as_P0_function
import numpy as np
import os
from branch_marking import color_branches
 # ------------------------------------------------------------------

if __name__ == '__main__':
    from xii import TangentCurve

    radius_ratio = 1.4
    mesh ,artery_radii, artery_roots = read_vtk_network("../mesh/networks/arteries_smooth.vtk", rescale_mm2m=False)
    radius_f = as_P0_function(artery_radii)
   
    assert mesh.topology().dim() == 1

    cell = mesh.ufl_cell()
    
    Velm = FiniteElement('Discontinuous Lagrange', cell, 0)
    Qelm = FiniteElement('Lagrange', cell, 1)
    Welm = MixedElement([Velm, Qelm])
    W = FunctionSpace(mesh, Welm)

    u, p = TrialFunctions(W)
    v, q = TestFunctions(W)
    f=Constant(0)
    g=Constant(0)
    
    r1 = radius_f
    r2 = r1*radius_ratio
    mu = Constant(0.697e-3) # dynamic visosity of water
    A_PVS = np.pi * (r2**2 - r1**2)
    # https://www.frontiersin.org/articles/10.3389/fphy.2022.882260/full, eq. 23
    # https://www.biorxiv.org/content/10.1101/2021.09.23.461519v2.full.pdf, e.q 8- 10
    kappa = 0.125 * (r2*r2 + r1*r1 - (r2*r2 - r1*r1)/ ln(r2/r1))  

    K = kappa / mu

    # This is stationary for simplicity, we would get evolution in time
    # by solving this with time dep boundary conditions. Later we might
    # want to add modifications based on time dependent momentum eq
    # and include viscuss term on the rhs?


    # Grab the tangent of xii; bottom line is that this is vector valued
    # function on the network describing tangent to the edge. Orientation
    # is arbitrary as long same tau is used throught the code
    
    tau = TangentCurve(mesh)
    ds = Measure('ds', domain = mesh, subdomain_data = artery_roots) 

    # A_pvs u = -K*A_pvs*grad(p)
    # div(A_pvs u ) = f 
    # robin boundary condition: - K*A_PVS*\partial_s(p) = alpha*(p - p_o): 
    # ds(1) picks up the outlet nodes

    alpha = Constant(0.0) 
    p_o   = Constant(0.0) 
    f     = Constant(-1e-3 * 0.19 / 0.7) 
    ## read branch length
    with XDMFFile("branch_length_read.xdmf") as file:
        DG = FunctionSpace(mesh, "DG", 0)
        branch_length = Function(DG)
        file.read_checkpoint(branch_length, "branch_length")
    #File('branch_length_test.pvd') << branch_length


   
    a = (inner(A_PVS*u, v)*dx + inner(v, K*A_PVS*dot(grad(p), tau))*dx
         + inner(A_PVS*u, dot(grad(q), tau))*dx) + inner((branch_length)*A_PVS*alpha*p, q)*ds(1)

    # NOTE: Flux bcs are part of weak form (g here represent -u.tau), pressure
    # bcs are set strongly
    L = -inner(A_PVS*f, q)*dx  + inner(g, q)*ds + inner((branch_length)*A_PVS*p_o, q)*ds(1)

    # NOTE: Flux bcs are part of weak form (g here represent -u.tau), pressure
    # bcs are set strongly

    bc_in = DirichletBC(W.sub(1), 0, artery_roots, 2)

    A, b = assemble_system(a, L, [bc_in])

    wh = Function(W)
    solve(A, wh.vector(), b)

    uh_mag, ph = wh.split()

   
    from IPython import embed
    active_cell_f = MeshFunction('size_t', mesh, 1, 0)
    dx_active = Measure('dx', domain=mesh, subdomain_data=active_cell_f)

    QQ = tau.function_space()
    # L2 projection to get the flux
    uh, qq = Function(QQ), TestFunction(QQ)
    hK = CellDiameter(mesh)
    assemble((1/hK)*inner(uh_mag*tau, qq)*dx, uh.vector())
    

     # Find branch and mark the cells that make it up by unique color
    marking_branch, branch_colors, loop_colors = color_branches(mesh)
    
    Q    = FunctionSpace(mesh, 'DG', 0)
    q    = TestFunction(Q)
    foo  = Function(Q)   
    foo_values = foo.vector().get_local()

  
    mean_vel = (1/hK)*(1/branch_length)*inner(uh_mag*tau,uh_mag*tau)*dx_active(1)
    active_cell_f_array = active_cell_f.array()
    active = q*dx_active(1)
    cell_colors = marking_branch.array()
    for color in branch_colors:
        active_cell_f_array[cell_colors == color] = 1
        mean_velocities = assemble(mean_vel)
        foo_values[assemble(active).get_local() > 0] = mean_velocities 
        active_cell_f_array *= 0
    foo.vector().set_local(foo_values)

    print(np.size(np.unique(foo.vector().get_local()))) 
    print(np.size(np.unique(branch_length.vector().get_local()))) 

    File('averaged_vel.pvd') << foo

    ph.rename("p", "p")
    uh.rename("u", "u")
    os.makedirs("../results/pvs_flow", exist_ok=True)
    
    with XDMFFile('results/pvs_flow/pvs_flow_robin_vis.xdmf') as xdmf:
        xdmf.write(uh)

    #import pyvista as pv
    #import matplotlib.pyplot as plt
    #from plotting_utils import set_plotting_defaults

    #grid = pv.read("results/pvs_flow/pvs_flow_vis.xdmf").compute_cell_sizes()
    #grid = grid.point_data_to_cell_data()
    #grid["umag"] = np.linalg.norm(grid["u"], axis=1) * 1e3
    #uavg = np.average(grid["umag"], weights=grid["Length"])
    #umax = grid['umag'].max()
    #umed = np.median(grid["umag"])
    #set_plotting_defaults()

    #fig, ax = plt.subplots()
    #plt.hist(grid["umag"], weights=grid["Length"], density=True, bins=50,
    #         range=(0, 0.1)
    #         )
    #plt.axvline(uavg, color="red")
    #plt.axvline(umed, color="cyan")
    #plt.xlabel("PVS flow velocity (mm/s)")
    #plt.ylabel("relative frequency")
    #plt.tight_layout()
    #plt.text(0.5, 0.5, f"max: {umax:.3f} mm/s", transform=ax.transAxes)
    #plt.text(0.5, 0.6, f"avg: {uavg:.3f} mm/s", transform=ax.transAxes)
    #plt.text(0.5, 0.7, f"median: {umed:.3f} mm/s",transform=ax.transAxes )
    #plt.savefig("results/pvs_flow/velocity_histo.png")