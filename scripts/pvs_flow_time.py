from dolfin import *
from solver import read_vtk_network, as_P0_function
import os
import numpy as np

def pvs_flow_system(radius_f, tau, radius_ratio, f=Constant(0), g=Constant(0)):
    '''The bilinear form corresponding to Darcy on the graph'''
    mesh = radius_f.function_space().mesh()
    assert mesh.topology().dim() == 1

    cell = mesh.ufl_cell()
    
    Velm = FiniteElement('Discontinuous Lagrange', cell, 0)
    Qelm = FiniteElement('Lagrange', cell, 1)
    Welm = MixedElement([Velm, Qelm])
    W = FunctionSpace(mesh, Welm)

    u, p = TrialFunctions(W)
    v, q = TestFunctions(W)

    r1 = radius_f
    r2 = r1*radius_ratio
    mu = Constant(0.697e-3) # dynamic visosity of water
    A_PVS = np.pi * (r2**2 - r1**2)
    # https://www.frontiersin.org/articles/10.3389/fphy.2022.882260/full, eq. 23
    # https://www.biorxiv.org/content/10.1101/2021.09.23.461519v2.full.pdf, e.q 8- 10
    kappa = 0.125 * (r2*r2 + r1*r1 - (r2*r2 - r1*r1)/ ln(r2/r1))  

    K = kappa*A_PVS / mu

    # This is stationary for simplicity, we would get evolution in time
    # by solving this with time dep boundary conditions. Later we might
    # want to add modifications based on time dependent momentum eq
    # and include viscuss term on the rhs?
    
    # A_pvs*u = -K*grad(p) where u is the cross-sectional averaged velocity
    # -div(A_pvs*u) = -f
    a = (inner((1/K)*A_PVS*u, v)*dx + inner(v, dot(grad(p), tau))*dx
         + inner(A_PVS*u, dot(grad(q), tau))*dx)

    # NOTE: Flux bcs are part of weak form (g here represent -u.tau), pressure
    # bcs are set strongly
    L = -inner(f, q)*dx  + inner(g, q)*ds

    return a, L, W

def pvs_flow_system_t(radius_f, tau, radius_ratio, pprev, dt, f=Constant(0), g=Constant(0)):
    '''The bilinear form corresponding to Darcy on the graph'''
    mesh = radius_f.function_space().mesh()
    assert mesh.topology().dim() == 1

    cell = mesh.ufl_cell()
    
    Velm = FiniteElement('Discontinuous Lagrange', cell, 0)
    Qelm = FiniteElement('Lagrange', cell, 1)
    Welm = MixedElement([Velm, Qelm])
    W = FunctionSpace(mesh, Welm)

    u, p = TrialFunctions(W)
    v, q = TestFunctions(W)

    r1 = radius_f
    r2 = r1*radius_ratio
    # https://www.frontiersin.org/articles/10.3389/fphy.2022.882260/full, eq. 23
    mu = Constant(0.697e-3) # dynamic visosity of water
    A_PVS = np.pi * (r2**2 - r1**2)
    # https://www.frontiersin.org/articles/10.3389/fphy.2022.882260/full, eq. 23
    # https://www.biorxiv.org/content/10.1101/2021.09.23.461519v2.full.pdf, e.q 8- 10
    kappa = 0.125 * (r2*r2 + r1*r1 - (r2*r2 - r1*r1)/ ln(r2/r1))  

    K = kappa*A_PVS / mu
    
    # u = -K*grad(p)
    # -div(u) = -f
    a = (1/dt)*inner(p,q)*dx +  (inner((1/K)*u, v)*dx + inner(v, dot(grad(p), tau))*dx
         + inner(u, dot(grad(q), tau))*dx)
    # NOTE: Flux bcs are part of weak form (g here represent -u.tau), pressure
    # bcs are set strongly
    L = (1/dt)*inner(pprev, q)*dx -inner(f, q)*dx  + inner(g, q)*ds

    return a, L, W
    
# --------------------------------------------------------------------

if __name__ == '__main__':
    from xii import TangentCurve

    radius_ratio = 1.4
    mesh ,artery_radii, artery_roots = read_vtk_network("mesh/networks/arteries_smooth.vtk")
    radius_f = as_P0_function(artery_radii)
    
    # Grab the tangent of xii; bottom line is that this is vector valued
    # function on the network describing tangent to the edge. Orientation
    # is arbitrary as long same tau is used throught the code
    tau = TangentCurve(mesh)

    a, L, W = pvs_flow_system(radius_f, tau, radius_ratio, f=Constant(-1e-3))

    bc_in = DirichletBC(W.sub(1), 0, artery_roots, 2)

    A, b = assemble_system(a, L, [bc_in])

    wh = Function(W)
    solve(A, wh.vector(), b)
   
    uh_mag, ph = wh.split()

    QQ = tau.function_space()
    # L2 projection to get the flux
    uh, qq = Function(QQ), TestFunction(QQ)
    hK = CellDiameter(mesh)
    assemble((1/hK)*inner(uh_mag*tau, qq)*dx, uh.vector())
    
    dt = 0.0001 
    Niter = 100
    n = 0 
            
    ph.rename("p","p")
    uh.rename("u", "u")
    whprev = Function(W)
    #uhprev, phprev = whprev.split()
    whprev.assign(wh)
    uhprev, phprev = whprev.split()
    uhprevve, qq = Function(QQ), TestFunction(QQ)

    assemble((1/hK)*inner(uhprev*tau, qq)*dx, uhprevve.vector())

    vtkfile_1  = File('pvs_flow_time.pvd')
    uhprevve.rename("uh","t")

    vtkfile_1 << (uhprevve, float(n))

    while n < Niter:
        hK = CellDiameter(mesh)
        a, L, W = pvs_flow_system_t(radius_f, tau, radius_ratio, phprev, dt, f=Constant(0), g=Constant(0))
        bc_in = DirichletBC(W.sub(1), 0, artery_roots, 2)

        A, b = assemble_system(a, L,bcs = None)

        wh = Function(W)
        solve(A, wh.vector(), b)

        uh_mag, ph = wh.split()
        whprev.assign(wh)
        uhprev, phprev = whprev.split()
        #uhprev.assign(uh_mag) 
        n+=1
    

        QQ = tau.function_space()
    # L2 projection to get the flux
        uh, qq = Function(QQ), TestFunction(QQ)
        hK = CellDiameter(mesh)
        assemble((1/hK)*inner(uh_mag*tau, qq)*dx, uh.vector())
        uh.rename("uh","t")
        vtkfile_1 << (uh, float(n))


        os.makedirs("../results/pvs_flow", exist_ok=True)
        with XDMFFile('results/pvs_flow/pvs_flow.xdmf') as xdmf:
            xdmf.write_checkpoint(uh, "velocity")
        File("results/pvs_flow/pvs_flow_t.pvd") << uh, n
