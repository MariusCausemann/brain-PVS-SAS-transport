from dolfin import *
from solver import read_vtk_network, as_P0_function
import os

def pvs_flow_system(radius_f, tau, f, g=Constant(0)):
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

    K = radius_f  # FIXME: correct dependence

    # This is stationary for simplicity, we would get evolution in time
    # by solving this with time dep boundary conditions. Later we might
    # want to add modifications based on time dependent momentum eq
    # and include viscuss term on the rhs?
    
    # u = -K*grad(p)
    # -div(u) = -f
    a = (inner((1/K)*u, v)*dx + inner(v, dot(grad(p), tau))*dx
         + inner(u, dot(grad(q), tau))*dx)

    # NOTE: Flux bcs are part of weak form (g here represent -u.tau), pressure
    # bcs are set strongly
    L = -inner(f, q)*dx  + inner(g, q)*ds

    return a, L, W
    
# --------------------------------------------------------------------

if __name__ == '__main__':
    from xii import TangentCurve

    mesh ,artery_radii, artery_roots = read_vtk_network("../mesh/networks/arteries_smooth.vtk")
    radius_f = as_P0_function(artery_radii)
    
    # Grab the tangent of xii; bottom line is that this is vector valued
    # function on the network describing tangent to the edge. Orientation
    # is arbitrary as long same tau is used throught the code
    tau = TangentCurve(mesh)

    a, L, W = pvs_flow_system(radius_f, tau, f=Constant(0))


    # For the bcs just make something z dependent here for illustration
    bc_out = DirichletBC(W.sub(1), 0, artery_roots, 1)
    bc_in = DirichletBC(W.sub(1), 1, artery_roots, 2)

    A, b = assemble_system(a, L, [bc_out, bc_in])

    wh = Function(W)
    solve(A, wh.vector(), b)

    uh_mag, ph = wh.split()

    QQ = tau.function_space()
    # L2 projection to get the flux
    uh, qq = Function(QQ), TestFunction(QQ)
    hK = CellDiameter(mesh)
    assemble((1/hK)*inner(uh_mag*tau, qq)*dx, uh.vector())
    
    ph.rename("p","")
    uh.rename("u", "")
    os.makedirs("../results/pvs_flow", exist_ok=True)
    File('../results/pvs_flow/flux.pvd') << uh
    File('../results/pvs_flow/pressure.pvd') << ph