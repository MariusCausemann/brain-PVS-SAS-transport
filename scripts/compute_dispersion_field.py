from fenics import *
import numpy as np
import typer


def get_dispersion_enhancement(csf_pressure_file:str, outfile:str):

    mesh = Mesh()
    with XDMFFile(csf_pressure_file) as file:
        file.read(mesh)
        DG = FunctionSpace(mesh, "DG", 0)
        pressure_csf = Function(DG)
        file.read_checkpoint(pressure_csf, "pressure")

    V = FunctionSpace(mesh, "CG", 1)
    p_cont = project(pressure_csf, V)
    gradp = sqrt(inner(grad(p_cont), grad(p_cont)))
    rho = 993 # kg/m^3
    nu = 7e-7 # m^2/s
    omega = 2*np.pi
    h = 3e-3 / 2
    P = gradp /(rho*omega*nu/h)
    alpha = np.sqrt(h**2 * omega / nu)
    R = project(P**2, DG)

    u,v = TrialFunction(V), TestFunction(V)
    R = Function(V)
    a = (Constant(1e-4)*inner(grad(u), grad(v)) + Constant(1)*u*v)*dx
    L = P**2*v*dx
    solve(a==L, R)
    R = interpolate(R, DG)
    assert R.vector().min() > 0

    with XDMFFile(outfile) as file:
        file.write_checkpoint(R, "R")


if __name__ == "__main__":
    typer.run(get_dispersion_enhancement)