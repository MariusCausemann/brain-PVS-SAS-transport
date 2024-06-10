# https://sci-hub.ru/https://doi.org/10.1029/92JB02858
import numpy as np


def filter_negative(f, tol=1E-10):
    '''Variant of Lenardic & Kaula filter. Ensuring positivity of f'''
    values = f.vector().get_local()
    # 1. The initial sum of all nodal C values is calculated a
    total_pre = np.sum(values)
    # 2. Nodal C values below 0 are set to 0 and the peak valu below 0 is assignedto the variable Cmin
    negative = values < (0 + tol)
    Cmin = np.min(values[negative])
    values[negative] = 0
    # 3. Nodal C values less than or equal to the absolute value
    # of Cmin  are set to 0.
    lt_Cmin = values < (abs(Cmin) + tol)
    values[lt_Cmin] = 0
    # 6. The sum of all nodal C values is calculated and as-
    total_post = np.sum(values)
    # 7. The number of nodal C values not I or 0 is assigned
    positive = values > 0
    values[positive] += (total_pre - total_post)/np.sum(positive)

    f.vector().set_local(values)

    return f

# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import *

    mesh = IntervalMesh(2048, 0, 10)

    # Define function spaces
    V = FunctionSpace(mesh, "DG", 1)

    b = Constant((1, ))
    # Parameters
    D = Constant(0)
    t_end = 10
    dt = 0.1

    # Define unknown and test function(s)
    v = TestFunction(V)
    u = TrialFunction(V)

    u0 = Function(V)
    g = Expression('(std::fabs(x[0]-1) < 0.2 + DOLFIN_EPS) ? 1: 0', degree=0)
    u0.interpolate(g)

    # STABILIZATION
    h = CellDiameter(mesh)
    n = FacetNormal(mesh)
    alpha = Constant(1)

    # ( dot(v, n) + |dot(v, n)| )/2.0
    bn = (dot(b, n) + abs(dot(b, n)))/2.0

    inlet_id = 2

    inlet = CompiledSubDomain("near(x[0], 0)")
    bm = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
    inlet.mark(bm, inlet_id)
    ds = Measure("ds", mesh, subdomain_data=bm)
    #g = Expression("(x[0] > 0.5 - w/2 && x[0] < 0.5 + w/2) ? -sin((x[0] - 0.5 - w/2) * pi/w) : 0", degree=2, w=0.3)


    def a(u,v) :
        # Bilinear form
        a_int = dot(grad(v), D*grad(u) - b*u)*dx

        upwindgrad = lambda d: conditional(bn("+") > 0, grad(d)("+"), grad(d)("-"))

        a_fac = D*(alpha/avg(h))*dot(jump(u, n), jump(v, n))*dS \
                - D*dot(upwindgrad(u), jump(v, n))*dS \
                - D*dot(jump(u, n), upwindgrad(v))*dS

        a_vel = dot(jump(v), bn('+')*u('+') - bn('-')*u('-') )*dS  + dot(v, bn*u)*ds

        a = a_int + a_fac + a_vel
        return a

    # Define variational forms
    A = (1/dt)*inner(u, v)*dx + a(u,v)
    b = (1/dt)*inner(u0,v)*dx + g*v*ds(inlet_id)

    # Create files for storing results
    file = File("temp/u.pvd")

    u = Function(V)
    problem = LinearVariationalProblem(A, b, u)
    solver  = LinearVariationalSolver(problem)
    solver.parameters['linear_solver'] = 'mumps'

    u.assign(u0)
    u.rename("u", "u")

    xh = V.tabulate_dof_coordinates().flatten()
    idx = np.argsort(xh)
    
    # Time-stepping
    t = 0.0    
    while t < t_end:

        # Compute
        solver.solve()
        # Save to file
        # Move to next time step
        u = filter_negative(u)
        
        u0.assign(u)
        t += dt
        file.write(u, t)
        print(u.vector().min(), u.vector().max())

        import matplotlib.pyplot as plt

        #u0_ = Function(V)
        #u0_.vector().axpy(1, u0.vector())
        #gh = filter_negative(u0_)
        
        plt.figure()
        plt.plot(xh[idx], u0.vector().get_local()[idx])
        # plt.plot(xh[idx], gh.vector().get_local()[idx])        
        plt.show()
