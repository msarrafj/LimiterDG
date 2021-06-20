# import pytest
from firedrake import *
import numpy as np
import matplotlib.pyplot as plt
from Limiter.flux_limiter import *
import subprocess
import sys


nx = 100  # three cases (nx, save_itr) = (20,1) (40,2) (80,4)
save_itr = 1
mesh = SquareMesh(nx, nx, 3)
iterations =nx/20 * 100
T = 0.5
dt = T / iterations
dtc = Constant(dt)
t = 0.

v = FunctionSpace(mesh, "DG", 1)
# ErrorSpace = FunctionSpace(mesh, "DG", 0)
Vt = FunctionSpace(mesh, "HDiv Trace", 0)
w = TestFunction(Vt)
fluxes = Function(Vt)
solMin = 0.0
solMax = 1.0
kSpace = FunctionSpace(mesh,"DG",0)
area_cell = Function(kSpace).interpolate(CellVolume(mesh)).vector().array()

# Material parameter
mu_o = 1
mu_w = 1
s_rw = 0.
s_ro = 0.

def s_star(s):
    return (s-s_rw)/(1-s_rw-s_ro)

def lmbda_w(s):
    return s_star(s)**2/mu_w

def lmbda_o(s):
    return (1-s_star(s))**2/mu_o

def F_w(D):
    return as_vector([ lmbda_w(D)/(lmbda_w(D) + lmbda_o(D)), (lmbda_w(D)/(lmbda_w(D) + lmbda_o(D)))*(1-5*(1-s_star(D))**2)  ])

def dFds(D):
    return as_vector([  ( (2*D)*(1-D) )/((2*D**2-2*D+1)**2) , (-2*D*(10*D**4-25*D**3+30*D**2-19*D+4))/((2*D**2-2*D+1)**2)  ])




# Dirichlet BC
bc = DirichletBC(v, Constant(0.), "on_boundary" ,method="geometric")
bcs = [bc]

# advection problem
x = SpatialCoordinate(mesh)
D_in =project(conditional(pow((x[0]-1.5),2)+pow((x[1]-1.5),2) < (0.5+1e-7) , 1., 0.),v)
D0 = interpolate(D_in, v)
DAvg = Function(kSpace).interpolate(D0)

phi = TestFunction(v)
D = Function(v)
D_sol = Function(v)

n = FacetNormal(mesh)
area = FacetArea(mesh)

n = FacetNormal(mesh)

h = CellDiameter(mesh)
h_avg = (h('+')+h('-'))/2
sigma = Constant(0.1)

F =  (1./dtc) * phi * (D-D0) * dx - dot(grad(phi) , F_w(D)) *  dx\
    + dot( jump(phi) , dot (avg( F_w(D)) , n('+')) + 0.5 * abs(dot( dFds(D0)('+') , n('+') )) * jump(D)  ) * dS



# Make slope limiter
slopelimiter = VertexBasedLimiter(v)
slopelimiter.apply(D0)

outfile = File("./Output/DG+FL+SL_nx%d.pvd"%nx)
D0.rename("Saturation+FL+SL","Saturation+FL+SL")
outfile.write(D0)

dD = TrialFunction(v)
counter = 0
DAvg0 = DAvg
while t < (T - dt / 2):

    t += dt
    counter += 1
    print("counter=", counter)

    DAvg0.assign(DAvg)

    J = derivative(F, D, dD)
    problem = NonlinearVariationalProblem(F,D,bcs,J)
    solver = NonlinearVariationalSolver(problem,solver_parameters=
                                            {
                                            'snes_type': 'newtonls',
                                            'snes_rtol': 1e-6,
                                            'snes_max_it': 200,
                                            # "snes_linesearch_type": "basic",
                                                # INNER LINEAR SOLVER
                                            'ksp_rtol': 1e-6,
                                            'ksp_max_it': 100,
                                            # Direct solver
                                            'ksp_type':'preonly',
                                            'mat_type': 'aij',
                                            'pc_type': 'lu',
                                            "pc_factor_mat_solver_type": "mumps",
                                            })
    solver.solve()
    D_sol.assign(D)
    DAvg = Function(kSpace).interpolate(D_sol)
    # ---------------------;
    # flux-limiter applied ;
    # ---------------------;
    F_flux = fluxes('+')*w('+')*dS + fluxes*w*ds - (\
        area * 1. * dot(w('+'),  dot(avg(F_w(D)) ,n('+')) + 0.5 * abs(dot(dFds(D0)('+'),n('+')) ) *jump(D) )*dS
    )

    solve(F_flux == 0, fluxes)
    # -----------------------------------------;
    # calculate local mass balance in each cell:
    # -----------------------------------------;
    FLUX = get_flux(mesh,fluxes).apply()
    # Apply flux limiter
    D_sol,convergedIter = flux_limiter(mesh,D0,D_sol,fluxes,solMax,solMin,dt).apply()
    # Apply slope limiter
    slopelimiter.apply(D_sol)
    D0.assign(D_sol)
    print("D_sol min is",D_sol.vector().array().min())
    print("D_sol max is",D_sol.vector().array().max())


    if counter % save_itr == 0.:
        outfile.write(D0)


