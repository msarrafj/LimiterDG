# import pytest
from firedrake import *
import numpy as np
import matplotlib.pyplot as plt
from Limiter.flux_limiter import *
import subprocess
import sys


nx = 25  # test cases are: (nx, save_itr) = (400,8) (200,4) (100,2) (50,1) (25,1/2)
save_itr = 0.5 * 0.75 # Save vtu every save_itr iterations
Length = 300
mesh = RectangleMesh(nx, 1 ,Length,1,quadrilateral = True)
iterations = nx * 672/200 * 0.75   # we compare results at 48;96;168 that corresponds to t=400 days, t=800 days and t=1400 days
T = 86400*1400.
dt = T / iterations
dtc = Constant(dt)
t = 0

# test function space
v = FunctionSpace(mesh, "DG", 1)
v_CG = FunctionSpace(mesh,"CG",1)
m = VectorFunctionSpace(mesh, "CG", 1)
Vt = FunctionSpace(mesh, "HDiv Trace", 0)
w = TestFunction(Vt)
fluxes = Function(Vt)
solMin = .10
solMax = .85
kSpace = FunctionSpace(mesh,"DG",0)
area_cell = Function(kSpace).interpolate(CellVolume(mesh)).vector().array()
mesh_val = mesh.coordinates.vector().array()
cell_map = v.cell_node_map().values
mesh_cell_map = mesh.coordinates.cell_node_map().values
row = np.size(mesh_cell_map)
row,col = mesh.cell_set.size , 2
CELL_x = np.zeros((row,col))
CELL_y = np.zeros((row,col))
for irow in range(row):
    for icol in range(col):
        CELL_x[irow,icol] = mesh_val[mesh_cell_map[irow,icol]*2] # x-coord
        CELL_y[irow,icol] = mesh_val[mesh_cell_map[irow,icol]*2+1] # y-coord

# Material parameter
mu_o = 1
mu_w = 1
s_rw = 0.1
s_ro = 0.15
# s_rw = 0.
# s_ro = 0.

def s_star(s):
    return (s-s_rw)/(1-s_rw-s_ro)

def lmbda_w(s):
    return s_star(s)**4/mu_w

def lmbda_o(s):
    return (1-s_star(s))**2*(1-s_star(s)**2)/mu_o

def f_w(D):
    return lmbda_w(D)/(lmbda_w(D) + lmbda_o(D))

def F_w(D):
    return as_vector([3e-7*lmbda_w(D)/(lmbda_w(D) + lmbda_o(D)) , 0 ])

def dFds(D):
    return as_vector([3e-7 * 1./(1-s_rw-s_ro) * (2 * D**3 * (D**3 - 3*D + 2)) / (2*D**3 -2*D+1)**2 , 0 ])
# advecting velocity
u0 = as_vector([3e-7, 0])
u = Function(m).interpolate(u0)
x = SpatialCoordinate(mesh)
# Initial Conditions

# DG_param
h_avg = Constant(Length/float(nx))
h = Constant(Length/float(nx))
SIG = 0.000001
D_L = Constant(1-s_ro)
bcs = []

# advection problem
phi = TestFunction(v)
D_in = Constant(s_rw)
D0 = interpolate( D_in, v)
DAvg = Function(kSpace).interpolate(D0)
D = Function(v)
D_sol = Function(v)
n = FacetNormal(mesh)
area = FacetArea(mesh)

a_mass = 0.2 * phi * D * dx
a_int = dot(grad(phi), -1 * F_w(D) ) * dx
a_flux = dot(jump(phi),  dot(avg(F_w(D)) ,n('+')) + 0.5 * abs(dot(dFds(D0)('+'),n('+')) ) *jump(D) )*dS  #best choice [standard upwind]
# BC imposed weakly
F = (1./dt) * a_mass + a_int + a_flux +\
       dot(F_w(D_L),n) * phi * ds(1)+\
       dot(F_w(D),n) * phi * ds(2) -\
       (1./dt) * 0.2 * phi * D0 * dx 
# Make slope limiter
slope_limiter = VertexBasedLimiter(v)
slope_limiter.apply(D0)

outfile = File("./Output/BL+FL+SL_2D_nx%d.pvd"%nx)
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
                                            'snes_rtol': 1e-5,
                                            'snes_max_it': 200,
                                            'ksp_rtol': 1e-5,
                                            'ksp_max_it': 100,
                                            # Direct solver
                                            'ksp_type':'preonly',
                                            'mat_type': 'aij',
                                            'pc_type': 'lu',
                                            "pc_factor_mat_solver_type": "mumps",
                                            # Iterative solvers
                                            })
    solver.solve()
    D_sol.assign(D)
    DAvg = Function(kSpace).interpolate(D_sol)
    # ---------------------;
    # flux-limiter applied ;
    # ---------------------;
    F_flux = fluxes('+')*w('+')*dS + fluxes*w*ds - (\
        area * 1. * dot(w('+'),  dot(avg(F_w(D)) ,n('+')) + 0.5 * abs(dot(dFds(D0)('+'),n('+')) ) *jump(D) )*dS + \
        area * dot(F_w(D_L),n) * w * ds(1)+\
        area * dot(F_w(D),n) * w * ds(2)\
    )

    solve(F_flux == 0, fluxes)

    # Applied flux limiter
    D_sol,convergedIter = flux_limiter(mesh,D0,D_sol,fluxes,solMax,solMin,dt).apply()

    slope_limiter.apply(D_sol)

    D0.assign(D_sol)
    if counter % save_itr == 0:
        D0.rename("Saturation+FL+SL","Saturation+FL+SL")
        outfile.write(D0)
        print("vtu saved!")

        D0_val= D0.vector().array()
        row , col = np.shape(cell_map)
        col = 2
        D0_cell = np.zeros((row,col))
        for irow in range(row):
            for icol in range(col):
                D0_cell[irow,icol] = D0_val[cell_map[irow,icol]]
        DG_val = np.zeros((2*nx,2))
        for j in range(nx):
            for i in range(2):
                DG_val[2*j+i,0] = CELL_x[j,i]
                DG_val[2*j+i,1] = D0_cell[j,i]
        if counter == 48*save_itr or counter == 96*save_itr or counter ==168*save_itr:
            np.savetxt('profile/BL_time_%d_nx_%d.dat'%(counter/save_itr,nx), DG_val, fmt=['%d','%1.5e'])



