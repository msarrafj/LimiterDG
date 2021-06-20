from firedrake import *
import numpy as np
import math,sys
from Limiter.flux_limiter import *
import time as tm

formulation = sys.argv[1]


save_itr = 1
T = 1000 * 450/750
num_steps =  1000
dt = T / num_steps # time step size
#=====================================;
#  Create mesh and identify boundary  ;
#=====================================;
Length = 200
Hight = 100
nx = 60
ny = 30
mesh = RectangleMesh(nx, ny ,Length,Hight,diagonal = 'crossed')
#==========================;
#  Define function spaces ;
#=========================;
pSpace = FunctionSpace(mesh,"DG" , 1)
sSpace = FunctionSpace(mesh,"DG" , 1)
wSpace = MixedFunctionSpace([pSpace,sSpace])
vSpace = VectorFunctionSpace(mesh, "DG", 1)

u0 = Function(wSpace)
kSpace = FunctionSpace(mesh,"DG",0)
solMax = .85
solMin = .20
Vt = FunctionSpace(mesh, "HDiv Trace", 0)
w = TestFunction(Vt)
fluxes = Function(Vt)
with u0.dat.vec as solution_vec:
  DoF = solution_vec.getSize()
area_cell = Function(kSpace).interpolate(CellVolume(mesh)).vector().array()
#===================================;
#  Define trial and test functions  ;
#===================================;
(p_0,s_0) = TrialFunction(wSpace) # for linear solvers
(z,v) = TestFunction(wSpace)
#=====================;
#  Define parameters  ;
#=====================;
K = Constant(1e-8)
x, y = SpatialCoordinate(mesh)

mu_w = Constant(2.5e-4)
mu_l = Constant(5e-3)

phi = Constant(0.2)
rhow = Constant(600)
rhol = Constant(1000)


# Gravity
g = Constant((0,+9.8))

time = 0
x = SpatialCoordinate(mesh)
p_an = Function(pSpace).interpolate(3e6*1-(2e6/100*x[0]))
s_an = Constant(0.2)
p0 = interpolate( p_an, pSpace)
s0 = interpolate( s_an, sSpace)
v0 = interpolate( as_vector((0,0)),vSpace)


s_rw = 0.20
s_ro = 0.15
R = 0.55
theta = 2
P_d = 1000

def s_star(s):
    return s

def Pc(s,s0):
    return conditional(gt(s_star(s0),R), P_d * s_star(s)**(-1./theta), P_d * R **(-1./theta) - P_d/theta * R **(-1-1./theta) * (s_star(s)-R))

def lmbda_w(s):
    return  K/mu_w *  s_star(s)**4
def lmbda_l(s):
    return ((1-s_star(s))**2)*(1-s_star(s)**(2.)) * K/mu_l



q_w = Constant(0)
q_l = Constant(0)
#=================================;
#  Dirichlet boundary conditions  ;
#=================================;
P_L = Constant(3e6)
P_R = Constant(1e6)
S_L = Constant(0.85)
bcs = []
#============================;
#   Define variational form  ;
#============================;
n = FacetNormal(mesh)
area = FacetArea(mesh)
# if triangular mesh is used
h = CellSize(mesh)
h_avg = (h('+')+h('-'))/2
# if quadrilateral mesh is used
# h_avg = Constant(Length/float(nx))
# h = Constant(Length/float(nx))
SIG = 1000
sigma = Constant(SIG)
sigma_bnd = Constant(10*SIG)
epsilon = Constant(1)
#====================;
#  March over time   ;
#====================;
du = TrialFunction(wSpace)
file_p = File("./Output/Delta400_FL+SL_pres.pvd")
file_s = File("./Output/Delta400_FL+SL_sat.pvd")
file_v = File("./Output/Delta400_FL+SL_vel.pvd")


a_s0 = s_0*v*dx
L_s0 = s_an * v * dx

a_p0 = p_0*z*dx
L_p0 = p_an * z * dx

a_init = a_p0 + a_s0
L_init = L_p0 + L_s0
params = {'ksp_type': 'preonly', 'pc_type':'lu',"pc_factor_mat_solver_type": "mumps" }
A = assemble(a_init, bcs=bcs, mat_type='aij')
b = assemble(L_init)
solve(A,u0,b,solver_parameters=params)
pSol,sSol = u0.split()
p0.assign(pSol)
s0.assign(sSol)
sAvg = Function(kSpace).interpolate(s0)



s0.rename("Saturation+FL+SL","Saturation")
file_s.write(s0)

p0.rename("Pressure","Pressure")
file_p.write(p0)

v0.rename("Velocity","Velocity")
file_v.write(v0)

# non-linear problem
u = Function(wSpace)
u.assign(u0)
(p,s) = (u[0],u[1]) # for non-linear solvers

# Method I, Using NIPG symmetry term. this is not correct as our test function
# will be non-linear
a_p1 = (1./dt)*( phi *  (1-s) * z) * dx +\
    lmbda_l(s) * inner((grad(p)+grad(Pc(s,s0))-rhol*g),grad(z)) * dx -\
    inner(avg(lmbda_l(s)* (grad(p)-rhol*g)) , jump(z,n)) * dS -\
    inner(lmbda_l(s) * (grad(p)-rhol*g) , n) * z * ds(1) -\
    inner(lmbda_l(s) * (grad(p)-rhol*g), n) * z * ds(2) -\
    inner(avg(lmbda_l(s)* (grad(Pc(s,s0))-rhol*g)) , jump(z,n)) * dS -\
    inner(lmbda_l(s) * (grad(Pc(s,s0))-rhol*g) , n) * z * ds(1) +\
    epsilon * inner(avg(lmbda_l(s) *  (grad(z))) ,jump(p,n)) * dS +\
    epsilon * inner(lmbda_l(s) * (grad(z)), n ) * p * ds(1) +\
    epsilon * inner(lmbda_l(s) *  (grad(z)), n ) * p * ds(2) +\
    epsilon * inner(avg(lmbda_l(s) *  (grad(z))) ,jump(Pc(s,s0),n)) * dS +\
    epsilon * inner(lmbda_l(s) * (grad(z)), n ) * Pc(s,s0) * ds(1) +\
    sigma/h_avg * jump(p) * jump(z)  * dS +\
    sigma_bnd/h * p * z * ds(1) +\
    sigma_bnd/h * p * z * ds(2) +\
    sigma/h_avg * jump(Pc(s,s0)) * jump(z)  * dS +\
    sigma_bnd/h * Pc(s,s0) * z * ds(1) -\
    epsilon * inner(lmbda_l(s) * (grad(z)) , n ) * P_L * ds(1) -\
    epsilon * inner(lmbda_l(s) * (grad(z)) , n ) * P_R * ds(2) -\
    epsilon * inner(lmbda_l(s) * (grad(z)) , n ) * Pc(S_L,S_L) * ds(1) 

L_p1 = (1./dt)*(phi *  (1-s0) * z) * dx +\
     q_l * z * dx +\
    sigma_bnd/h * P_L * z * ds(1) +\
    sigma_bnd/h * P_R * z * ds(2) +\
    sigma_bnd/h * Pc(S_L,S_L) * z * ds(1) 


# Mehtod 2: Upwinding
a_p2 = (1./dt)*(phi(p) *  (1-s) * z) * dx +\
    lmbda_l(s) * inner((grad(p) +grad(Pc(s,s0)) -rhol*g) ,grad(z)) * dx -\
    conditional( gt(inner(avg(grad(p0)+grad(Pc(s,s0))-rhol*g),n('+')),0) ,\
        inner(lmbda_l(s)('+') * (avg(grad(p0)+grad(Pc(s,s0))-rhol*g)) , jump(z,n)),\
        inner(lmbda_l(s)('-') * (avg(grad(p0)+grad(Pc(s,s0))-rhol*g)) , jump(z,n))\
                       ) * dS -\
    inner(lmbda_l(s) * (grad(p)+grad(Pc(s,s0))-rhol*g), n) * z * ds(1) -\
    inner(lmbda_l(s) * (grad(p)+grad(Pc(s,s0))-rhol*g), n) * z * ds(2) +\
    sigma/h_avg * jump(p) * jump(z)  * dS +\
    sigma_bnd/h * p * z * ds(1) +\
    sigma_bnd/h * p * z * ds(2)
L_p2 = (1./dt)*(phi *  (1-s0) * z) * dx +\
    q_l * z * dx +\
    sigma_bnd/h * P_L * z * ds(1) +\
    sigma_bnd/h * P_R * z * ds(2)

a_s = (1./dt) * phi*  s * v * dx  + \
    lmbda_w(s) * inner((grad(p)-rhow*g) , grad(v)) * dx -\
    conditional( gt(inner(avg(grad(p0)-rhow*g),n('+')),0) ,\
        inner(lmbda_w(s)('+') * avg(grad(p0)-rhow*g) , jump(v,n)),\
        inner(lmbda_w(s)('-') * avg(grad(p0)-rhow*g) , jump(v,n))\
                       ) * dS -\
    lmbda_w(S_L)* inner((grad(p)-rhow*g),n) * v * ds(1) -\
    lmbda_w(s)* inner((grad(p)-rhow*g),n) * v * ds(2) +\
    sigma/h_avg * jump(s) * jump(v)  * dS +\
    sigma_bnd/h * s * v * ds(1) 

L_s = (1./dt) * phi*  s0 * v * dx +\
    q_w * v * dx + \
    sigma_bnd/h * S_L * v * ds(1) 

Slope_limiter = VertexBasedLimiter(sSpace)
Slope_limiter.apply(s0)

if formulation == 'NIPG':
    F =  a_p1 - L_p1 + a_s - L_s
elif formulation == 'Upwind':
    F =  a_p2 - L_p2 + a_s - L_s
else:
    sys.exit('***%s does not exist *******\n'
          '****************************\n'
          '****************************\n'
          '****************************\n'
          '***************************'%formulation)


sAvg0 = sAvg
counter = 1
print('Step \t solve \t FL(no_fluxCost)\t FL(all) \t SL \t convergedIter\n ')
for nIter in range(num_steps):
    print ("=====================")
    print ('Time_step =%d'%counter)
    print ("=====================")
    sAvg0.assign(sAvg)
    #update time
    time += dt
    counter += 1

    solve_start = tm.time()
    q_w.time = time
    q_l.time = time

    J = derivative(F, u, du)
    #initial guess for solution is set to 1.0
    # u.assign(Constant(1))
    problem = NonlinearVariationalProblem(F,u,bcs,J)
    solver = NonlinearVariationalSolver(problem,solver_parameters=
                                            {
                                                #OUTER NEWTON SOLVER
                                            'snes_type': 'newtonls',
                                            'snes_rtol': 1e-5,
                                            'snes_max_it': 200,
                                            # "snes_linesearch_type": "basic",
                                            # 'snes_monitor': None,
                                            ## 'snes_view': None,
                                            # 'snes_converged_reason': None,
                                                # INNER LINEAR SOLVER
                                            'ksp_rtol': 1e-5,
                                            'ksp_max_it': 100,
                                            # Direct solver
                                            'ksp_type':'preonly',
                                            'mat_type': 'aij',
                                            'pc_type': 'lu',
                                            "pc_factor_mat_solver_type": "mumps",
                                            # Iterative solvers
                                            # 'pc_type': 'hypre',
                                            # 'hypre_type': 'boomeramg',
                                            # 'ksp_type' : 'fgmres',
                                            # 'ksp_gmres_restart': '100',
                                            # 'ksp_initial_guess_non_zero': True,
                                            # 'ksp_converged_reason': None,
                                            # 'ksp_monitor_true_residual': None,
                                            ## 'ksp_view':None,
                                            # 'ksp_monitor':None
                                            })
    solver.solve()

    pSol,sSol = u.split()
    solve_end = tm.time()

    sAvg = Function(kSpace).interpolate(sSol)

    # ---------------------;
    # flux-limiter applied ;
    # ---------------------;
    FL_start = tm.time()
    F_flux = fluxes('+')*w('+')*dS + fluxes*w*ds - (\
        area * -1. *conditional( gt(inner(avg(grad(p0)-rhow*g),n('+')),0.) ,\
        lmbda_w(sSol)('+') * inner(avg(grad(p0)-rhow*g) ,n('+')) * w('+'),\
        lmbda_w(sSol)('-') * inner(avg(grad(p0)-rhow*g) ,n('+')) * w('+')\
                       )* dS -\
    area * lmbda_w(S_L)* inner((grad(pSol)-rhow*g),n) * w * ds(1) -\
    area * lmbda_w(sSol)* inner((grad(pSol)-rhow*g),n) * w * ds(2) +\
    area * sigma/h_avg * jump(sSol) * w('+') * dS +\
    area * sigma_bnd/h * sSol * w * ds(1) -\
    area * sigma_bnd/h * S_L * w * ds(1) 
    )

    solve(F_flux == 0, fluxes)
    # Applied flux limiter
    FL_solve_start = tm.time()
    sSol,convergedIter = flux_limiter(mesh,s0,sSol,fluxes,solMax,solMin,dt).apply()
    FL_end = tm.time()

    SL_start = tm.time()
    # Slope-limiter applied
    Slope_limiter.apply(sSol)
    SL_end = tm.time()

    p0.assign(pSol)
    s0.assign(sSol)
    sAvg = Function(kSpace).interpolate(s0)
    # True velocity and pressure
    rhow_v = 1000
    g_v = Constant((0,-9.8))

    v0 = Function(vSpace).interpolate(-1*lmbda_w(s0)*(grad(p0)-rhow_v*g_v))



    v0.rename("Velocity","Velocity")
    p0.rename("Pressure","Pressure")
    s0.rename("Saturation+FL+SL","Saturation")
    if counter % 10 == 0:
        file_p.write(p0)
        file_s.write(s0)
        file_v.write(v0)

