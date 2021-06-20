from firedrake import *
import numpy as np
import math,sys,time
from matplotlib import pyplot as plt
from Limiter.flux_limiter_well import *
from scipy.io import loadmat

formulation = sys.argv[1]
VertexLim = "Kuzmin"


T = 86400.*2.5           # final time(second) 2.5 day
num_steps = 600     # number of time steps

dt = T / num_steps # time step size
#=====================================;
#  Create mesh and identify boundary  ;
#=====================================;
domain_length = 1000.
nx = 50
layer = '13' # layer 13 and layer 73
mesh = RectangleMesh(nx, nx ,domain_length,domain_length,diagonal='crossed')
#==========================;
#  Define function spaces ;
#=========================;
pSpace = FunctionSpace(mesh,"DG" , 1)
sSpace = FunctionSpace(mesh,"DG" , 1)
wSpace = MixedFunctionSpace([pSpace,sSpace])
vSpace = VectorFunctionSpace(mesh,"DG" , 1)
u0 = Function(wSpace)
kSpace = FunctionSpace(mesh,"DG",0)
solMax = .85
solMin = .2
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
x = SpatialCoordinate(mesh)
q_I = interpolate(conditional(And(And(x[0]< 120,x[0]>20),And(x[1]<120,x[1]>20)) , 0.28e-4 , 0.0 ), kSpace)

q_P = interpolate(conditional(And(And(x[0]< 980,x[0]>880),And(x[1]<980,x[1]>880)) , 0.28e-4 , 0.0 ), kSpace)

# SPE10 loaded from file
class MyExpression(Expression):
    def Perm(self,x_val,y_val):
        data = loadmat('./K_Input/all/Kx_%s.mat'%layer)
        data_m = data['Kx_layer']
        x_new , y_new = x_val*(60./domain_length) , y_val*(220./domain_length)
        # x_new , y_new = x_val*(2./domain_length) , y_val*(2./domain_length)
        return data_m[math.floor(x_new),math.floor(y_new)]

    def eval(self, value, x):
        value[0] = self.Perm(x[0],x[1])

    def value_shape(self):
        return ()

K = interpolate(MyExpression(),kSpace)

File('perm.pvd').write(K)


mu_w = Constant(5e-4)
mu_l = Constant(2e-3)

phi_0 = Constant(0.2)
rhow_0 = Constant(1000)
rhol_0 = Constant(850)

# Gravity
g = Constant((0,0))

time = 0
p_an = Constant(0.)

s_an = Constant(0.2)
p0 = interpolate( p_an, pSpace)
s0 = interpolate( s_an, sSpace)
v0 = interpolate( as_vector((0,0)),vSpace)

s_rw = 0.2
s_ro = 0.15
R = 0.05
theta = 2
P_d = 4500

def s_star(s):
    return s

def Pc(s,s0):
    return conditional(gt(s_star(s0),R), P_d * s_star(s)**(-1./theta), P_d * R **(-1./theta) - P_d/theta * R **(-1-1./theta) * (s_star(s)-R))

def lmbda_w(s):
    return  s**(5) * 1./mu_w

def lmbda_l(s):
    return (1-s)*(1-s)*(1-s**(5)) * 1./mu_l

def phi(p):
    return phi_0

def rhow(p):
    return rhow_0

def rhol(p):
    return rhol_0
def f_w(s):
    return lmbda_w(s)/(lmbda_w(s)+lmbda_l(s))

def f_o(s):
    return lmbda_l(s)/(lmbda_w(s)+lmbda_l(s))
# q_l = Constant(0)

#=================================;
#  Dirichlet boundary conditions  ;
#=================================;
bcs = []
#============================;
#   Define variational form  ;
#============================;
n = FacetNormal(mesh)
area = FacetArea(mesh)
h = CellSize(mesh)
h_avg = (h('+')+h('-'))/2
# h = Constant(Length/float(nx))
# h_avg = Constant(Length/float(nx))
SIG = 10
sigma = Constant(SIG)
sigma_bnd = Constant(10*SIG)
epsilon = Constant(1)
#====================;
#  March over time   ;
#====================;
du = TrialFunction(wSpace)
file_p = File("./Output/pres+FL+SL.pvd")
file_s = File("./Output/sat+FL+SL.pvd")
file_v = File("./Output/vel+FL+SL.pvd")
file_error = File("./Output/err+FL+SL.pvd")


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
error = Function(kSpace).interpolate(Constant(0.0))

error.rename("Local mass balance error","error")
file_error.write(error)

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
a_p1 = (1./dt)*(phi_0 *   (1-s) * z) * dx +\
     lmbda_l(s) * inner(( K * grad(p) - K * rhol_0*g),grad(z)) * dx -\
    inner(avg(lmbda_l(s)* (K * grad(p) - K * rhol_0*g)) , jump(z,n)) * dS +\
    epsilon * inner(avg(lmbda_l(s) *  (K * grad(z))) ,jump(p,n)) * dS +\
    sigma/h_avg * jump(p) * jump(z)  * dS 

L_p1 = (1./dt) * phi_0 *  (1-s0) * z * dx +\
    ( f_o(0.85) * q_I -  f_o(s0) * q_P) * z * dx 


# Mehtod 2: Upwinding
a_p2 = (1./dt)*(phi_0 *  (1-s) * z) * dx +\
     lmbda_l(s) * inner((K * grad(p) - K * rhol_0*g),grad(z)) * dx -\
    conditional( gt(inner(avg((K * grad(p0)- K * rhol_0*g)),n('+')),0) ,\
        inner(lmbda_l(s)('+') * avg((K * grad(p0)- K * rhol_0*g)) , jump(z,n)),\
        inner(lmbda_l(s)('-') * avg((K * grad(p0)- K * rhol_0*g)) , jump(z,n))\
                       ) * dS +\
    sigma/h_avg * jump(p) * jump(z)  * dS

L_p2 = (1./dt)*(phi_0 *  (1-s0) * z) * dx +\
    ( f_o(0.85) * q_I -  f_o(s0) * q_P) * z * dx


a_s = (1./dt) * phi_0*  s * v * dx  + \
     lmbda_w(s) * inner(( K * grad(p) - K * rhow_0*g) , grad(v)) * dx -\
    conditional( gt(inner(avg((K * grad(p0)- K * rhow_0*g)),n('+')),0) ,\
        inner(lmbda_w(s)('+') * avg(( K * grad(p0)- K * rhow_0*g)) , jump(v,n)),\
        inner(lmbda_w(s)('-') * avg(( K * grad(p0)- K * rhow_0*g)) , jump(v,n))\
                       ) * dS +\
    sigma/h_avg * jump(s) * jump(v)  * dS 

L_s = + (1./dt) * phi_0*  s0 * v *  dx+\
      (f_w(0.85) * q_I - f_w(s0)*q_P) * v * dx 

if VertexLim == 'Kuzmin':
    # Make slope limiter
    limiter = VertexBasedLimiter(sSpace)
    limiter.apply(s0)
elif VertexLim == 'None':
    print('NO LIMITER Used')
else:
    sys.exit('***Not a valid limiter***\n'
          '****************************\n'
          '****************************\n'
          '****************************\n'
          '***************************')

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
for nIter in range(num_steps):
    print ("=====================")
    print ('Time_step =%d'%counter)
    print ("=====================")
    sAvg0.assign(sAvg)
    #update time
    time += dt
    counter += 1
    # q_w.time = time
    # q_l.time = time

    # p_an.time=time
    # s_an.time=time
    J = derivative(F, u, du)
    #initial guess for solution is set to 1.0
    # u.assign(Constant(1))
    problem = NonlinearVariationalProblem(F,u,bcs,J)
    solver = NonlinearVariationalSolver(problem,solver_parameters=
                                            {
                                                #OUTER NEWTON SOLVER
                                            'snes_type': 'newtonls',
                                            'snes_atol': 1e-6,
                                            'snes_max_it': 200,
                                            "snes_linesearch_type": "basic",
                                            # "snes_linesearch_norms":True,
                                            # "snes_monitor_lg_residualnorm":None,
                                            # 'snes_monitor_short': None,
                                            # 'snes_view': None,
                                            # 'snes_converged_reason': None,
                                            # 'snes_monitor_residual': "ascii:res.h5::append", # print residual vector
                                            # 'snes_monitor_residual':None, # print residual vector
                                            # 'snes_monitor_lg_range':None,
                                                # INNER LINEAR SOLVER
                                            'ksp_rtol': 1e-6,
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
    sAvg = Function(kSpace).interpolate(sSol)
    wells_avg = Function(kSpace).interpolate(( f_w(0.85) * q_I - f_w(s0)*q_P ))


    F_flux = fluxes('+')*w('+')*dS + fluxes*w*ds - (\
        area * -1. *conditional(\
        gt(inner(avg((K * grad(p0)- K * rhow_0*g)),n('+')),0.) ,\
        lmbda_w(sSol)('+') * inner(avg(( K * grad(p0)- K * rhow_0*g)) ,n('+')) * w('+'),\
        lmbda_w(sSol)('-') * inner(avg(( K * grad(p0)- K * rhow_0*g)) ,n('+')) * w('+')\
                       )* dS +\
    area * sigma/h_avg * jump(sSol) * w('+')  * dS 
    )
    solve(F_flux == 0, fluxes)
    # -----------------------------------------;
    # calculate local mass balance in each cell:
    # -----------------------------------------;
    FLUX = get_flux(mesh,fluxes).apply()
    # print('FLUX',FLUX)
    error_val =   sAvg.vector().array() -\
            ( sAvg0.vector().array() -\
            (1./(0.2)) * (dt/area_cell) * FLUX.sum(axis=1)) -(1./(0.2)) * dt * wells_avg.vector().array()

    u_error = Function(kSpace)
    u_error.vector().set_local(error_val)
    error = Function(kSpace).interpolate(u_error)

    # flux limiter applied
    sSol,convergedIter = flux_limiter(mesh,s0,sSol,wells_avg,q_I,q_P,fluxes,solMax,solMin,dt).apply()

    # Slope-limiter applied
    limiter.apply(sSol)

    p0.assign(pSol)
    s0.assign(sSol)
    sAvg = Function(kSpace).interpolate(s0)
    v0 = Function(vSpace).interpolate(-1*K*lmbda_w(s0)*(grad(p0)-rhow_0*g))

    p0.rename("Pressure","Pressure")
    s0.rename("Saturation+FL+SL","Saturation")
    v0.rename("Velocity","Velocity")
    error.rename("Local mass balance error","error")

    file_p.write(p0)
    file_s.write(s0)
    file_v.write(v0)
    file_error.write(error)


print('dofs=', DoF)
