from firedrake import *
import numpy as np
import math,sys
from Limiter.flux_limiter_well import *

formulation = "Upwind"


T = 1
nx = int(sys.argv[1])
num_steps = (nx)**2
dt = T / num_steps # time step size
dtc = Constant(dt)
#=====================================;
#  Create mesh and identify boundary  ;
#=====================================;
mesh = UnitSquareMesh(nx,nx, diagonal="right")
#==========================;
#  Define function spaces ;
#=========================;
pSpace = FunctionSpace(mesh,"DG" , 1)
sSpace = FunctionSpace(mesh,"DG" , 1)
wSpace = MixedFunctionSpace([pSpace,sSpace])
u0 = Function(wSpace)
kSpace = FunctionSpace(mesh,"DG",0)
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
K = Constant(1)
x, y = SpatialCoordinate(mesh)

mu_w = Constant(1)
mu_l = Constant(1)

phi = Constant(0.2)
rhow = Constant(1000)
rhol = Constant(1000)


# Gravity
g = Constant((0,0))

time = 0
x = SpatialCoordinate(mesh)
s_an =  Function(sSpace).interpolate(0.2*(2.+2.*x[0]*x[1]+cos(time+x[0]))) 
p_an =  Function(pSpace).interpolate(2.+pow(x[0],2)*x[1]-pow(x[1],2)+pow(x[0],2)*sin(time+x[1])\
           -1./6. * ( 2*cos(time) -2*cos(time+1) + 11 ))
p0 = interpolate( p_an, pSpace)
s0 = interpolate( s_an, sSpace)


def Pc(s,s0):
    return  50 * s **(-0.5)


def lmbda_w(s):
    return  s*s*K/mu_w
def lmbda_l(s):
    return (1-s)*(1-s)*K/mu_l



time_ = Constant(0)

#=================================;
#  Dirichlet boundary conditions  ;
#=================================;
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
SIG = 100
sigma = Constant(SIG)
sigma_bnd = Constant(10*SIG)
epsilon = Constant(1)
#====================;
#  March over time   ;
#====================;
du = TrialFunction(wSpace)

# s_ex = Function(sSpace)
# p_ex = Function(pSpace)
# file_p = File("./Output/pres.pvd")
# file_s = File("./Output/sat+FL+SL_Avg.pvd")
# file_p_ex = File("./Output/pres_ex.pvd")
# file_s_ex = File("./Output/sat_ex_Avg.pvd")
# file_q_w = File("./Output/q_w.pvd")

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

# p0.rename("Pressure","Pressure")
# file_p.write(p0,time=time)

p_ex = interpolate(p_an,pSpace)
# p_ex.rename("Pressure","Pressure")
# file_p_ex.write(p_ex,time=time)

# s0.rename("Saturation","Saturation")
# file_s.write(s0,time=time)

s_ex = interpolate(s_an,sSpace)
s_ex_avg = Function(kSpace).interpolate(s_ex)
# s_ex_avg.rename("Saturation","Saturation")
# file_s_ex.write(s_ex_avg,time=time)


# non-linear problem
u = Function(wSpace)
u.assign(u0)
(p,s) = (u[0],u[1]) # for non-linear solvers


Slope_limiter = VertexBasedLimiter(sSpace)


counter = 1
print('Step \t solve \t FL(no_fluxCost)\t FL(all) \t SL \t convergedIter\n ')
for nIter in range(num_steps):
    print ("=====================")
    print ('Time_step =%d'%counter)
    print ("=====================")

    q_w = Function(sSpace).interpolate(-1* phi* (0.2*sin(time_+x[0]))\
     -( 0.04*pow(2. + 2.*x[0]*x[1] + cos(time_ + x[0]),2)*(2*x[1] + 2*sin(time_ + x[1])) + \
     0.08*(2. + 2.*x[0]*x[1] + cos(time_ + x[0]))*(2.*x[1] - sin(time_ + x[0]))*(2*x[0]*x[1] + 2*x[0]*sin(time_ + x[1]))\
     +.16*x[0]*(2. + 2.*x[0]*x[1] + cos(time_ + x[0]))*(pow(x[0],2) - 2*x[1] + pow(x[0],2)*cos(time_ + x[1])) +\
     0.04*pow(2. + 2.*x[0]*x[1] + cos(time_ + x[0]),2)*(-2 - pow(x[0],2)*sin(time_ + x[1]))\
     ))

    q_l = Function(sSpace).interpolate(phi*0.2*sin(time_+x[0])  \
     -(pow(1 - 0.2*(2. + 2.*x[0]*x[1] + cos(time_ + x[0])),2)*(2*x[1] + 2*sin(time_ + x[1])) -\
         0.4*(1 - 0.2*(2. + 2.*x[0]*x[1] + cos(time_ + x[0])))*(2.*x[1] - sin(time_ + x[0]))*(2*x[0]*x[1] + 2*x[0]*sin(time_ + x[1]))\
         -0.8*x[0]*(1 - 0.2*(2. + 2.*x[0]*x[1] + cos(time_ + x[0])))*(pow(x[0],2) - 2*x[1] + pow(x[0],2)*cos(time_ + x[1])) +\
         pow(1 - 0.2*(2. + 2.*x[0]*x[1] + cos(time_ + x[0])),2)*(-2 - pow(x[0],2)*sin(time_ + x[1]))\
         )\
     +(\
         (1.118033988749895*(0. - cos(time_ + x[0]))*pow(1 - 0.2*(2. + 2.*x[0]*x[1] + cos(time_ + x[0])),  2))/ \
         pow(2. + 2.*x[0]*x[1] + cos(time_ + x[0]),1.5) - \
         (0.447213595499958*(1 - 0.2*(2. + 2.*x[0]*x[1] + cos(time_ + x[0])))*pow(2.*x[1] - sin(time_ + x[0]),2))/ \
         pow(2. + 2.*x[0]*x[1] + cos(time_ + x[0]),1.5) - \
         (1.6770509831248424*pow(1 - 0.2*(2. + 2.*x[0]*x[1] + cos(time_ + x[0])),2)*pow(2.*x[1] - \
         sin(time_ + x[0]),2))/ \
         pow(2. + 2.*x[0]*x[1] + cos(time_ + x[0]),2.5) \
         +(-1.788854381999832*pow(x[0],2)*(1 - 0.2*(2. + 2.*x[0]*x[1] + cos(time_ + x[0]))))/ \
         pow(2. + 2.*x[0]*x[1] + cos(time_ + x[0]),1.5) - \
         (6.708203932499369*pow(x[0],2)*pow(1 - 0.2*(2. + 2.*x[0]*x[1] + cos(time_ + x[0])),2))/ \
         pow(2. + 2.*x[0]*x[1] + cos(time_ + x[0]),2.5) \
         )*50)

    # q_w.rename("q_w","q_w")
    # file_q_w.write(q_w,time=time)
    #update time
    time += dt
    counter += 1
    time_.assign(time)

    s_an =  Function(sSpace).interpolate(0.2*(2.+2.*x[0]*x[1]+cos(time_+x[0]))) 
    p_an =  Function(pSpace).interpolate(2.+pow(x[0],2)*x[1]-pow(x[1],2)+pow(x[0],2)*sin(time_+x[1])\
           -1./6. * ( 2*cos(time_) -2*cos(time_+1) + 11 ))

    p_ex = interpolate(p_an,pSpace)
    # p_ex.rename("Pressure","Pressure")
    s_ex = interpolate(s_an,sSpace)
    # print("s_exact max",s_ex.vector().array().max())
    # print("s_exact min",s_ex.vector().array().min())
    # file_p_ex.write(p_ex,time=time)
    p_ex_avg = Function(kSpace).interpolate(p_ex)
    s_ex_avg = Function(kSpace).interpolate(s_ex)
    # s_ex_avg.rename("Saturation","Saturation")
    # file_s_ex.write(s_ex_avg,time=time)

    solMin = s_ex.vector().array().min()
    solMax = s_ex.vector().array().max()


    # Mehtod 2: Upwinding
    a_p2 = (1./dtc)*(phi(p) *  (1-s) * z) * dx +\
        lmbda_l(s) * inner((grad(p) +grad(Pc(s,s0)) -rhol*g) ,grad(z)) * dx -\
        conditional( gt(inner(avg(grad(p0)+grad(Pc(s,s0))-rhol*g),n('+')),0) ,\
            inner(lmbda_l(s)('+') * (avg(grad(p0)+grad(Pc(s,s0))-rhol*g)) , jump(z,n)),\
            inner(lmbda_l(s)('-') * (avg(grad(p0)+grad(Pc(s,s0))-rhol*g)) , jump(z,n))\
                           ) * dS -\
        inner(lmbda_l(s) * (grad(p)+grad(Pc(s,s0))-rhol*g), n) * z * ds +\
        sigma/h_avg * jump(p) * jump(z)  * dS +\
        sigma_bnd/h * p * z * ds 

    L_p2 = (1./dtc)*(phi *  (1-s0) * z) * dx +\
        q_l * z * dx +\
        sigma_bnd/h * p_an * z * ds

    a_s = (1./dtc) * phi*  s * v * dx  + \
        lmbda_w(s) * inner((grad(p)-rhow*g) , grad(v)) * dx -\
        conditional( gt(inner(avg(grad(p0)-rhow*g),n('+')),0) ,\
            inner(lmbda_w(s)('+') * avg(grad(p0)-rhow*g) , jump(v,n)),\
            inner(lmbda_w(s)('-') * avg(grad(p0)-rhow*g) , jump(v,n))\
                           ) * dS -\
        lmbda_w(s_an)* inner((grad(p)-rhow*g),n) * v * ds +\
        sigma/h_avg * jump(s) * jump(v)  * dS +\
        sigma_bnd/h * s * v * ds 

    L_s = (1./dtc) * phi*  s0 * v * dx +\
        q_w * v * dx + \
        sigma_bnd/h * s_an * v * ds

    F =  a_p2 - L_p2 + a_s - L_s


    J = derivative(F, u, du)
    #initial guess for solution is set to 1.0
    # u.assign(Constant(1))
    problem = NonlinearVariationalProblem(F,u,bcs,J)
    solver = NonlinearVariationalSolver(problem,solver_parameters=
                                            {
                                                #OUTER NEWTON SOLVER
                                            'snes_type': 'newtonls',
                                            'snes_rtol': 1e-6,
                                            'snes_max_it': 200,
                                            "snes_linesearch_type": "basic",
                                            # 'snes_monitor': None,
                                            ## 'snes_view': None,
                                            # 'snes_converged_reason': None,
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


    # ---------------------;
    # flux-limiter applied ;
    # ---------------------;
    wells_avg = Function(kSpace).interpolate(q_w)
    F_flux = fluxes('+')*w('+')*dS + fluxes*w*ds - (\
        area * -1. *conditional( gt(inner(avg(grad(p0)-rhow*g),n('+')),0.) ,\
       lmbda_w(sSol)('+') * inner(avg(grad(p0)-rhow*g) ,n('+')) * w('+'),\
        lmbda_w(sSol)('-') * inner(avg(grad(p0)-rhow*g) ,n('+')) * w('+')\
                      )* dS -\
    area * lmbda_w(s_an)* inner((grad(pSol)-rhow*g),n) * w * ds +\
    area * sigma/h_avg * jump(sSol) * w('+') * dS +\
    area * sigma_bnd/h * sSol * w * ds -\
    area * sigma_bnd/h * s_an * w * ds
    )

    solve(F_flux == 0, fluxes)
    p0.assign(pSol)
    s0.assign(sSol)

    p0_avg = Function(kSpace).interpolate(p0)
    s0_avg = Function(kSpace).interpolate(s0)

    # p0_avg.rename("Pressure","Pressure")
    # s0_avg.rename("Saturation","Saturation")
    # file_p.write(p0_avg)
    # file_s.write(s0_avg)

L2_error_s = errornorm(s_ex_avg,s0_avg,'L2')
L2_error_p = errornorm(p_ex_avg,p0_avg,'L2')
print("L2_error_of_s %e" % L2_error_s)
print("L2_error_of_p %e" % L2_error_p)

print('dofs', DoF)
