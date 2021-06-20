from firedrake import *
import numpy as np
import math,sys
import time as tm
from matplotlib import pyplot as plt
# from Limiter.hsign_backup import *
from Limiter.hsign import *
np.set_printoptions(precision=9)
from warnings import *

class get_flux(SignFlip):
    def __init__(self,mesh,fluxes):
        self.mesh = mesh
        self.fluxes = fluxes

    def apply(self):
        hsign,cell_map,TwoCell,boundary_facet = SignFlip(self.mesh).HSIGN()
        #========================;
        # Step 0: Define H_{E,E'}:
        #========================;
        flux_val= self.fluxes.vector().array()
        row , col = np.shape(cell_map)
        old_FLUX = np.zeros((row,col))
        for irow in range(row):
            for icol in range(col):
                old_FLUX[irow,icol] = flux_val[cell_map[irow,icol]]
        # print('Old flux:\n',old_FLUX)
        FLUX = np.multiply(old_FLUX,hsign)
        return FLUX


class flux_limiter(SignFlip):
    def __init__(self,mesh,sol_0,sol,fluxes,solMax,solMin,dt):
    # def __init__(self,mesh,sol_0,sol,wells_avg,q_I,q_P,fluxes,solMax,solMin,dt):
        self.mesh = mesh
        self.sol_0 = sol_0
        self.sol = sol
        # self.wells_avg = wells_avg.vector().array()
        # self.q_I = q_I.vector().array()
        # self.q_P = q_P.vector().array()
        self.fluxes = fluxes
        # self.rPAvg_0 = rPAvg_0.vector().array()
        # self.rPAvg = rPAvg.vector().array()
        self.solMax = solMax
        self.solMin = solMin
        self.dt = dt


    # def wells_update(self,s):
    #     fw = lambda t:\
    #     ( ((t-0.15)/0.7) * ((t-0.15)/0.7) * 1./(5e-4) )\
    #     /( ( ((t-0.15)/0.7) * ((t-0.15)/0.7) * 1./(5e-4) ) +\
    #     ( (1-((t-0.15)/0.7)) *(1-((t-0.15)/0.7) ) * 1./(2e-3) ))
    #     return (1000* self.q_I - 1000 * fw(s) * self.q_P) 
    #     # self.wells_update(s) is when we want to use this function in this class
    #     # forexample self.wells_update(solAvg)


    def apply(self):
        hsign,cell_map,TwoCell,boundary_facet = SignFlip(self.mesh).HSIGN()
        # print('Hsign', signing)
        # print('cell_map\n', cell_map)
        # print('TwoCell', TwoCell)
        # print('boundary_facet', boundary_facet)

        #========================;
        # Step 0: Define H_{E,E'}:
        #========================;
        flux_val= self.fluxes.vector().array()
        row , col = np.shape(cell_map)
        old_FLUX = np.zeros((row,col))
        for irow in range(row):
            for icol in range(col):
                old_FLUX[irow,icol] = flux_val[cell_map[irow,icol]]
        # print('Old flux:\n',old_FLUX)
        FLUX = np.multiply(old_FLUX,hsign)
        # print('original FLUX:\n',FLUX)


        #========================================;
        # Calculate element correction factors:
        #========================================;

        # container to store alphaEplus and alphaEminu
        Vdg0 = FunctionSpace(self.mesh, "DG", 0)
        V0 = Function(Vdg0)
        # print(type(self.solMax))
        deltaBound = 10*np.finfo(self.solMax).eps #Arithmetic perturbations may violate the bounds.
        epsBound = np.finfo(self.solMax).eps
        # print('delatabound',deltaBound)
        area = V0.interpolate(CellVolume(self.mesh)).vector().array()
        # print(area)
        solAvgDG0_0 = V0.interpolate(self.sol_0).vector().array() 
        # print('sol DG0 avg previous time step:\n',solAvgDG0_0)
        # print('sol DG0 avg previous time step: min:%f,max:%f'%(solAvgDG0_0.min(),solAvgDG0_0.max()))

        # assert (np.all((solAvgDG0_0 <= self.solMax+deltaBound ) & (solAvgDG0_0 >= self.solMin-deltaBound) ) )\
        #             ,"solAvg of previous time step is not  between min and max!"

        #==================;
        # Start iteration  ;
        #==================;
        maxIter = 10000
        epsFLUX = 1e-6
        epsDelta = 1e-6
        solAvg = solAvgDG0_0 # limiter solution from last time step s0 on DG0 space
        for numIter in range(maxIter):
            # print('numIter of Flux-limiter',numIter)
           
            # Save suppressed fluxes.
            FLUX_0 = FLUX
            alphaEplus= np.zeros((row,1))
            # alphaEplus  = min(1, g.areaT .* max(0, (umax - deltaBound) - lowOrderMeans) ./ (-sum( supprFluxes .* (supprFluxes < 0), 2) * tau + eps) );

            # Qplus = (area * np.maximum(0,np.subtract((self.solMax-deltaBound),solAvg)) )
            # incompressible patch test
            Qplus = area*np.maximum(0, (0.2)*\
            np.subtract((self.solMax-deltaBound),solAvg)) 
            
            # Qplus = (area *(0.2 * 1000)* np.maximum(0,np.subtract((self.solMax-deltaBound),solAvg)) )
            # incompressible well problem
            # Qplus = area*np.maximum(0, (0.2*1000)*\
            # np.subtract((self.solMax-deltaBound),solAvg)-self.wells_avg*self.dt) 
            # Qplus = area*np.maximum(0, (0.2*1000)*\
            # np.subtract((self.solMax-deltaBound),solAvg)+\
            # self.wells_update(solAvg)*self.dt) 

            # compressible patch test
            # Qplus = area * np.maximum(0,np.subtract\
            #         ((self.rPAvg * self.solMax-deltaBound),self.rPAvg_0 * solAvg)) 
            # print('Qplus:\n',Qplus)
            assert (np.all(Qplus >= 0))\
                    ,"Qplus is less than zero!"
            Pplus = -1 *np.where(FLUX<0,FLUX,0).sum(axis=1) * self.dt + epsBound
            # print('Pplus:\n',Pplus)
            alphaEplus = np.minimum(1,np.divide(Qplus,Pplus))
            # print('alphaEplus at iteration %d:\n'%numIter,alphaEplus)
            # if alphaEplus is 1  means that no limiting is needed. (100% of unlimited flux is allowed
            # without introducing a mean-value overshoot); If alphaEplus is 0 this means that:
            # means Qplus = 0 which means  no mass is allowed to be stored in E without introducing
            # mean-value overshoot.
          
            # incompressible patch test
            Qminus = area *(np.minimum(0,(0.2)*\
            np.subtract((self.solMin+deltaBound),solAvg )))
            #
            # Qminus = (area * np.maximum(0,np.subtract(solAvg,self.solMin+deltaBound) ))
            # Qminus = (area *(0.2*1000)* np.maximum(0,np.subtract(solAvg,self.solMin+deltaBound) ))
            # incompressible well problem
            # Qminus = area *(np.minimum(0,(0.2*1000)*\
            # np.subtract((self.solMin+deltaBound),solAvg)-self.wells_avg*self.dt) )
            # Qminus = area *(np.minimum(0,(0.2*1000)*\
            # np.subtract((self.solMin+deltaBound),solAvg ) -\
            # self.wells_update(solAvg)*self.dt) )

            # compressible patch test
            # Qminus = area * np.minimum(0,np.subtract\
            #         ((self.rPAvg * self.solMin+deltaBound),self.rPAvg_0 * solAvg ))
            #
            # Qminus = (area * np.subtract(solAvg,(self.solMin-deltaBound) ) )
            # Qminus = (area * np.subtract(solAvg,(self.solMin) ) )
            # assert (np.all(Qminus >= 0))\
            #         ,"Qminus is less than zero!"
            assert (np.all(Qminus <= 0))\
                    ,"Qminus is greater than zero!"
            # Qminus = (area * np.subtract(solAvg,(self.solMin+deltaBound)) )
            # print('Qminus:\n',Qminus)
            # Pminus = 1 *np.where(FLUX>0,FLUX,0).sum(axis=1) * self.dt + epsBound
            Pminus = -1 *np.where(FLUX>0,FLUX,0).sum(axis=1) * self.dt - epsBound
            # print('Pminus:\n',Pminus)
            alphaEminus = np.minimum(1,np.divide(Qminus,Pminus))
            # print('alphaEminus at iteration %d:\n' %numIter,alphaEminus)
            # alphaEminus shows the percentage of howmuch of mass (Pminus) is allowd to exit element E
            # if alphaEminus = 1 no limiting is needed. and alphaEminus = 0 
            # it means that no mass is allowed to exit and hence 100% of flux should be limited.
            # np.savetxt('Qminus.out',Qminus)

            #============================================;
            # Compute edge correction factors alpha_E,E' :
            #============================================;
            # met1_Start = tm.time()
            alphaEface = np.ones((row,col))
            for irow in range(row):
                for icol in range(col):
                    facet = cell_map[irow,icol]
                    # Handling boundary terms
                    if facet in boundary_facet:
                        if FLUX[irow,icol] < 0:
                            alphaEface[irow,icol] = alphaEplus[irow] 

                        elif FLUX[irow,icol] > 0:
                            alphaEface[irow,icol] = alphaEminus[irow] 


                    # Handling interior edges
                    else:
                        if FLUX[irow,icol] < 0:
                            b0 = TwoCell[facet] # cellID of irow and the opposite cell
                            oppCell_ID = int(b0[np.where( b0 != irow )]) # includes only oppoiste cell ID
                            alphaEface[irow,icol] = np.minimum(alphaEplus[irow] , alphaEminus[oppCell_ID])

                        elif FLUX[irow,icol] > 0:
                            b0 = TwoCell[facet] # cellID of irow and the opposite cell
                            oppCell_ID = int(b0[np.where( b0 != irow )]) # includes only oppoiste cell ID
                            alphaEface[irow,icol] = np.minimum(alphaEminus[irow] , alphaEplus[oppCell_ID])


            # met1_End = tm.time()
            # print('alphaEface at numIter %d:\n'%numIter,alphaEface,type(alphaEface))

            # met2_Start = tm.time()
            # # use list comprehension in python which is very fast
            # alphaEface_fast=[alphaEplus[irow]*np.where(FLUX[irow,icol]<0,1,0) +
            #             alphaEminus[irow]*np.where(FLUX[irow,icol]>0,1,0)
            #             if cell_map[irow,icol] in boundary_facet
            #             else np.minimum(alphaEplus[irow] ,
            #     alphaEminus[ int(TwoCell[cell_map[irow,icol]][np.where( TwoCell[cell_map[irow,icol]] != irow )]) ])
            #             *np.where(FLUX[irow,icol]<0,1,0) + np.minimum(alphaEminus[irow] ,
            #     alphaEplus[ int(TwoCell[cell_map[irow,icol]][np.where( TwoCell[cell_map[irow,icol]] != irow )]) ])
            #             *np.where(FLUX[irow,icol]>0,1,0) 
            #         for irow in range(row)
            #         for icol in range(col)
            #                 ]

            # alphaEface_fast = np.asarray(alphaEface_fast).reshape((row,col))
            # met2_End = tm.time()

            # print('alphaEface_fast:\n',alphaEface_fast,type(alphaEface_fast))
            # # print('For-loop took: %f and comprehension took: %f'%(met1_End-met1_Start,met2_End-met2_Start))
            # # Comment: loop-comprehension(even withour asarray conversion and reshaping) 
            # # for some reason takes more time than for-loop!
            #
            # Verify that all correction factors are within [0,1].
            assert (np.all((alphaEface <= 1 ) & (alphaEface >= 0) ) )\
                    ,"alphaEface are not between 0 and 1!"
            #=========================================;
            # Compute the updated solution and fluxes ;
            #=========================================;
            # Incompressible patch test
            solAvg =  solAvg - self.dt/area * \
                    (1./(0.2))* np.multiply(alphaEface,FLUX).sum(axis=1)
            # solAvg =  solAvg - self.dt/area * \
            #         (1./(0.2*1000))* np.multiply(alphaEface,FLUX).sum(axis=1)
            # incompressible well problem
            # solAvg =  solAvg - self.dt/area * \
            #         (1./(0.2*1000))* np.multiply(alphaEface,FLUX).sum(axis=1) + \
            #         self.dt * (1./(0.2*1000)) * self.wells_avg
            # solAvg =  solAvg - self.dt/area * \
            #         (1./(0.2*1000))* np.multiply(alphaEface,FLUX).sum(axis=1) - \
            #         self.dt * (1./(0.2*1000)) * self.wells_update(solAvg)
            
            # compressible patch test
            # solAvg =  (self.rPAvg_0/self.rPAvg) * solAvg - self.dt/area * \
                    # (1./self.rPAvg)* np.multiply(alphaEface,FLUX).sum(axis=1)
            # print("rpAvg_0/rPAvg\n",self.rPAvg_0/self.rPAvg)
            # solAvg =   solAvg - self.dt/area * \
            #         (1./self.rPAvg)* np.multiply(alphaEface,FLUX).sum(axis=1)
            # solAvg =   solAvg - self.dt/area * \
                    # (1./self.rPAvg_0)* np.multiply(alphaEface,FLUX).sum(axis=1)
            # print('alphaEface',alphaEface[150])
            # last one seems to work the better than others 

            # print('solAvgUpdated at iteration %d is:\n'%numIter,solAvg)
            # print(alphaEface*FLUX)
            FLUX = FLUX *  np.subtract(1.,alphaEface)
            # print('FLUX0 at iteration %d is:\n'%numIter,FLUX_0)
            # print('updatedFLUX at iteration %d is:\n'%numIter,FLUX)

            #=========================;
            # Check stopping criteria ;
            #=========================;
            # Compute new errors:
            # Criterion 1 
            # method 1(the maximum absolute row sum)
            # normFLUX =  np.linalg.norm(FLUX,np.inf)
            # print('normFLUX',normFLUX)
            #
            #method 2   
            normFLUX = np.abs(FLUX).max()
            # print('normFLUX',normFLUX)
        
            # Criterion 2
            # method 1
            # normDelta = np.linalg.norm(np.subtract(FLUX_0,FLUX),np.inf)
            # method 2
            normDelta = np.abs(np.subtract(FLUX_0,FLUX)).max()
            # print('normDelta',normDelta)

            # Check stopping criteria.
            # if (normFLUX < epsFLUX)|(normDelta < epsDelta) :
            # if (normFLUX < epsFLUX)&(normDelta < epsDelta) :
            if normFLUX < epsFLUX :
                flag = 0
                break
            elif normDelta < epsDelta:
                flag = 1;
                break
            elif numIter == maxIter:
                flag = 2;
                break

        # print('Flux-limiter converged in %d iterations'%numIter)
        # print('normFLUX',normFLUX)
        # print('normDelta',normDelta)
        # print('Exit flag is:', flag)


        #==================================;
        # Compute new reconstructed values ;
        #==================================;
        V = self.sol.function_space()
        sol_value = self.sol.vector().array()
        # print('sol_value_old:\n',sol_value)
        solPost = solAvg
        solAvgDG_current = V0.interpolate(self.sol).vector().array() # unlimited solution current step
        # print('sol_value current time (not limited)',solAvgDG_current)
        Diff = solPost-solAvgDG_current
        # print('solPost - solAvgDG_current:\n',Diff)
        sol_cell_map = V.cell_node_map().values 
        # print('sol_cell_map',sol_cell_map)
        row , col = np.shape(sol_cell_map)
        # Add Diff to our nodal solution
        for irow in range(row):
            for icol in range(col):
                sol_value[sol_cell_map[irow,icol]] = sol_value[sol_cell_map[irow,icol]] + Diff[irow]

        # print('sol_value_new:\n',sol_value)
        u_sol = Function(V)
        u_sol.vector().set_local(sol_value)
        # solFinal = interpolate(u_sol,V)

        u_solAvgDG0 = V0.interpolate(u_sol).vector().array() 
        print('limited sol avg after FL:\n min:%f,max:%f'%(u_solAvgDG0.min(),u_solAvgDG0.max()))
        # print('u_sol DG0 avg:\n',solAvgDG0)
         
        # assert (np.all((u_solAvgDG0 <= (self.solMax+1e-7) ) & (u_solAvgDG0 >= (self.solMin-1e-7) ) ) )\
                # ,"constructed avgs are not inbound"
               
        if np.all((u_solAvgDG0 < self.solMax ) & (u_solAvgDG0 > self.solMin)):
            warn('WARNING* averages of constructed sol are not in the range of [solMin,solMax] ')

        return u_sol,numIter

