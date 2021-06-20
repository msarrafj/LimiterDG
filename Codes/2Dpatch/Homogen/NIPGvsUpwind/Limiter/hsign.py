from firedrake import *
import numpy as np
import collections
from ufl.geometry import *

class SignFlip:
    def __init__(self,mesh):
        self.mesh = mesh

    def HSIGN(self):
        Vt = FunctionSpace(self.mesh, "HDiv Trace", 0)
        v = TestFunction(Vt)
        x = SpatialCoordinate(self.mesh)
        np.set_printoptions(precision=9)

        n = FacetNormal(self.mesh)
        normals = Function(Vt)
        # =======================;
        # Extract cell midpoints
        # =======================;
        node_per_cell = self.mesh.coordinates.cell_node_map().arity
        mesh_type = self.mesh.ufl_cell().cellname()
        mesh_val = self.mesh.coordinates.vector().array()
        mesh_cell_map = self.mesh.coordinates.cell_node_map().values
        # File('mesh.pvd').write(self.mesh.coordinates)
        # print('mesh_cell_map',mesh_cell_map)
        row = np.size(mesh_cell_map)
        row,col = self.mesh.cell_set.size , node_per_cell 
        CELL_x = np.zeros((row,col))
        CELL_y = np.zeros((row,col))
        for irow in range(row):
            for icol in range(col):
                CELL_x[irow,icol] = mesh_val[mesh_cell_map[irow,icol]*2] # x-coord
                CELL_y[irow,icol] = mesh_val[mesh_cell_map[irow,icol]*2+1] # y-coord
        cellx_midPt = np.average(CELL_x,axis=1)
        celly_midPt = np.average(CELL_y,axis=1)
        # print('cell_x',CELL_x)
        # =======================;
        # Extract facet midpoints
        # =======================;
        cell_map = Vt.cell_node_map().values
        # np.savetxt('cell_map.out', cell_map,fmt='%10.7g')
        # print('cell_map:\n',cell_map)
        # Fplus = normals('+')*v('+')*dS + normals * v*ds -\
        #         (inner(as_vector([conditional(x[0]<7.5,1,2),0]),n))('+')*v('+') * dS - inner(as_vector([conditional(x[0]<7.5,1,2),0]),n)*v *ds
        #         # (inner(as_vector([0,conditional(x[1]<2.5,1,2)]),n))('+')*v('+') * dS - inner(as_vector([0,conditional(x[1]<2.5,1,2)]),n)*v *ds

        # solve(Fplus == 0, normals)
        # normals_val_x = normals.vector().array()  
        # row , col = np.shape(cell_map)
        # NPLUS = np.zeros((row,col))
        # for irow in range(row):
        #     for icol in range(col):
        #             NPLUS[irow,icol] = normals_val_x[cell_map[irow,icol]]
        # print('Nplus:\n',NPLUS)

        facetx_midPt = np.zeros((row,col))
        facety_midPt = np.zeros((row,col))
        # if element is Quadrilateral
        if mesh_type == "quadrilateral":
            for irow in range(row):
                for icol in range(col):
                    if icol == 0:
                        facetx_midPt[irow,icol] = (CELL_x[irow,0]+CELL_x[irow,1])/2.
                        facety_midPt[irow,icol] = (CELL_y[irow,0]+CELL_y[irow,1])/2.
                    elif icol == 1:
                        facetx_midPt[irow,icol] = (CELL_x[irow,2]+CELL_x[irow,3])/2.
                        facety_midPt[irow,icol] = (CELL_y[irow,2]+CELL_y[irow,3])/2.
                    elif icol == 2:
                        facetx_midPt[irow,icol] = (CELL_x[irow,0]+CELL_x[irow,2])/2.
                        facety_midPt[irow,icol] = (CELL_y[irow,0]+CELL_y[irow,2])/2.
                    else:
                        facetx_midPt[irow,icol] = (CELL_x[irow,1]+CELL_x[irow,3])/2.
                        facety_midPt[irow,icol] = (CELL_y[irow,1]+CELL_y[irow,3])/2.

                
        # if element is Triangular
        if mesh_type == "triangle":
            for irow in range(row):
                for icol in range(col):
                    facetAvgX = 0
                    facetAvgY = 0
                    for k in [x for x in range(col) if x != icol]:
                        facetAvgX = facetAvgX + CELL_x[irow,k]
                        facetAvgY = facetAvgY + CELL_y[irow,k]
                    facetx_midPt[irow,icol] = facetAvgX/2.
                    facety_midPt[irow,icol] = facetAvgY/2.

        # print('facetx_midPt:\n',facetx_midPt)
        # np.savetxt('facetx_midPt.out', facetx_midPt,fmt='%10.7g')
        # print('facety_midPt:\n',facety_midPt)
        # np.savetxt('facety_midPt.out', facety_midPt,fmt='%10.7g')
        # mark boundaries
        left = mesh_val[::2].min()
        right = mesh_val[::2].max()
        bottom = mesh_val[1::2].min()
        top = mesh_val[1::2].max()

        left_facet_idx = (np.where(facetx_midPt == left))
        left_facet = cell_map[left_facet_idx]
        
        right_facet_idx = (np.where(facetx_midPt == right))
        right_facet = cell_map[right_facet_idx]

        bottom_facet_idx = (np.where(facety_midPt == bottom))
        bottom_facet = cell_map[bottom_facet_idx]

        top_facet_idx = (np.where(facety_midPt == top))
        top_facet = cell_map[top_facet_idx]

        # Crack boundaries extracted manually and painfully
        # crack_working.mesh
        # left_crack_facet = np.array([776, 759, 808, 827, 845, 866])
        # right_crack_facet = np.array([630, 649, 704, 746, 789, 811])
        #
        # left_crack_facet = np.array([367,327,267,220,247,289])
        # right_crack_facet = np.array([438,443,446,415,371,348])


        def find_repeat(numbers):
            seen = set()
            for num in numbers:
                if num in seen:
                    return num
                seen.add(num)


        Left_cell = np.array([])
        Right_cell = np.array([])
        Bottom_cell = np.array([])
        Top_cell = np.array([])
        for irow in range(row):
            numbers = CELL_x[irow,:] #each cell
            if find_repeat(numbers) == left:
                Left_cell =  np.append(Left_cell, irow )
            elif find_repeat(numbers) == right:
                Right_cell =  np.append(Right_cell, irow )

            numbers = CELL_y[irow,:] #each cell
            if find_repeat(numbers) == bottom:
                Bottom_cell =  np.append(Bottom_cell, irow )
            elif find_repeat(numbers) == top:
                Top_cell =  np.append(Top_cell, irow )
        #==================================;
        # Figure out two adjacent neighbors;
        # =================================;
        TwoCell = np.ones((Vt.dim(),2)) * -999 
        TwoCell = TwoCell.astype('int32') 
        # all boundary facets have only one neighbor cell, Instead of NA we place -999
        for itrn in range(Vt.dim()):
            neigh = np.concatenate(np.where(cell_map == itrn)[::2] )
            for itrnj in range(len(neigh)):
                TwoCell[itrn,itrnj] = neigh[itrnj]
        # print('TwoCell',TwoCell)
        # Mark boundary facets:
        boundary_facet = np.concatenate(np.where(TwoCell<0)[::2])
        # marking left/right/top/bottom boundaries
        # old implementation of marking boundaries
        # it has a bug for 'right' and 'left' mesh
        # only works for 'crossed'
        # identify index where first column is Left_cell and last column is -999
        # left_facet= np.array([])
        # for cellidx in Left_cell:
        #     print('cellidx',cellidx)
        #     fct_index = np.concatenate( np.where( (TwoCell[:,-1] == -999)&(TwoCell[:,0] == cellidx)) )
        #     left_facet = np.append(left_facet,fct_index)
        # print('left_facet',left_facet)

        # right_facet= np.array([])
        # for cellidx in Right_cell:
        #     fct_index = np.concatenate( np.where( (TwoCell[:,-1] == -999)&(TwoCell[:,0] == cellidx)) )
        #     right_facet = np.append(right_facet,fct_index)
        # print('right_face',right_facet)

        # top_facet= np.array([])
        # for cellidx in Top_cell:
        #     fct_index = np.concatenate( np.where( (TwoCell[:,-1] == -999)&(TwoCell[:,0] == cellidx)) )
        #     top_facet = np.append(top_facet,fct_index)
        # print('top_facet',top_facet)

        # bottom_facet= np.array([])
        # for cellidx in Bottom_cell:
        #     fct_index = np.concatenate( np.where( (TwoCell[:,-1] == -999)&(TwoCell[:,0] == cellidx)) )
        #     bottom_facet = np.append(bottom_facet,fct_index)
        # print('bottom_facet',bottom_facet)

        #======================================;
        #determine + - edges at each cell map
        #======================================;
        Fplus = normals('+')*v('+')*dS + normals * v*ds -\
                (inner(Constant((1,0)),n) )('+')*v('+') * dS - inner(Constant((1,0)),n)*v *ds
        solve(Fplus == 0, normals)
        normals_val_x = normals.vector().array()  

        Fminus = normals('+')*v('+')*dS + normals * v*ds -\
                (inner(Constant((0,1)),n) )('+')*v('+') * dS - inner(Constant((0,1)),n)*v *ds
        solve(Fminus == 0, normals)
        normals_val_y = normals.vector().array()  

        IdxVertic = np.concatenate(np.where(normals_val_x == 0) )
        # IdxVertic = np.concatenate(np.where((normals_val_x >-1e-5)|(normals_val_x<1e-5) ))
        # IdxVertic = np.concatenate(np.where((normals_val_x >-10)|(normals_val_x<10) ))

        row , col = np.shape(cell_map)
        NPLUS = np.zeros((row,col))
        for irow in range(row):
            for icol in range(col):
                if cell_map[irow,icol] in IdxVertic:
                    NPLUS[irow,icol] = normals_val_y[cell_map[irow,icol]]

                else:
                    NPLUS[irow,icol] = normals_val_x[cell_map[irow,icol]]

        # print("NPLUS",NPLUS[151,:])
        # np.savetxt('NPLUS',NPLUS)

        Hsign = np.zeros((row,col))
        # We have three levels: single-valued Facet; two cells that are adjacent to
        # Facet; We need to compare midpoint of two cells adjacent to facet
        # facet number is calculated from cell_map[irow,icol] and is one value
        # Twocell[cell_map[irow,icol]] gives two cells adjacent to a facet
        # print(arr[np.where(arr != 2)])] gives all elements of array except element 2
        # Twocell[np.where( Twocell[cell_map[irow,icol] != irow )]]
        # gives all elements of Two cell except the irow of current cell. In other words,
        # it is giving the ID cell adjacent to facet by not the current cell.
        # cellx_midP[irow] gives the midpoint of cell irow; 
        # cellx_midP[Twocell[cell_map[irow,icol]]    ] gives the
        # cellx_midP[ Twocell[cell_map[irow,icol]]    ] gives the
        #
        def stepfunc(x):
            if x>0:
                return int(1)
            elif x<0: 
                return int(-1)
            else:
                return int(0)

        for irow in range(row):
            for icol in range(col):
                facet = cell_map[irow,icol]
                # if facet in np.concatenate((left_facet,bottom_facet,right_crack_facet)): #if facet belongs to bottom/left boundary
                if facet in np.concatenate((left_facet,bottom_facet)): #if facet belongs to bottom/left boundary
                    Hsign[irow,icol] = -1 *stepfunc(NPLUS[irow,icol])
                # elif facet in np.concatenate((top_facet,right_facet,left_crack_facet)): #if facet belongs to top/right boundary
                elif facet in np.concatenate((top_facet,right_facet)): #if facet belongs to top/right boundary
                    Hsign[irow,icol] =  1 *stepfunc(NPLUS[irow,icol])
                else: # Inner facets
                    if facet in IdxVertic: # first we address vertical facets
                        a = celly_midPt[irow] # midpoint of current cell
                        b0 = TwoCell[facet] # cell includes irow and the opposite cell
                        b1 = int(b0[np.where( b0 != irow )]) # Get the cell number; includes only oppoiste cell
                        b = celly_midPt[b1] # midpoint of opposite cell
                        # Now we compare a and b:
                        # if NPLUS is positive ==> if a>b Hsign= -1  else a<b = +1
                        # if NPLUS is negative ==> if a>b Hsign= +1  else a<b = -1
                        # this is equivalent to -1 * stepfunc(NPLUS*(a-b) )
                        Hsign[irow,icol] = -1 * stepfunc(NPLUS[irow,icol]*(a-b))


                    else: # Now we address all inner facets except vertical
                        a = cellx_midPt[irow] # midpoint of current cell
                        b0 = TwoCell[facet] # cell includes irow and the opposite cell
                        b1 = int(b0[np.where( b0 != irow )]) # includes only oppoiste cell
                        b = cellx_midPt[b1] # midpoint of opposite cell
                        # Now we compare a and b:
                        # if NPLUS is positive ==> if a>b Hsign= -1  else a<b = +1
                        # if NPLUS is negative ==> if a>b Hsign= +1  else a<b = -1
                        # this is equivalent to -1 * stepfunc(NPLUS*(a-b) )
                        Hsign[irow,icol]  = -1 * stepfunc(NPLUS[irow,icol]*(a-b))

        return Hsign,cell_map,TwoCell,boundary_facet


# mesh = RectangleMesh(2, 2 ,1,1,diagonal='crossed')
# signing = SignFlip(mesh)
# Hsign = signing.HSIGN()

# print('Hsign', Hsign)
