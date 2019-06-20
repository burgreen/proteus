# msu MeshFileDomain.py

import sys
#import numpy as np
#from proteus.Profiling import logEvent
#from proteus import MeshTools
from proteus.mprans import BoundaryConditions as bc
from proteus import Domain
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus import default_p

class MeshFileDomain(Domain.D_base):
    """
    domains from ordinary mesh files
    """
    def __init__( self, filename, dim=3 ):

        Domain.D_base.__init__(self,dim)
        self.filename  = filename
        self.dim       = dim
        self.bc_zone   = {}
        self.bc_zones  = []
        self.ele_zone  = {}
        self.ele_zones = []
        self.meshtag_bcIdx = {}
        self.interiorEdgeTags = []

    def get_bcIdx( self, flag ):
        print( 'flag=',flag )
        return self.meshtag_bcIdx[flag]

    def bc_zone_define( self, name=None, meshtag=None, condition={'method':None} ):

        if name == None:                raise RuntimeError('name=None')
        if name in self.bc_zone.keys(): raise RuntimeError('name already exists: {0}'.format(name) )
        if condition['method'] == None: raise RuntimeError('bc method is not assigned')

        print('++++++ bc defined:', name, condition['method'] )

        self.bc_zone[name] = { 'name':      name, 
                               'meshtag':   meshtag, 
                               'condition': condition,
                               'bc_rans':   bc.BC_RANS(nd=self.dim,name=name) }

        self.bc_zones.append( self.bc_zone[name] )

        self.meshtag_bcIdx[meshtag] = len(self.bc)

        self.bc.append( self.bc_zone[name]['bc_rans'] )

    def ele_zone_define( self, name=None, meshtag=None, condition={type:'undefined'} ):

        if name == None:                 raise RuntimeError('name=None')
        if name in self.ele_zone.keys(): raise RuntimeError('name already exists: {0}'.format(name) )
        
        self.ele_zone[name] = { 'name':name, 'meshtag':meshtag, 'condition':condition }

        self.ele_zones.append( self.ele_zone[name] )

    def set_interiorEdgeTags( self, tagList ):

        self.interiorEdgeTags = np.array(tagList)

    def print_bc_attrs( self, mesh ):

        print( 'mesh.num_attr:', mesh.num_attrF )
        for i in range(mesh.num_attrF): print( 'attrF:', mesh.attrF[i]  )

    def set_bc_definitions( self, mesh ):

        # called after mesh is read and processed in NumericalSolution.py

        for zone in self.bc_zones: 
            
            bc        = zone['bc_rans']
            condition = zone['condition']
            meshtag   = zone['meshtag']
            method    = condition['method']
            
            normal = [0,0,0]
            used = False
            for i in range(mesh.num_attrF): 
               if mesh.attrF[i] == meshtag:
                  used = True
                  normal[0] = mesh.nx_mean[i]
                  normal[1] = mesh.ny_mean[i]
                  normal[2] = mesh.nz_mean[i]
            if not used: 
               #print( 'mesh.num_attr:', mesh.num_attrF )
               #for i in range(mesh.num_attrF): print( 'attrF:', mesh.attrF[i]  )
               #raise RuntimeError("meshtag not found: " + str(meshtag))
               continue

            condition['_normal'] = normal

            method( bc, condition )
