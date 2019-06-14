# msu MeshFileDomain.py

import sys
#import numpy as np
#from proteus.Profiling import logEvent
#from proteus import MeshTools
from proteus.mprans import BoundaryConditions as bc
from proteus import Domain
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus import default_p
from proteus.msu import msu_bc as msu_bc

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

    def bc_zone_define( self, name=None, meshtag=None, condition={'type':'undefined'}, custom=None ):

        if name == None:                raise RuntimeError('name=None')
        if name in self.bc_zone.keys(): raise RuntimeError('name already exists: {0}'.format(name) )

        print('++++++ bc defined:', name, condition['type'], custom )

        self.bc_zone[name] = { 'name':      name, 
                               'meshtag':   meshtag, 
                               'condition': condition,
                               'custom':    custom,
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
            custom    = zone['custom']
            typename  = condition['type']
            
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

            used = False;

            print('====== bc set:', zone['name'], condition['type'], custom )

            if custom and not used: 
               custom( bc, condition ); used = True;
               continue

            if 'freeSlip' in typename and not used: 
               msu_bc.freeSlip( bc ); used = True;
               continue

            if 'noSlip' in typename and not used: 
               msu_bc.noSlip( bc ); used = True;
               continue

            if 'processorBoundary' in typename and not used: 
               msu_bc.nonMaterial( bc ); used = True
               continue

            if 'tank' in typename and not used:  
               if 'AxisOriented' in typename:  
                 for i in range(3): normal[i] = int(round(normal[i])) # nint
                 msu_bc.tank( bc, normal ); used = True
                 print( zone['meshtag'], self.meshtag_bcIdx[zone['meshtag']], zone['name'], normal )
               else:
                 raise RuntimeError("curved Tank not supported yet")
               continue
  
            if 'openAir' in typename and not used:  
               if 'AxisOriented' in typename:  
                 for i in range(3): normal[i] = int(round(normal[i])) # nint
                 msu_bc.openAir( bc, normal ); used = True
               else:
                 raise RuntimeError("curved OpenAir not supported yet")
               continue

            if 'velocityInlet' in typename and not used: 
               if 'AxisOriented' in typename:  
                 for i in range(3): normal[i] = int(round(normal[i])) # nint
               if 'twoPhase' in typename:  
                 msu_bc.velocityInlet_rans2p( bc, condition, normal ); used = True
               elif 'rans2p' in typename:  
                 msu_bc.velocityInlet_rans2p( bc, condition, normal ); used = True
               elif 'rans3p' in typename:  
                 msu_bc.velocityInlet_rans3p( bc, condition, normal ); used = True
               else:
                 msu_bc.velocityInlet( bc, condition, normal ); used = True
               continue

            if 'outflow' in typename and not used: 
               if 'twoPhase' in typename:  
                 msu_bc.outflow_rans2p( bc, condition ); used = True
               elif 'rans2p' in typename:  
                 msu_bc.outflow_rans2p( bc, condition ); used = True
               elif 'rans3p' in typename:  
                 msu_bc.outflow_rans3p( bc, condition ); used = True
               else:
                 msu_bc.outflow( bc, condition ); used =  True
               continue

            if 'interior' in typename and not used: 
               msu_bc.interior( bc, condition ); used = True
               continue

            if 'open' in typename and not used: 
               msu_bc.open( bc, condition ); used = True
               continue

        if not used:
          print('typename =', typename )
          sys.exit('bc not handled in msu/MeshFileDomain.py')
