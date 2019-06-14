# msu/msu_bc.py

import sys
#import numpy as np
#from proteus.Profiling import logEvent
#from proteus import MeshTools
#from proteus.mprans import BoundaryConditions as bc
#from proteus import Domain
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus import default_p

def tank( bc, normal ):
    '''
    Sets tank wall bc for sloshing type problems
    '''
    bc.u_stress.uOfXT = 0.
    bc.v_stress.uOfXT = 0.
    bc.w_stress.uOfXT = 0.
    bc.hx_dirichlet.resetBC()
    bc.hy_dirichlet.resetBC()
    bc.hz_dirichlet.resetBC()
    if normal[0] == 1 or normal[0] == -1:
       bc.hx_dirichlet.setConstantBC(0.)
       bc.u_stress.uOfXT = None
    elif normal[1] == 1 or normal[1] == -1:
       bc.hy_dirichlet.setConstantBC(0.)
       bc.v_stress.uOfXT = None
    elif len(normal) > 2 and (normal[2] == 1 or normal[2] == -1):
       bc.hz_dirichlet.setConstantBC(0.)
       bc.w_stress.uOfXT = None

def openAir( bc, normal, vof_air=1. ):
    """
    Sets atmosphere boundary conditions (water can come out)
    (!) pressure dirichlet set to 0 for this BC
    Parameters
    ----------
    normal: 
    normal of the boundary. Optional if orientation was already
    passed when creating the BC_RANS class instance.
    vof_air: Optional[float]
    VOF value of air (default is 1.)
    """
    bc.BC_type = 'OpenAir'
    bc.reset()
    bc.p_dirichlet.setConstantBC(0.)
    bc.pInc_dirichlet.setConstantBC(0.)
    bc.pInit_dirichlet.setConstantBC(0.)
    bc.u_dirichlet.setConstantBC(0.)
    bc.v_dirichlet.setConstantBC(0.)
    bc.w_dirichlet.setConstantBC(0.)
    bc.us_dirichlet.setConstantBC(0.)
    bc.vs_dirichlet.setConstantBC(0.)
    bc.ws_dirichlet.setConstantBC(0.)
    bc.vof_dirichlet.setConstantBC(vof_air)  # air
    bc.vos_dirichlet.setConstantBC(0.)  # air
    if normal[0] == 1. or normal[0] == -1.:
        bc.u_diffusive.setConstantBC(0.)
        bc.us_diffusive.setConstantBC(0.)
    if normal[1] == 1. or normal[1] == -1.:
        bc.v_diffusive.setConstantBC(0.)
        bc.vs_diffusive.setConstantBC(0.)
    if normal[2] == 1. or normal[2] == -1.:
        bc.w_diffusive.setConstantBC(0.)
        bc.ws_diffusive.setConstantBC(0.)
    bc.k_diffusive.setConstantBC(0.)
    bc.dissipation_diffusive.setConstantBC(0.)

def nonMaterial( bc ):
    """
    Sets non-material boundary conditions (diffusive flux and advective vof to 0.).
    """
    bc.reset()
    bc.BC_type = 'NonMaterial'
    bc.vof_advective.setConstantBC(0.)
    bc.vos_advective.setConstantBC(0.)
    bc.u_diffusive.setConstantBC(0.)
    bc.v_diffusive.setConstantBC(0.)
    bc.w_diffusive.setConstantBC(0.)
    bc.us_diffusive.setConstantBC(0.)
    bc.vs_diffusive.setConstantBC(0.)
    bc.ws_diffusive.setConstantBC(0.)
    bc.k_diffusive.setConstantBC(0.)
    bc.dissipation_diffusive.setConstantBC(0.)
    bc.pInc_diffusive.setConstantBC(0.)

def outflow( bc, condition ):
    """
    Sets outflow bc
    """
    bc.reset()

    def constant(c):
        def function(x,t): return c
        return function

    bc.p_dirichlet.setConstantBC(0.)
    bc.u_diffusive.setConstantBC(0.)
    bc.v_diffusive.setConstantBC(0.)
    bc.w_diffusive.setConstantBC(0.)
    #bc.u_advective.setConstantBC(0.)
    #bc.v_advective.setConstantBC(0.)
    #bc.w_advective.setConstantBC(0.)

    bc.k_dirichlet.uOfXT = constant(0.)
    bc.k_advective.uOfXT = None
    bc.k_diffusive.uOfXT = constant(0.)

    bc.dissipation_dirichlet.uOfXT = constant(0.)
    bc.dissipation_advective.uOfXT = None
    bc.dissipation_diffusive.uOfXT = constant(0.)

    #Vmag = 0.
    #if 'Vmag' in condition: Vmag = condition['Vmag']
    #bc.p_advective.setConstantBC(Vmag)

    if 'p'    in condition: bc.p_dirichlet.setConstantBC(condition['p'])
    if 'p_xt' in condition: bc.p_dirichlet.uOfXT = condition['p_xt']
    #set pInit and Pinc??

def velocityInlet( bc, condition, normal ):
    """
    Sets velocity inlet bc
    """
    bc.reset()

    Vmag = 0.
    if 'Vmag' in condition: Vmag = condition['Vmag']

    turb_k = 0.
    turb_e = 0.
    if 'turb_k' in condition: turb_k = condition['turb_k']
    if 'turb_e' in condition: turb_e = condition['turb_e']

    def constant(c):
        def function(x,t): return c
        return function

    bc.u_dirichlet.setConstantBC(Vmag*normal[0])
    bc.v_dirichlet.setConstantBC(Vmag*normal[1])
    bc.w_dirichlet.setConstantBC(Vmag*normal[2])

    if 'u_xt' in condition: bc.u_dirichlet.uOfXT = condition['u_xt']
    if 'v_xt' in condition: bc.v_dirichlet.uOfXT = condition['v_xt']
    if 'w_xt' in condition: bc.w_dirichlet.uOfXT = condition['w_xt']

    bc.k_dirichlet.uOfXT = constant(turb_k)
    #bc.k_advective.
    bc.k_diffusive.uOfXT = constant(0.)

    bc.dissipation_dirichlet.uOfXT = constant(turb_e)
    #bc.dissipation_advective.
    bc.dissipation_diffusive.uOfXT = constant(0.)

    bc.p_advective.setConstantBC(-Vmag)

    if 'flux_xt' in condition: bc.p_advective.uOfXT = condition['flux_xt']
    
    #set pInit and Pinc??

def velocityInlet_rans2p( bc, condition, normal ):
    """
    Sets velocity inlet bc for twoPhase flows
    """
    bc.reset()

    Vmag = 0.
    if 'Vmag' in condition: Vmag = condition['Vmag']
    Vmag_air   = Vmag*0
    Vmag_water = Vmag

    turb_k = 0.
    turb_e = 0.
    if 'turb_k' in condition: turb_k = condition['turb_k']
    if 'turb_e' in condition: turb_k = condition['turb_e']

    sdf   = condition['sdf']
    fluid = condition['fluid']

    phi_air   = fluid['air']['phi']
    phi_water = fluid['water']['phi']

    def vel_dirichlet(i):
            def dirichlet(x,t):
                phi = sdf(x)
                #if phi <= 0.:
                #    H = 0.0
                #elif 0 < phi <= smoothing:
                #    H = smoothedHeaviside(smoothing / 2., phi - smoothing / 2.)
                #else:
                #    H = 1.0
                H = 1.
                if phi <= 0.: H = 0.
                vel = H * Vmag_air*normal[i] + (1-H) * Vmag_water*normal[i]
                return vel
            return dirichlet

    def vof_dirichlet(x,t):
            phi = sdf(x)
            #if phi >= smoothing:
            #    H = 1.
            #elif smoothing > 0 and -smoothing < phi < smoothing:
            #    H = smoothedHeaviside(smoothing, phi)
            #elif phi <= -smoothing:
            #    H = 0.
            H = 1.
            if phi <= 0.: H = 0.
            vof = H * phi_air + (1-H) * phi_water
            return vof

    def p_advective(x,t):
            phi = sdf(x)
            #if phi <= 0.:
            #    H = 0.0
            #elif 0 < phi <= smoothing:
            #    H = smoothedHeaviside(smoothing / 2., phi - smoothing / 2.)
            #else:
            #    H = 1.0
            H = 1.
            if phi <= 0.: H = 0.
            u = H * Vmag_air + (1-H) * Vmag_water
            # This is the normal velocity, based on the inwards boundary
            # orientation -b_or
            #u_p = np.sum(u * np.abs(b_or))
            #return -u_p
            return -u

    def constant(c):
        def function(x,t): return c
        return function

    bc.u_dirichlet.uOfXT = vel_dirichlet(0)
    bc.u_advective.uOfXT = None
    bc.u_diffusive.uOfXT = None

    bc.v_dirichlet.uOfXT = vel_dirichlet(1)
    bc.w_advective.uOfXT = None
    bc.w_diffusive.uOfXT = None

    bc.w_dirichlet.uOfXT = vel_dirichlet(2)
    bc.w_advective.uOfXT = None
    bc.w_diffusive.uOfXT = None

    if 'u_xt' in condition: bc.u_dirichlet.uOfXT = condition['u_xt']
    if 'v_xt' in condition: bc.v_dirichlet.uOfXT = condition['v_xt']
    if 'w_xt' in condition: bc.w_dirichlet.uOfXT = condition['w_xt']

    bc.k_dirichlet.uOfXT = constant(turb_k)
    #bc.k_advective.
    #bc.k_diffusive.

    bc.dissipation_dirichlet.uOfXT = constant(turb_e)
    #bc.dissipation_advective.
    #bc.dissipation_diffusive.

    bc.vof_dirichlet.uOfXT = vof_dirichlet
    bc.vof_advective.uOfXT = constant(0.)
    ##.vof_diffusive.uOfXT = does not exist

    bc.p_dirichlet.uOfXT = None
    bc.p_advective.uOfXT = p_advective
    ##.p_diffusive.uOfXT = does not exist

    if 'flux_xt' in condition: bc.p_advective.uOfXT = condition['flux_xt']
    
def outflow_rans2p( bc, condition ):
    """
    Sets outflow bc for twoPhase flows
    """
    bc.reset()

    sdf        = condition['sdf']
    fluid      = condition['fluid']
    gravity    = condition['gravity']
    waterLevel = condition['waterLevel']
    smoothing  = 0.05 
    fac        = 0.1 

    rho_air   = fluid['air']['rho']
    rho_water = fluid['water']['rho']
    mu_air    = fluid['air']['mu']
    mu_water  = fluid['water']['mu']
    phi_air   = fluid['air']['phi']
    phi_water = fluid['water']['phi']

    def p_dirichlet(x,t):  # assumes a constant fixed fluid height at the outlet
        g_component = gravity[2]
        p_air   = rho_air    * g_component * ( default_p.L[2] - waterLevel )
        p_water = rho_water  * g_component * ( waterLevel - x[2] )
        p_hydrostatic = p_air
        #if sdf(x) < 0: p_hydrostatic = p_water
        if x[2] < waterLevel: p_hydrostatic = p_water
        return -p_hydrostatic

    def u_diffusive(x,t):
        g = p_dirichlet(x,t)
        phi = sdf(x)
        if phi >= smoothing: H = 1.
        elif -smoothing < phi < smoothing: H = smoothedHeaviside(smoothing, phi)
        elif phi <= -smoothing: H = 0.
        #H = 1.
        #if phi <= 0.: H = 0.
        diff =  H*(mu_air)*g + (1-H)*(mu_water)*g
        return fac*diff

    def vof_dirichlet(x,t):
        phi = sdf(x)
        if phi >= smoothing: H = 1.
        elif -smoothing < phi < smoothing: H = smoothedHeaviside(smoothing, phi)
        elif phi <= -smoothing: H = 0.
        #H = 1.
        #if phi <= 0.: H = 0.
        return H * phi_air + (1-H) * phi_water

    def constant(c):
        def function(x,t): return c
        return function

    bc.u_dirichlet.uOfXT = constant(0.)
    bc.u_advective.uOfXT = None
    bc.u_diffusive.uOfXT = u_diffusive # some level is needed for stability

    bc.v_dirichlet.uOfXT = constant(0.)
    bc.v_advective.uOfXT = None
    bc.v_diffusive.uOfXT = constant(0.)

    bc.w_dirichlet.uOfXT = constant(0.)
    bc.w_advective.uOfXT = None
    bc.w_diffusive.uOfXT = constant(0.)

    bc.k_dirichlet.uOfXT = constant(0.)
    bc.k_advective.uOfXT = None
    bc.k_diffusive.uOfXT = constant(0.)

    bc.dissipation_dirichlet.uOfXT = constant(0.)
    bc.dissipation_advective.uOfXT = None
    bc.dissipation_diffusive.uOfXT = constant(0.)

    bc.vof_dirichlet.uOfXT = vof_dirichlet
    bc.vof_advective.uOfXT = None
    ##.vof_diffusive.uOfXT = does not exist

    bc.p_dirichlet.uOfXT = p_dirichlet
    bc.p_advective.uOfXT = None
    ##.p_diffusive.uOfXT = does not exist

    '''
    if U is not None:
            def get_inlet_ux_dirichlet(i):
                def ux_dirichlet(x, t):
                    phi = x[vert_axis] - seaLevel
                    if phi <= 0.:
                        H = 0.0
                    elif 0 < phi <= smoothing:
                        H = smoothedHeaviside(smoothing / 2., phi - smoothing / 2.)
                    else:
                        H = 1.0
                    return H * Uwind[i] + (1 - H) * U[i]
                return ux_dirichlet

            if Uwind is None:
                Uwind = np.zeros(3)
            U = np.array(U)
            Uwind = np.array(Uwind)
            self.u_dirichlet.uOfXT = get_inlet_ux_dirichlet(0)
            self.v_dirichlet.uOfXT = get_inlet_ux_dirichlet(1)
            self.w_dirichlet.uOfXT = get_inlet_ux_dirichlet(2)
            self.u_diffusive.resetBC()
    '''

def freeSlip( bc ):
    '''
    Sets free slip conditions at the boundary
    '''
    bc.reset()
    bc.BC_type = 'freeSlip'

    bc.u_dirichlet.uOfXT = None
    bc.u_advective.uOfXT = None
    bc.u_diffusive.setConstantBC(0.) 

    bc.v_dirichlet.uOfXT = None
    bc.v_advective.uOfXT = None
    bc.v_diffusive.setConstantBC(0.)

    bc.w_dirichlet.uOfXT = None
    bc.w_advective.uOfXT = None
    bc.w_diffusive.setConstantBC(0.)

    bc.vof_dirichlet.uOfXT = None
    bc.vof_advective.setConstantBC(0.)
    ##.vof_diffusive.uOfXT = does not exist

    bc.k_dirichlet.setConstantBC(0.) 
    #bc.k_advective.
    bc.k_diffusive.setConstantBC(0.)

    #bc.dissipation_dirichlet.
    #bc.dissipation_advective.
    bc.dissipation_diffusive.setConstantBC(0.)  
    
    bc.p_dirichlet.uOfXT = None
    bc.p_advective.setConstantBC(0.)
    #bc.p_diffusive.uOfXT = None

    bc.pInc_dirichlet.uOfXT = None
    bc.pInc_advective.setConstantBC(0.)  
    bc.pInc_diffusive.setConstantBC(0.)

    bc.pInit_dirichlet.uOfXT = None
    bc.pInit_advective.setConstantBC(0.)
    bc.pInit_diffusive.setConstantBC(0.)

    # dirichlet
    # advective        
    bc.us_advective.setConstantBC(0.)
    bc.vs_advective.setConstantBC(0.)
    bc.ws_advective.setConstantBC(0.)
    bc.vos_advective.setConstantBC(0.)
    # diffusive
    bc.us_diffusive.setConstantBC(0.)
    bc.vs_diffusive.setConstantBC(0.)
    bc.ws_diffusive.setConstantBC(0.)

def velocityInlet_rans3p( bc, condition, normal ):
    """
    Sets velocity inlet bc for rans3p flows
    """
    bc.reset()

    Vmag = 0.
    if 'Vmag' in condition: Vmag = condition['Vmag']
    Vmag_air   = Vmag*0
    Vmag_water = Vmag

    sdf   = condition['sdf']
    fluid = condition['fluid']
    smoothing = condition['smoothing']

    phi_air   = fluid['air']['phi']
    phi_water = fluid['water']['phi']

    def vel_dirichlet(i):
            def dirichlet(x,t):
                phi = sdf(x)
                if phi >= smoothing: H = 1.
                elif -smoothing < phi < smoothing: H = smoothedHeaviside(smoothing, phi)
                elif phi <= -smoothing: H = 0.
                #H = 1.
                #if phi <= 0.: H = 0.
                vel = H * Vmag_air*normal[i] + (1-H) * Vmag_water*normal[i]
                return vel
            return dirichlet

    def vof_dirichlet(x,t):
            phi = sdf(x)
            if phi >= smoothing: H = 1.
            elif -smoothing < phi < smoothing: H = smoothedHeaviside(smoothing, phi)
            elif phi <= -smoothing: H = 0.
            #H = 1.
            #if phi <= 0.: H = 0.
            vof = H * phi_air + (1-H) * phi_water
            return vof

    def p_advective(x,t):
            phi = sdf(x)
            if phi >= smoothing: H = 1.
            elif -smoothing < phi < smoothing: H = smoothedHeaviside(smoothing, phi)
            elif phi <= -smoothing: H = 0.
            #H = 1.
            #if phi <= 0.: H = 0.
            u = H * Vmag_air + (1-H) * Vmag_water
            # This is the normal velocity, based on the inwards boundary
            # orientation -b_or
            #u_p = np.sum(u * np.abs(b_or))
            #return -u_p
            return -u

    def constant(c):
            def function(x,t): return c
            return function

    bc.u_dirichlet.uOfXT = vel_dirichlet(0)
    bc.u_advective.uOfXT = None
    bc.u_diffusive.uOfXT = constant(0.)

    bc.v_dirichlet.uOfXT = vel_dirichlet(1)
    bc.v_advective.uOfXT = None
    bc.v_diffusive.uOfXT = constant(0.)

    bc.w_dirichlet.uOfXT = vel_dirichlet(2)
    bc.w_advective.uOfXT = None
    bc.w_diffusive.uOfXT = constant(0.)

    if 'u_xt' in condition: bc.u_dirichlet.uOfXT = condition['u_xt']
    if 'v_xt' in condition: bc.v_dirichlet.uOfXT = condition['v_xt']
    if 'w_xt' in condition: bc.w_dirichlet.uOfXT = condition['w_xt']

    # which dirichlet??
    bc.vof_dirichlet.uOfXT = vof_dirichlet
    ##.vof_dirichlet.uOfXT = None
    bc.vof_advective.uOfXT = constant(0.)
    ##.vof_diffusive.uOfXT = does not exist

    bc.p_dirichlet.uOfXT = None
    bc.p_advective.uOfXT = constant(0.)
    ##.p_diffusive.uOfXT = does not exist

    if 'flux_xt' in condition: bc.p_advective.uOfXT = condition['flux_xt']

    bc.pInc_dirichlet.uOfXT = None
    bc.pInc_advective.uOfXT = p_advective
    bc.pInc_diffusive.uOfXT = constant(0.)

    bc.pInit_dirichlet.uOfXT = None
    bc.pInit_advective.uOfXT = constant(0.)
    bc.pInit_diffusive.uOfXT = constant(0.)
    
def outflow_rans3p( bc, condition ):
    """
    Sets outflow bc for twoPhase flows
    """
    bc.reset()

    sdf        = condition['sdf']
    fluid      = condition['fluid']
    gravity    = condition['gravity']
    waterLevel = condition['waterLevel']
    smoothing  = condition['smoothing']

    rho_air   = fluid['air']['rho']
    rho_water = fluid['water']['rho']
    mu_air    = fluid['air']['mu']
    mu_water  = fluid['water']['mu']
    phi_air   = fluid['air']['phi']
    phi_water = fluid['water']['phi']

    def p_dirichlet(x,t):  # assumes a constant fixed fluid height at the outlet
        g_component = gravity[2]
        p_air   = rho_air    * g_component * ( default_p.L[2] - waterLevel )
        p_water = rho_water  * g_component * ( waterLevel - x[2] )
        p_hydrostatic = p_air
        #if phi >= smoothing: H = 1.
        #elif -smoothing < phi < smoothing: H = smoothedHeaviside(smoothing, phi)
        #elif phi <= -smoothing: H = 0.
        if sdf(x) < 0: p_hydrostatic += p_water
        return -p_hydrostatic

    def u_diffusive(x,t):
        g = p_dirichlet(x,t)
        phi = sdf(x)
        if phi >= smoothing: H = 1.
        elif -smoothing < phi < smoothing: H = smoothedHeaviside(smoothing, phi)
        elif phi <= -smoothing: H = 0.
        #H = 1.
        #if phi <= 0.: H = 0.
        return H*(mu_air)*g + (1-H)*(mu_water)*g

    def vof_dirichlet(x,t):
        phi = sdf(x)
        if phi >= smoothing: H = 1.
        elif -smoothing < phi < smoothing: H = smoothedHeaviside(smoothing, phi)
        elif phi <= -smoothing: H = 0.
        #H = 1.
        #if phi <= 0.: H = 0.
        return H * phi_air + (1-H) * phi_water

    def constant(c):
        def function(x,t): return c
        return function

    bc.u_dirichlet.uOfXT = constant(0.)
    bc.u_advective.uOfXT = None
    bc.u_diffusive.uOfXT = u_diffusive

    bc.v_dirichlet.uOfXT = constant(0.)
    bc.v_advective.uOfXT = None
    bc.v_diffusive.uOfXT = constant(0.)

    bc.w_dirichlet.uOfXT = constant(0.)
    bc.w_advective.uOfXT = None
    bc.w_diffusive.uOfXT = constant(0.)

    bc.vof_dirichlet.uOfXT = vof_dirichlet
    bc.vof_advective.uOfXT = None
    ##.vof_diffusive.uOfXT = does not exist

    bc.p_dirichlet.uOfXT = p_dirichlet
    bc.p_advective.uOfXT = None
    ##.p_diffusive.uOfXT = does not exist

    bc.pInc_dirichlet.uOfXT = constant(0.)
    bc.pInc_advective.uOfXT = None
    bc.pInc_diffusive.uOfXT = None
    
    bc.pInit_dirichlet.uOfXT = p_dirichlet
    bc.pInit_advective.uOfXT = None
    bc.pInit_diffusive.uOfXT = None

def interior( bc, condition ):
    """
    Sets interior bc for rans3p flows
    """
    bc.reset()

    def constant(c):
        def function(x,t): return c
        return function

    bc.p_dirichlet.uOfXT = None

    bc.u_diffusive.setConstantBC(0.)
    bc.v_diffusive.setConstantBC(0.)
    bc.w_diffusive.setConstantBC(0.)
    bc.pInc_diffusive.setConstantBC(0.)
    bc.us_diffusive.setConstantBC(0.)
    bc.vs_diffusive.setConstantBC(0.)
    bc.ws_diffusive.setConstantBC(0.)
    bc.k_diffusive.setConstantBC(0.)
    bc.dissipation_diffusive.setConstantBC(0.)

    bc.vof_advective.setConstantBC(0.)
    bc.vos_advective.setConstantBC(0.)

def open( bc, condition ):
    """
    Sets open air bc for rans3p flows
    """
    bc.reset()

    def constant(c):
        def function(x,t): return c
        return function

    bc.u_dirichlet.uOfXT = None
    bc.u_advective.uOfXT = None
    bc.u_diffusive.uOfXT = constant(0.)

    bc.v_dirichlet.uOfXT = None
    bc.v_advective.uOfXT = None
    bc.v_diffusive.uOfXT = constant(0.)

    bc.w_dirichlet.uOfXT = None
    bc.w_advective.uOfXT = None
    bc.w_diffusive.uOfXT = constant(0.)

    bc.vof_dirichlet.uOfXT = constant(0.)
    bc.vof_advective.uOfXT = None
    ##.vof_diffusive.uOfXT = does not exist

    bc.p_dirichlet.uOfXT = constant(0.)
    bc.p_advective.uOfXT = None
    ##.p_diffusive.uOfXT = does not exist

    bc.pInc_dirichlet.uOfXT = constant(0.)
    bc.pInc_advective.uOfXT = None
    bc.pInc_diffusive.uOfXT = None
    
    bc.pInit_dirichlet.uOfXT = constant(0.)
    bc.pInit_advective.uOfXT = None
    bc.pInit_diffusive.uOfXT = None

def noSlip( bc ):
    """
    Sets no slip conditions at the boundary
    """
    bc.reset()
    bc.BC_type = 'NoSlip'

    bc.u_dirichlet.setConstantBC(0.)
    bc.u_advective.uOfXT = None
    #bc.u_diffusive.setConstantBC(0.) # if enabled, no BL forms

    bc.v_dirichlet.setConstantBC(0.)
    bc.v_advective.uOfXT = None
    #bc.v_diffusive.setConstantBC(0.) # if enabled, not BL forms

    bc.w_dirichlet.setConstantBC(0.)
    bc.w_advective.uOfXT = None
    #bc.w_diffusive.setConstantBC(0.) # if enabled, not BL forms

    bc.vof_dirichlet.uOfXT = None
    bc.vof_advective.setConstantBC(0.)
    ##.vof_diffusive.uOfXT = does not exist

    bc.k_dirichlet.setConstantBC(0.0)
    #bc.k_advective.
    #bc.k_diffusive.setConstantBC(0.)

    bc.dissipation_dirichlet.setConstantBC(0.0375)
    #bc.dissipation_advective.
    #bc.dissipation_diffusive.setConstantBC(0.)   

    bc.p_dirichlet.uOfXT = None
    bc.p_advective.setConstantBC(0.)
    #bc.p_diffusive.uOfXT = None

    bc.pInc_dirichlet.uOfXT = None
    bc.pInc_advective.setConstantBC(0.)  
    bc.pInc_diffusive.setConstantBC(0.)

    bc.pInit_dirichlet.uOfXT = None
    bc.pInit_advective.setConstantBC(0.)
    bc.pInit_diffusive.setConstantBC(0.)

    bc.vos_advective.setConstantBC(0.)
    bc.us_dirichlet.setConstantBC(0.)
    bc.vs_dirichlet.setConstantBC(0.)
    bc.ws_dirichlet.setConstantBC(0.)  

