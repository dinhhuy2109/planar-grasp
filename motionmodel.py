import numpy as np


class MotionModel(object):
    """
    This class is the motion model of pushing operation under quasi-static assumption.
    """

    def __init__(self, setup=None):


    def compute_friction_cone_edges(self,contact_points,contact_normals,contact_mu):
    """
    return [fl,fr,tau_l,tau_r]
    """
        num_contacts = len(contact_points)
        fc_edges = []
        for i in range(num_contacts):
            fn = -contact_normals[i]
            ft = contact_mu*np.array(-fn[1],fn[0])
            f1 = fn+ft
            f2 = fn-ft
            fl = f1/np.linalg.norm(f1) if np.linalg.norm(f1) != 0 else f1
            fr = f2/np.linalg.norm(f2) if np.linalg.norm(f2) != 0 else f2
            tau_l = np.cross(contact_points[i],fl)
            tau_r = np.cross(contact_points[i],fr)
            fc_edges.append([fl,fr,tau_l,tau_r])
        return fc_edges

    def compute_motion_cones(self,fc_edges,ls_coeff):
        vc_edges = []
        for edges in fc_edges:
            v_lx = ls_coeff**2*(fc_edges[0][0]/fc_edges[2]) # v_lx/omega
            v_ly = ls_coeff**2*(fc_edges[0][1]/fc_edges[2]) # v_ly/omega
            v_rx = ls_coeff**2*(fc_edges[1][0]/fc_edges[3])
            v_ry = ls_coeff**2*(fc_edges[1][1]/fc_edges[3])
            vl = np.array([v_lx,v_ly])
            vr = np.array([v_rx,v_ry])
            vc_edges.append([vl,vr])
        return vc_edges

    def determine_contact_model(self,v_push, point, contact_normal, contact_mu, ls_coeff):
        fc_edges = compute_friction_cone_edges([point],[contact_normal],contact_mu)
        [vc_edges] = compute_motion_cones(fc_edges,ls_coeff)
        # Determine which side of vc_edges vp belongs to
        k1 = v_push[0]*vc_edges[1][1] - v_push[1]*vc_edges[1][0] # right vc + vpush
        k2 = v_push[0]*vc_edges[0][1] - v_push[1]*vc_edges[0][0] # left vc + vpush
        k3 = vc_edges[1][1]*vc_edges[0][0] - vc_edges[1][0]*vc_edges[0][1] # right vc + left vc
        
        if (k1*k2 <= 0):
            contact_mode = 'sticking'
        elif (k1*k3 >= 0): # motion is to the left of the left motion cone edge.
            contact_mode = 'leftsliding'
        else: # motion is to the right of the right motion cone edge.
            contact_mode = 'rightsliding'
        # print contact_mode
        return contact_mode

    def compute_vel_single_contact(self,v_push,point,contact_normal,contact_mu,ls_coeff,contact_mode):
        if contact_mode = 'sticking':
            
            elif contact_mode = 'leftsliding':
                elif contact_mode = 'rightsliding':
                    
