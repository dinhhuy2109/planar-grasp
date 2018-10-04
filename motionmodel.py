import numpy as np
import lemkelcp as lcp
import IPython

class MotionModel(object):
    """
    This class is the motion model of pushing operation under quasi-static assumption.
    """

    def __init__(self):
        print 'Motion model init DONE'

    def compute_friction_cone_edges(self,contact_points,contact_normals,contact_mu):
        num_contacts = len(contact_points)
        fc_edges = []
        for i in range(num_contacts):
            fn = -contact_normals[i]
            ft = contact_mu*np.array([-fn[1],fn[0]])
            f1 = fn+ft
            f2 = fn-ft
            fl = f1/np.linalg.norm(f1) if np.linalg.norm(f1) != 0 else f1
            fr = f2/np.linalg.norm(f2) if np.linalg.norm(f2) != 0 else f2
            tau_l = contact_points[i][0]*fl[1] - contact_points[i][1]*fl[0]
            tau_r = contact_points[i][0]*fr[1] - contact_points[i][1]*fr[0]
            fc_edges.append([fl,fr,tau_l,tau_r])
        return fc_edges

    def compute_vel_multi_contacts(self,v_pushes, contact_points,contact_normals,contact_mu,A):
        num_cts = len(contact_points)
        N = np.zeros((num_cts,3))
        L = np.zeros((2*num_cts,3))
        a = np.zeros(num_cts)
        b = np.zeros(2*num_cts)
        E = np.zeros((2*num_cts,num_cts))
        Mus = np.diag(contact_mu*np.ones(num_cts))
        F = np.zeros(3)
        V = np.zeros(3)

        for i in range(num_cts):
            E[2*i:2*i+2,i] = [1,1]
            Jp = np.array([[1,0,-contact_points[i][1]],
                           [0,1,contact_points[i][0]]])
            N[i,:] = np.dot(np.transpose(-contact_normals[i]),Jp)
            D = np.array([[-contact_normals[i][1], contact_normals[i][1]],
                          [contact_normals[i][0] ,-contact_normals[i][0]]])
            L[2*i:2*i+2,:] = np.dot(np.transpose(D),Jp)
            a[i] = np.dot(np.transpose(contact_normals[i]),v_pushes[i])
            b[2*i:2*i+2] = np.dot(-np.transpose(D),v_pushes[i])

        q = np.zeros(num_cts*4)
        q[:num_cts] = a
        q[num_cts:3*num_cts] = b
        
        lcp_M = np.zeros((4*num_cts,4*num_cts))
        
        lcp_M[:num_cts,:num_cts] = np.dot(np.dot(N,A),np.transpose(N))
        lcp_M[num_cts:3*num_cts,:num_cts] = np.dot(np.dot(L,A),np.transpose(N))
        lcp_M[3*num_cts:,:num_cts] = Mus

        lcp_M[:num_cts,num_cts:3*num_cts] =  np.dot(np.dot(N,A),np.transpose(L))
        lcp_M[num_cts:3*num_cts,num_cts:3*num_cts] = np.dot(np.dot(L,A),np.transpose(L))
        lcp_M[3*num_cts:,num_cts:3*num_cts] = -np.transpose(E)

        lcp_M[num_cts:3*num_cts,3*num_cts:] = E

        sol = lcp.lemkelcp(lcp_M,q)

        if sol[2] == 'Solution Found':
            z = sol[0]
            w = np.dot(lcp_M,z) + q
            if (np.sum(w < -1e-3)==0 and np.sum(z<-1e-3)==0 and np.abs(np.dot(np.transpose(w),z)) <=1e-3):
                flag = True
            else:
                flag = False
        else:
            flag = False
            
        if flag:
            fns = z[:num_cts]
            fts = z[num_cts:3*num_cts]
            V = np.dot(A,np.dot(np.transpose(N),fns)+np.dot(np.transpose(L),fts))
            F = np.linalg.solve(A,V)
        return V
    
    def compute_vel_single_contact(self,v_push,contact_point,contact_normal,contact_mu,ls_coeff):
        c = ls_coeff
        [fc_edges] = self.compute_friction_cone_edges([contact_point],[contact_normal],contact_mu)
        vc_edges = self.compute_motion_cone(fc_edges,ls_coeff)
        contact_mode = self.determine_contact_model(v_push, vc_edges)
        x_c = contact_point[0]
        y_c = contact_point[1]
        if contact_mode is 'sticking':
            v_o = v_push
            v_x = ((c**2+x_c**2)*v_o[0] + x_c*y_c*v_o[1])/(c**2+x_c**2+y_c**2)
            v_y = ((c**2+y_c**2)*v_o[1] + x_c*y_c*v_o[0])/(c**2+x_c**2+y_c**2)
            v_obj = np.array([v_x,v_y])
            omega = (x_c*v_y-y_c*v_x)/(c**2)
            v_slip = np.zeros(2)
        elif contact_mode is 'leftsliding':
            k = np.dot(v_push,-contact_normal)/np.dot(vc_edges[1],-contact_normal)
            v_o = k*vc_edges[1]
            v_x = ((c**2+x_c**2)*v_o[0] + x_c*y_c*v_o[1])/(c**2+x_c**2+y_c**2)
            v_y = ((c**2+y_c**2)*v_o[1] + x_c*y_c*v_o[0])/(c**2+x_c**2+y_c**2)
            v_obj = np.array([v_x,v_y])
            omega = (x_c*v_o[1]-y_c*v_o[0])/(c**2)
            v_slip = v_push - v_o
        elif contact_mode is 'rightsliding':
            k = np.dot(v_push,-contact_normal)/np.dot(vc_edges[0],-contact_normal)
            v_o = k*vc_edges[0]
            v_x = ((c**2+x_c**2)*v_o[0] + x_c*y_c*v_o[1])/(c**2+x_c**2+y_c**2)
            v_y = ((c**2+y_c**2)*v_o[1] + x_c*y_c*v_o[0])/(c**2+x_c**2+y_c**2)
            v_obj = np.array([v_x,v_y])
            omega = (x_c*v_y-y_c*v_x)/(c**2)
            v_slip = v_push - v_o
        vel = [v_obj[0], v_obj[1], omega]
        return vel, contact_mode

    def compute_motion_cone(self,fc_edges,ls_coeff):
        v_lx = ls_coeff**2*(fc_edges[0][0]/fc_edges[2])*np.sign(fc_edges[2]) # v_lx/omega
        v_ly = ls_coeff**2*(fc_edges[0][1]/fc_edges[2])*np.sign(fc_edges[2])# v_ly/omega
        v_rx = ls_coeff**2*(fc_edges[1][0]/fc_edges[3])*np.sign(fc_edges[3])
        v_ry = ls_coeff**2*(fc_edges[1][1]/fc_edges[3])*np.sign(fc_edges[3])
        vl = np.array([v_lx,v_ly])
        vr = np.array([v_rx,v_ry])
        vc_edges = [vl,vr]
        # import IPython
        # IPython.embed()
        return vc_edges

    def determine_contact_model(self,v_push, vc_edges):
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
