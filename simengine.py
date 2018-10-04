import numpy as np
import copy
from shapely.geometry import Polygon, Point, LineString
from shapely import affinity

import matplotlib.pyplot as plt
from descartes.patch import PolygonPatch
from scipy import spatial

import IPython

import time

TIME = 4
IN_CONTACT = 3
V_PUSHES = 2
NORMALS = 1
POINTS = 0

class SimulationEngine(object):
    """
    This class is the motion model of pushing operation under quasi-static assumption.
    """

    def __init__(self, obj, pusher1, pusher2, motionmodel, contact_mu):
        self.obj = obj
        self.pusher1 = pusher1
        self.pusher2 = pusher2
        self.motion = motionmodel
        self.contact_mu = contact_mu
        self.fig = plt.figure(1)
        self.ax = self.fig.add_subplot(111)
        self.plot_pusher(self.ax,self.pusher1.poly)
        self.plot_pusher(self.ax,self.pusher2.poly)
        print 'Simulation Engine init DONE'

    def plot_pusher(self,ax,pusher):
        x, y = pusher.xy
        ax.plot(x, y, 'o', color='#6699cc', zorder=1)
        ax.plot(x, y, color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
        return True

    def roll_out(self, init_obj_pose, show_final_only = True, show_init = False, discrtimestep= 2e-1, t_max = 5):
        self.obj.update_pose(init_obj_pose)
        pose = init_obj_pose
        pusher1 = copy.copy(self.pusher1)
        pusher2 = copy.copy(self.pusher2)
        
        if show_init:
            patch = PolygonPatch(self.obj.poly_world, alpha=0.25, zorder=2)
            self.ax.add_patch(patch)

        # print ('--- start simulation ---')
        start_time = time.time()
        contacts_in_world = self.obj.eval_contact_with_2pushers(pusher1,pusher2,pose)

        if contacts_in_world[IN_CONTACT]:
            contact_points = []
            contact_normals = []
            v_pushes = []
            for i in range(len(contacts_in_world[POINTS])):
                contact_points.append(np.dot(np.linalg.inv(pose),np.hstack((contacts_in_world[POINTS][i],1.)))[:2])
                contact_normals.append(np.dot(np.linalg.inv(pose[:2,:2]),contacts_in_world[NORMALS][i]))
                v_pushes.append(np.dot(np.linalg.inv(pose[:2,:2]),contacts_in_world[V_PUSHES][i]))
            contacts = contact_points, contact_normals, v_pushes, contacts_in_world[IN_CONTACT]
            if len(contacts[POINTS]) == 1: # single contact
                v_obj,contact_mode = self.motion.compute_vel_single_contact(contacts[V_PUSHES][0], contacts[POINTS][0],contacts[NORMALS][0],self.contact_mu,self.obj.ls_coeff)
            else: # multi contacts
                v_obj = self.motion.compute_vel_multi_contacts(contacts[V_PUSHES],contacts[POINTS],contacts[NORMALS],self.contact_mu,self.obj.ls_A)
                if np.linalg.norm(v_obj) < 1e-3:
                    print 'JAMMING or STOP MOVING'
        else:
            print 'NO CONTACT'

        t_contact = contacts_in_world[TIME]
        pusher1.update_pusher(pusher1.eval_pusher(t_contact))
        pusher2.update_pusher(pusher2.eval_pusher(t_contact))
        # IPython.embed()
        t=0
        previous_pose = self.obj.pose
        while t<t_max:
            pose = self.obj.eval_pose(discrtimestep, v_obj)
    
            new_poly = self.obj.compute_poly(pose)
            pusher1.update_pusher(pusher1.eval_pusher(discrtimestep))
            pusher2.update_pusher(pusher2.eval_pusher(discrtimestep))
            if pusher1.poly.intersects(new_poly) or pusher2.poly.intersects(new_poly):
                if discrtimestep > 1e-5:
                    print discrtimestep
                    pusher1.update_pusher(pusher1.eval_pusher(-discrtimestep))
                    pusher2.update_pusher(pusher2.eval_pusher(-discrtimestep))
                    discrtimestep= discrtimestep/10.
                    continue
                else:
                    break
            t +=discrtimestep
            contacts_in_world = self.obj.eval_contact_with_2pushers(pusher1,pusher2,pose,tolerance = 10.*discrtimestep)
            if contacts_in_world[IN_CONTACT]:
                self.obj.update_pose(pose)
                if not show_final_only:
                    # self.plot_pusher(self.ax,pusher1.poly)
                    # self.plot_pusher(self.ax,pusher2.poly)
                    patch = PolygonPatch(self.obj.poly_world, alpha=0.25, zorder=2)
                    for pt in contacts_in_world[POINTS]:
                        self.ax.add_patch(patch)
                        self.ax.plot(pt[0],pt[1], 'o', color='#999999')
            else:
                if len(contacts_in_world[POINTS])==0:
                    print 'No contact'
                    break
            contact_points = []
            contact_normals = []
            v_pushes = []

            for i in range(len(contacts_in_world[POINTS])):
                contact_points.append(np.dot(np.linalg.inv(pose),np.hstack((contacts_in_world[POINTS][i],1.)))[:2])
                contact_normals.append(np.dot(np.linalg.inv(pose[:2,:2]),contacts_in_world[NORMALS][i]))
                v_pushes.append(np.dot(np.linalg.inv(pose[:2,:2]),contacts_in_world[V_PUSHES][i]))
            contacts = contact_points, contact_normals, v_pushes, contacts_in_world[IN_CONTACT]

            if len(contacts[POINTS]) == 1: # single contact
                v_obj,contact_mode = self.motion.compute_vel_single_contact(contacts[V_PUSHES][0], contacts[POINTS][0],contacts[NORMALS][0],self.contact_mu,self.obj.ls_coeff)
                # print 'single contact'
            else: # multi contacts
                v_obj = self.motion.compute_vel_multi_contacts(contacts[V_PUSHES],contacts[POINTS],contacts[NORMALS],self.contact_mu,self.obj.ls_A)
                # print 'multiple contacts'
            # print contacts_in_world[POINTS]
            # print v_obj
            # raw_input()
            
            if np.linalg.norm(previous_pose-pose)<1e-4:
                print 'STOP MOVING'
                break
            else:
                previous_pose = pose
            if np.linalg.norm(v_obj[:2])+np.abs(v_obj[2]) < 1e-5:
                print 'STOP MOVING'
                print 'at t=',t,'s'
                break
            
        print ("--- %s seconds ---" % (time.time() - start_time))
        # print 'Final obj pose: ', self.obj.pose
        patch = PolygonPatch(self.obj.poly_world, color = 'b', alpha=0.2, zorder=2)
        self.ax.add_patch(patch)
        # self.plot_pusher(self.ax,pusher1.poly)
        # self.plot_pusher(self.ax,pusher2.poly)
        # IPython.embed()
        return self.obj.pose

    def show(self,SHOW):
        plt.grid(True)
        self.ax.set_xlim([-70,70])
        self.ax.set_ylim([-50,50])
        self.ax.set_title('Visualization')
        plt.axis('equal')
        plt.show(True)
        raw_input()
        return True
