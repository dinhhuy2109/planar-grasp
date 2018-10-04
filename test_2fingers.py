import numpy as np
import copy
from shapely.geometry import Polygon, Point, LineString
from shapely import affinity

import matplotlib.pyplot as plt
from descartes.patch import PolygonPatch
from scipy import spatial

import motionmodel
import pushedobject
import pusher
import IPython

import time

def plot_pusher(ax,pusher):
    x, y = pusher.xy
    ax.plot(x, y, 'o', color='#6699cc', zorder=1)
    ax.plot(x, y, color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
    return True
TIME = 4
IN_CONTACT = 3
V_PUSHES = 2
NORMALS = 1
POINTS = 0

poly = Polygon([(63./2.,26./2.),(-63./2.,26./2),(-63./2.,-26./2.),(63./2.,-26./2.)])
pose = [0.,0.,-np.pi/6.] # x,y,theta in world frame
pose = np.array([[np.cos(pose[2]),-np.sin(pose[2]),pose[0]],
              [np.sin(pose[2]), np.cos(pose[2]),pose[1]],
              [0,0,1]])
num_support_points = 10
pressure_option = 'uniform'

obj = pushedobject.PushedObject('polygon',poly,pose,num_support_points,pressure_option)

contact_mu = 0.5

motion = motionmodel.MotionModel()

v_push1 = np.array([0.,10.]) # mm/s (in world)
pusher1 = pusher.Pusher(poly=LineString([(-50,-40),(50,-40)]),pose=np.eye(3),v_push=v_push1)

v_push2 = np.array([0.,-10.]) # mm/s (in world)
pusher2 = pusher.Pusher(poly=LineString([(-50,40),(50,40)]),pose=np.eye(3),v_push=v_push2)

fig = plt.figure(1)
ax = fig.add_subplot(111)
plot_pusher(ax,pusher1.poly)
plot_pusher(ax,pusher2.poly)
patch = PolygonPatch(obj.poly_world, alpha=0.25, zorder=2)
ax.add_patch(patch)

print ('--- start simulation ---')
start_time = time.time()
contacts_in_world = obj.eval_contact_with_2pushers(pusher1,pusher2,pose)

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
        v_obj,contact_mode = motion.compute_vel_single_contact(contacts[V_PUSHES][0], contacts[POINTS][0],contacts[NORMALS][0],contact_mu,obj.ls_coeff)
    else: # multi contacts
        v_obj = motion.compute_vel_multi_contacts(contacts[V_PUSHES],contacts[POINTS],contacts[NORMALS],contact_mu,obj.ls_A)
        if np.linalg.norm(v_obj) < 1e-3:
            print 'JAMMING or STOP MOVING'
else:
    print 'NO CONTACT'

t_contact = contacts_in_world[TIME]
pusher1.update_pusher(pusher1.eval_pusher(t_contact))
pusher2.update_pusher(pusher2.eval_pusher(t_contact))

discrtimestep= 2e-1
t=0
previous_pose = obj.pose
while t<5:
    t +=discrtimestep
    pose = obj.eval_pose(discrtimestep, v_obj)
    contacts_in_world = obj.eval_contact_with_2pushers(pusher1,pusher2,pose,tolerance = 10.*discrtimestep)
    
    if contacts_in_world[IN_CONTACT]:
        obj.update_pose(pose)
        pusher1.update_pusher(pusher1.eval_pusher(discrtimestep))
        pusher2.update_pusher(pusher2.eval_pusher(discrtimestep))
        plot_pusher(ax,pusher1.poly)
        plot_pusher(ax,pusher2.poly)
        patch = PolygonPatch(obj.poly_world, alpha=0.25, zorder=2)
        for pt in contacts_in_world[POINTS]:
            ax.add_patch(patch)
            ax.plot(pt[0],pt[1], 'o', color='#999999')
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
        v_obj,contact_mode = motion.compute_vel_single_contact(contacts[V_PUSHES][0], contacts[POINTS][0],contacts[NORMALS][0],contact_mu,obj.ls_coeff)
        # print 'single contact'
    else: # multi contacts
        v_obj = motion.compute_vel_multi_contacts(contacts[V_PUSHES],contacts[POINTS],contacts[NORMALS],contact_mu,obj.ls_A)
        # print 'multiple contacts'
    # print contacts_in_world[POINTS]
    # print v_obj
    # raw_input()
    if np.linalg.norm(previous_pose-pose)<1e-5:
        print 'STOP MOVING'
        break
    else:
        previous_pose = pose
                
    if np.linalg.norm(v_obj[:2])+np.abs(v_obj[2]) < 1e-3:
        if discrtimestep > 1e-2: # to have a fine simulation
            discrtimestep= discrtimestep**-2.
            continue
        else:
            print 'STOP MOVING'
            print 'at t=',t,'s'
            break

print ("--- %s seconds ---" % (time.time() - start_time))
print 'Final obj pose: ', obj.pose
patch = PolygonPatch(obj.poly_world, color = 'r', alpha=0.7, zorder=2)
ax.add_patch(patch)

plt.grid(True)
ax.set_xlim([-70,70])
ax.set_ylim([-50,50])
ax.set_title('Visualization')
plt.show(False)
plt.axis('equal')
raw_input()
