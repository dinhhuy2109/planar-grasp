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

poly = Polygon([(100./2.,45./2.),(-100./2.,45./2),(-100./2.,-45./2.),(100./2.,-45./2.)])
pose_cov = np.array([[16.,3,0],
                     [3,36,0],
                     [0,0,0.04]])
pose_mean = np.array([-5.,4.,0]) # x,y,theta in world frame
pose = np.random.multivariate_normal(pose_mean, pose_cov,1)[0]
print pose 
pose = np.array([[np.cos(pose[2]),-np.sin(pose[2]),pose[0]],
              [np.sin(pose[2]), np.cos(pose[2]),pose[1]],
              [0,0,1.]])
num_support_points = 10
pressure_option = 'uniform'

obj = pushedobject.PushedObject('polygon',poly,pose,num_support_points,pressure_option)

contact_mu = 0.5

motion = motionmodel.MotionModel()

v_push1 = np.array([0.,20.]) # mm/s (in world)
pusher1 = pusher.Pusher(poly=LineString([(-11,-42.5),(11,-42.5)]),pose=np.eye(3),v_push=v_push1)

v_push2 = np.array([0.,-20.]) # mm/s (in world)
pusher2 = pusher.Pusher(poly=LineString([(-11,42.5),(11,42.5)]),pose=np.eye(3),v_push=v_push2)

fig = plt.figure(1)
ax = fig.add_subplot(111)
plot_pusher(ax,pusher1.poly)
plot_pusher(ax,pusher2.poly)
patch = PolygonPatch(obj.poly_world, alpha=0.25, zorder=2)
ax.add_patch(patch)
# raw_input()
print ('--- start simulation ---')
start_time = time.time()
contacts = obj.eval_contact_with_2shortfingers(pusher1,pusher2,pose)
print contacts

if contacts[IN_CONTACT]:
    if len(contacts[POINTS]) == 1: # single contact
        v_obj,contact_mode = motion.compute_vel_single_contact(contacts[V_PUSHES][0], contacts[POINTS][0],contacts[NORMALS][0],contact_mu,obj.ls_coeff)
    else: # multi contacts
        v_obj = motion.compute_vel_multi_contacts(contacts[V_PUSHES],contacts[POINTS],contacts[NORMALS],contact_mu,obj.ls_A)
    print v_obj
    if np.linalg.norm(v_obj) < 1e-3:
        print 'JAMMING or STOP MOVING'
        
t_contact = contacts[TIME]
pusher1.update_pusher(pusher1.eval_pusher(t_contact))
pusher2.update_pusher(pusher2.eval_pusher(t_contact))
plot_pusher(ax,pusher1.poly)
plot_pusher(ax,pusher2.poly)
discrtimestep= 1e-1
t=0
prv_num_contact = len(contacts[POINTS])
min_dicrtimestep = 8e-3
while t<5:
    # TODO check pose
    delta = (np.arccos(obj.pose[0,0]))
    # if discrtimestep>1e-3:
    #     if delta < 0.1:
    #         discrtimestep = discrtimestep/10.
    #         print discrtimestep
    t +=discrtimestep
    pose = obj.eval_pose(discrtimestep, v_obj)
    obj.update_pose(pose)
    
    contacts = obj.eval_contact_with_2shortfingers(pusher1,pusher2,pose)
    num_contact = len(contacts[POINTS])
    if num_contact > prv_num_contact and num_contact!=4:
        pusher1.update_pusher(pusher1.eval_pusher(-discrtimestep))
        pusher2.update_pusher(pusher2.eval_pusher(-discrtimestep))
        obj.eval_pose(-discrtimestep, v_obj)
        obj.update_pose(pose)
        t -=discrtimestep
        if discrtimestep>min_dicrtimestep:
            discrtimestep = discrtimestep/2.
        # print 'here'
        prv_num_contact = num_contact
        continue
    if num_contact == 1: # single contact
        v_obj,contact_mode = motion.compute_vel_single_contact(contacts[V_PUSHES][0], contacts[POINTS][0],contacts[NORMALS][0],contact_mu,obj.ls_coeff)
        # print 'single contact'
    else: # multi contacts
        v_obj = motion.compute_vel_multi_contacts(contacts[V_PUSHES],contacts[POINTS],contacts[NORMALS],contact_mu,obj.ls_A)
        # print 'multiple contacts', len(contacts[POINTS])
    if  np.linalg.norm(v_obj[:2])<1e-4 and np.abs(v_obj[2]) < 1e-4:
        if discrtimestep<=min_dicrtimestep:
            print pusher1.points[0], pusher2.points[0]
            print obj.pose
            # import IPython
            # IPython.embed()
            print 'STOP MOVING at t=',t,'s'
            break
        pusher1.update_pusher(pusher1.eval_pusher(-discrtimestep))
        pusher2.update_pusher(pusher2.eval_pusher(-discrtimestep))
        obj.eval_pose(-discrtimestep, v_obj)
        obj.update_pose(pose)
        t -=discrtimestep
        # print 'there'
        continue
    if np.linalg.norm(pusher1.points[0] - pusher2.points[0]) <= 45:
        # print np.linalg.norm(pusher1.points[0] - pusher2.points[0])
        # pusher1.update_pusher(pusher1.eval_pusher(-discrtimestep))
        # pusher2.update_pusher(pusher2.eval_pusher(-discrtimestep))
        # print np.linalg.norm(pusher1.points[0] - pusher2.points[0])
        # if discrtimestep<=min_dicrtimestep:
        print 'Finger closed'
        break
        continue
    # print num_contact, prv_num_contact
    prv_num_contact = num_contact
    pusher1.update_pusher(pusher1.eval_pusher(discrtimestep))
    pusher2.update_pusher(pusher2.eval_pusher(discrtimestep))
    # plot_pusher(ax,pusher1.poly)
    # plot_pusher(ax,pusher2.poly)
    # patch = PolygonPatch(obj.poly_world, alpha=0.25, zorder=2)
    # ax.add_patch(patch)
    
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
