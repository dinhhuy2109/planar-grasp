import numpy as np
from shapely.geometry import Polygon, Point, LineString
from shapely import affinity

import matplotlib.pyplot as plt
from descartes.patch import PolygonPatch
from scipy import spatial

import motionmodel
import pushedobject

import IPython

poly = Polygon([(63./2.,26./2.),(-63./2.,26./2),(-63./2.,-26./2.),(63./2.,-26./2.)])
pose = [0.,0.,0.] #x,y,theta in world frame

num_support_points = 100
pressure_option = 'uniform'


obj = pushedobject.PushedObject('polygon',poly,pose,num_support_points,pressure_option)

contact_mu = 0.5

motion = motionmodel.MotionModel(obj, contact_mu)
v_push_world = np.array([0.,20.]) # mm/s (in world)
pusher_starting_point = np.array([63/2.5,-30./2.]) # mm (in world)
pusher_line = LineString([pusher_starting_point,pusher_starting_point+1000*v_push_world/np.linalg.norm(v_push_world)])


contact_point_world, contact_normal_world, in_contact = obj.eval_contact(pusher_line,pose,v_push_world)
x,y,theta = obj.pose
T = np.array([[np.cos(theta),-np.sin(theta),x],
              [np.sin(theta), np.cos(theta),y],
              [0,0,1]])
contact_point_world = [contact_point_world[0],contact_point_world[1],1]
contact_point = np.dot(np.linalg.inv(T),contact_point_world)[:2]
contact_normal = np.dot(np.linalg.inv(T[:2,:2]),contact_normal_world)
v_push = np.dot(np.linalg.inv(T[:2,:2]),v_push_world)
v_o, omega, v_slip, contact_mode = motion.compute_vel_single_contact(v_push,contact_point,contact_normal,contact_mu)


fig = plt.figure(1)
ax = fig.add_subplot(111)
patch = PolygonPatch(obj.poly_world, alpha=0.25, zorder=2)
ax.add_patch(patch)
ax.plot(contact_point_world[0],contact_point_world[1], 'o', color='#999999')

discrtimestep= 1e-3
t=0
while t<2:
    t +=discrtimestep
    
    # IPython.embed()
    print np.linalg.norm(np.asarray(obj.poly_world.exterior.coords[3]).reshape(2) - np.asarray(contact_point_world[:2]).reshape(2))
    
    pose = obj.eval_pose(discrtimestep, v_o, omega)
    contact_point_world, contact_normal_world, in_contact = obj.eval_contact(pusher_line,pose,v_push_world)
    if in_contact:
        patch = PolygonPatch(obj.poly_world, alpha=0.25, zorder=2)
        ax.add_patch(patch)
        ax.plot(contact_point_world[0],contact_point_world[1], 'o', color='#999999')
        obj.update_pose(pose)
    else:
        print 'No contact'
        break
        
    x,y,theta = obj.pose
    T = np.array([[np.cos(theta),-np.sin(theta),x],
                  [np.sin(theta), np.cos(theta),y],
                  [0.,0.,1.]])
    contact_point_world = [contact_point_world[0],contact_point_world[1],1]
    contact_point = np.dot(np.linalg.inv(T),contact_point_world)[:2]
    contact_normal = np.dot(np.linalg.inv(T[:2,:2]),contact_normal_world)
    v_push = np.dot(np.linalg.inv(T[:2,:2]),v_push_world)

    
    # IPython.embed()
    v_o, omega, v_slip, contact_mode = motion.compute_vel_single_contact(v_push,contact_point,contact_normal,contact_mu)
    print contact_point, contact_normal, v_push, contact_mode, contact_point_world
    if contact_mode is not 'sticking':
        # IPython.embed()
        break

ax.set_xlim([-70,70])
ax.set_ylim([-50,50])
ax.set_title('Visualization')
plt.show(False)
plt.axis('equal')
raw_input()
