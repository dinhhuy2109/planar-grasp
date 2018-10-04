import numpy as np
from scipy.integrate import ode
from shapely.geometry import Polygon, Point, LineString
from shapely import affinity

import matplotlib.pyplot as plt
from descartes.patch import PolygonPatch
from scipy import spatial

import motionmodel
import pushedobject
import pusher
import IPython

def plot_pusher(ax,pusher):
    x, y = pusher.xy
    ax.plot(x, y, 'o', color='#6699cc', zorder=1)
    ax.plot(x, y, color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
    return True

poly = Polygon([(63./2.,26./2.),(-63./2.,26./2),(-63./2.,-26./2.),(63./2.,-26./2.)])
pose = [0.,0.,-np.pi/6.] #x,y,theta in world frame
pose = np.array([[np.cos(pose[2]),-np.sin(pose[2]),pose[0]],
              [np.sin(pose[2]), np.cos(pose[2]),pose[1]],
              [0,0,1]])
num_support_points = 100
pressure_option = 'uniform'

obj = pushedobject.PushedObject('polygon',poly,pose,num_support_points,pressure_option)
contact_mu = 0.5
motion = motionmodel.MotionModel(obj, contact_mu)
v_push_world = np.array([0.,10.]) # mm/s (in world)

pusher = pusher.Pusher(poly=LineString([(-50,-40),(50,-40)]),pose=np.eye(3),v_push=v_push_world)

contact_point_world, contact_normal_world, IN_CONTACT = obj.eval_contact_pusher_obj(pusher,pose)
if IN_CONTACT and len(contact_point_world)<2:
    contact_point_world = contact_point_world[0]
    contact_normal_world = contact_normal_world[0]
    contact_point_world = [contact_point_world[0],contact_point_world[1],1]
    contact_point = np.dot(np.linalg.inv(pose),contact_point_world)[:2]
    contact_normal = np.dot(np.linalg.inv(pose[:2,:2]),contact_normal_world)
    v_push = np.dot(np.linalg.inv(pose[:2,:2]),v_push_world)

v_obj, omega, v_slip, contact_mode = motion.compute_vel_single_contact(v_push,contact_point,contact_normal,contact_mu)

fig = plt.figure(1)
ax = fig.add_subplot(111)
plot_pusher(ax,pusher.poly)
patch = PolygonPatch(obj.poly_world, alpha=0.25, zorder=2)
ax.add_patch(patch)

discrtimestep= 1e-2
t=0
while t<5:
    t +=discrtimestep
    pose = obj.eval_pose(discrtimestep, v_obj, omega)
    contact_point_world, contact_normal_world, IN_CONTACT= obj.eval_contact_pusher_obj(pusher,pose,tolerance = 10.*discrtimestep)
    # print contact_point_world, contact_normal
    if IN_CONTACT and len(contact_point_world)<2:
        contact_point_world = contact_point_world[0]
        contact_normal_world = contact_normal_world[0]
        obj.update_pose(pose)
        patch = PolygonPatch(obj.poly_world, alpha=0.25, zorder=2)
        ax.add_patch(patch)
        ax.plot(contact_point_world[0],contact_point_world[1], 'o', color='#999999')
    else:
        if len(contact_point_world)<2:
            print 'No contact'
        else:
            print 'Multiple contact'
        break
    contact_point_world = [contact_point_world[0],contact_point_world[1],1]
    contact_point = np.dot(np.linalg.inv(pose),contact_point_world)[:2]
    contact_normal = np.dot(np.linalg.inv(pose[:2,:2]),contact_normal_world)
    v_push = np.dot(np.linalg.inv(pose[:2,:2]),v_push_world)

    # IPython.embed()
    v_obj, omega, v_slip, contact_mode = motion.compute_vel_single_contact(v_push,contact_point,contact_normal,contact_mu)
    print contact_mode, contact_point
    #, contact_normal_world,pusher.normal,v_push_world# v_obj, v_push
    # raw_input()
    
plt.grid(True)
ax.set_xlim([-70,70])
ax.set_ylim([-50,50])
ax.set_title('Visualization')
plt.show(False)
plt.axis('equal')
raw_input()
