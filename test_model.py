import numpy as np
from scipy.integrate import ode
from shapely.geometry import Polygon, Point, LineString
from shapely import affinity

import matplotlib.pyplot as plt
from descartes.patch import PolygonPatch
from scipy import spatial

import motionmodel
import pushedobject
import IPython

def object_hand_motion(t, x, pusher):
    # velocity of the hand trajectory and the object. 
    # x is the combined state of object and hand.
    dx = np.zeros(np.size(x))
    # dx[3:] = v_push_world
    obj_pose = x
    # contact_p = x[3:]
    v_push_world, pusher_line = pusher
    # x,y,theta = obj_pose
    contact_point_world, contact_normal_world, in_contact = obj.eval_contact(pusher_line,obj_pose,v_push_world)
    if not in_contact:
        print 'No contact'
        return False
    
    T = np.array([[np.cos(obj_pose[2]),-np.sin(obj_pose[2]),obj_pose[0]],
                  [np.sin(obj_pose[2]), np.cos(obj_pose[2]),obj_pose[1]],
                  [0,0,1]])
    contact_point_world = [contact_point_world[0],contact_point_world[1],1.]
    contact_point = np.array([25.2,-13])#np.dot(np.linalg.inv(T),contact_point_world)[:2]
    contact_normal = np.array([0,-1.])#np.dot(np.linalg.inv(T[:2,:2]),contact_normal_world)
    
    print contact_point, np.dot(np.linalg.inv(T[:2,:2]),contact_point_world[:2]-T[:2,2]),np.dot(np.linalg.inv(T),contact_point_world)[:2], x
    v_push = np.dot(np.linalg.inv(T[:2,:2]),v_push_world)
    v_o, omega, v_slip, contact_mode = motion.compute_vel_single_contact(v_push,contact_point,contact_normal,contact_mu)
    if contact_mode is not 'sticking':
        return False
    R = T[:2,:2]
    Adj = np.zeros((3,3))
    Adj[:2,:2] = R
    Adj[:2,2] = obj_pose[:2]
    Adj[2,2] = 1
    dx_local = [v_o[0], v_o[1], omega]
    dx = np.dot(Adj,dx_local)
    # IPython.embed()
    return dx

# import IPython

poly = Polygon([(63./2.,26./2.),(-63./2.,26./2),(-63./2.,-26./2.),(63./2.,-26./2.)])
pose = [0.,0.,0.] #x,y,theta in world frame
pose = np.array([[np.cos(pose[2]),-np.sin(pose[2]),pose[0]],
              [np.sin(pose[2]), np.cos(pose[2]),pose[1]],
              [0,0,1]])
num_support_points = 100
pressure_option = 'uniform'

obj = pushedobject.PushedObject('polygon',poly,pose,num_support_points,pressure_option)

contact_mu = 0.3

motion = motionmodel.MotionModel(obj, contact_mu)
v_push_world = np.array([0.,20.]) # mm/s (in world)
pusher_starting_point = np.array([63/2.5,-30./2.]) # mm (in world)
pusher_line = LineString([pusher_starting_point,pusher_starting_point+1000*v_push_world/np.linalg.norm(v_push_world)])
pusher = (v_push_world, pusher_line)
contact_point_world, contact_normal_world, in_contact = obj.eval_contact(pusher_line,pose,v_push_world)

# solver = ode(object_hand_motion)
# solver.set_integrator('dopri5')
# solver.set_f_params(pusher)
# t0 = 0
# x0 = pose
# solver.set_initial_value(x0,t0)

# t1 = 0.7
# sols = []
# sols.append(x0)
# dt = 1e-1
# t = t0
# # Repeatedly call the `integrate` method to advance the
# # solution to time t[k], and save the solution in sol[k].
# k = 1
# while solver.successful() and solver.t < t1:
#     t = t+dt
#     solver.integrate(t)
#     sols.append(solver.y)
#     k += 1

# # Plot the solution...
# fig = plt.figure(1)
# ax = fig.add_subplot(111)
# for k in range(len(sols)):
#     new_poly = affinity.rotate(poly,sols[k][2], use_radians=True)
#     new_poly = affinity.translate(new_poly,sols[k][0],sols[k][1])
#     patch = PolygonPatch(new_poly, alpha=0.25, zorder=2)
#     contact_point_world, contact_normal_world, in_contact = obj.eval_contact(pusher_line,sols[k],v_push_world)
#     # T = np.array([[np.cos(sols[k][2]),-np.sin(sols[k][2]),sols[k][0]],
#     #               [np.sin(sols[k][2]), np.cos(sols[k][2]),sols[k][1]],
#     #               [0,0,1]])
#     # contact_point_world = [contact_point_world[0],contact_point_world[1],1.]
#     # print np.dot(np.linalg.inv(T),contact_point_world)
#     ax.plot(contact_point_world[0],contact_point_world[1], 'o', color='#999999')
#     ax.add_patch(patch)
# plt.grid(True)
# ax.set_xlim([-70,70])
# ax.set_ylim([-50,50])
# ax.set_title('Visualization')
# plt.show(False)
# plt.axis('equal')
# raw_input()


contact_point_world = [contact_point_world[0],contact_point_world[1],1]
contact_point = np.dot(np.linalg.inv(pose),contact_point_world)[:2]
contact_normal = np.dot(np.linalg.inv(pose[:2,:2]),contact_normal_world)
v_push = np.dot(np.linalg.inv(pose[:2,:2]),v_push_world)
v_obj, omega, v_slip, contact_mode = motion.compute_vel_single_contact(v_push,contact_point,contact_normal,contact_mu)


fig = plt.figure(1)
ax = fig.add_subplot(111)
patch = PolygonPatch(obj.poly_world, alpha=0.25, zorder=2)
ax.add_patch(patch)
ax.plot(contact_point_world[0],contact_point_world[1], 'o', color='#999999')

discrtimestep= 1e-2
t=0
while t<3:
    t +=discrtimestep
    pose = obj.eval_pose(discrtimestep, v_obj, omega)
    contact_point_world, contact_normal_world, in_contact = obj.eval_contact(pusher_line,pose,v_push_world)
    if in_contact:
        obj.update_pose(pose)
        patch = PolygonPatch(obj.poly_world, alpha=0.25, zorder=2)
        ax.add_patch(patch)
        ax.plot(contact_point_world[0],contact_point_world[1], 'o', color='#999999')
    else:
        print 'No contact'
        break
    contact_point_world = [contact_point_world[0],contact_point_world[1],1]
    contact_point = np.dot(np.linalg.inv(pose),contact_point_world)[:2]
    contact_normal = np.dot(np.linalg.inv(pose[:2,:2]),contact_normal_world)
    v_push = np.dot(np.linalg.inv(pose[:2,:2]),v_push_world)

    # IPython.embed()
    v_obj, omega, v_slip, contact_mode = motion.compute_vel_single_contact(v_push,contact_point,contact_normal,contact_mu)
    # print contact_point, contact_normal, v_slip#, v_push, contact_mode, contact_point_world
    # if contact_mode is not 'sticking':
        # print contact_mode
        # IPython.embed()
        # break
print obj.pose
plt.grid(True)
ax.set_xlim([-70,70])
ax.set_ylim([-50,50])
ax.set_title('Visualization')
plt.show(False)
plt.axis('equal')
raw_input()
