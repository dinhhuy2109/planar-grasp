import numpy as np
from shapely.geometry import Polygon, Point, LineString
from shapely import affinity

import matplotlib.pyplot as plt
from descartes.patch import PolygonPatch
from scipy import spatial

import motionmodel
import pushedobject


poly = Polygon([(63./2.,26./2.),(-63./2.,26./2),(-63./2.,-26./2.),(63./2.,-26./2.)])
pose = [10.,0.,0.] #x,y,theta in world frame
num_support_points = 100
pressure_option = 'uniform'
obj = pushedobject.PushedObject('polygon',poly,pose,num_support_points,pressure_option)

contact_mu = 0.5

motion = motionmodel.MotionModel(obj, contact_mu)
v_push = np.array([0.,20.]) # mm/s
pusher_starting_point = np.array([63/2.5+10,-30./2.])

pusher_line = LineString([pusher_starting_point,pusher_starting_point+1000*v_push/np.linalg.norm(v_push)])

contact_point_world, contact_normal_world = obj.eval_contact(pusher_line,pose,v_push)
x,y,theta = obj.pose
T = np.array([[np.cos(theta),-np.sin(theta),x],
              [np.sin(theta), np.cos(theta),y],
              [0,0,1]])
contact_point_world = [contact_point_world[0],contact_point_world[1],1]
contact_point = np.dot(np.linalg.inv(T),contact_point_world)
contact_normal = np.dot(np.linalg.inv(T[:2,:2]),contact_normal_world)
v_o, omega, v_slip, contact_mode = motion.compute_vel_single_contact(v_push,contact_point,contact_normal,contact_mu)


fig = plt.figure(1)
ax = fig.add_subplot(111)

# # new_poly = affinity.rotate(obj.poly,obj.pose[2],origin='centroid',use_radians=True)
# # new_poly = affinity.translate(obj.poly,obj.pose[0],obj.pose[1])
# # patch = PolygonPatch(new_poly, alpha=0.25, zorder=2)
# # ax.add_patch(patch)
# # ax.plot(contact_point[0],contact_point[1], 'o', color='#999999')
# # print v_o, omega, v_slip, contact_mode


discrtimestep= 1e-2
t=0
while t<2:
    t +=discrtimestep
    obj.update_pose(discrtimestep, v_o, omega)
    contact_point_world, contact_normal_world = obj.eval_contact(pusher_line,obj.pose,v_push)
    x,y,theta = obj.pose
    T = np.array([[np.cos(theta),-np.sin(theta),x],
                  [np.sin(theta), np.cos(theta),y],
                  [0,0,1]])
    contact_point_world = [contact_point_world[0],contact_point_world[1],1]
    contact_point = np.dot(np.linalg.inv(T),contact_point_world)
    contact_normal = np.dot(np.linalg.inv(T[:2,:2]),contact_normal_world)

    new_poly = affinity.affine_transform(obj.poly,[np.cos(theta),-np.sin(theta),np.sin(theta),np.cos(theta),obj.pose[0],obj.pose[1]])
    patch = PolygonPatch(new_poly, alpha=0.25, zorder=2)
    ax.add_patch(patch)
    ax.plot(contact_point_world[0],contact_point_world[1], 'o', color='#999999')
    
    v_o, omega, v_slip, contact_mode = motion.compute_vel_single_contact(v_push,contact_point,contact_normal,contact_mu)
    print contact_point, contact_normal, contact_mode
    if contact_mode is not 'sticking':
        break
    ax.set_xlim([-70,70])
    ax.set_ylim([-20,50])
    ax.set_title('Visualization')
    plt.show(False)
    raw_input()




raw_input()
