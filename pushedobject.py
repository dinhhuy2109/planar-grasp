import numpy as np
from shapely.geometry import Polygon, Point, LineString
from shapely import affinity
import IPython

class PushedObject(object):
    """
    Pushed object
    """
    def __init__(self,object_type,poly,pose,num_support_points,pressure_option):
        self.object_type = object_type
        if self.object_type is 'polygon':
            self.pose = pose
            self.poly = poly
            self.poly_world = affinity.rotate(self.poly,self.pose[2], origin='centroid', use_radians=True)
            self.poly_world = affinity.translate(self.poly_world,self.pose[0],self.pose[1])
            self.support_points = []
            self.num_support_points = num_support_points
            self.grid_support_points(num_support_points,'polygon')
            self.pressure_option = pressure_option
            self.presure_weights = self.assign_pressure(self.support_points)
            self.ls_coeff = self.compute_limit_surface_coeffs()

        print 'Object init DONE'


    def rotation(self,pt,angle):
        return np.array([np.cos(angle)*pt[0] - np.sin(angle)*pt[1],np.sin(angle)*pt[0] + np.cos(angle)*pt[1]])

    def compute_pose(self,timestep,v_o,omega):
        d = timestep*np.asarray([v_o[0],v_o[1],omega])
        pts = list(self.poly_world.exterior.coords)[:-1]
        init_pts = list(self.poly.exterior.coords)[:-1]
        new_pts = []
        self.pose = self.pose + d
        for i in range(len(pts)):
            new_pts.append(d[:2]+np.array([d[2]*init_pts[i][-1],-d[2]*init_pts[i][0]])+(pts[i]-np.asarray(self.poly_world.centroid.xy).reshape(2)) + np.asarray(self.poly_world.centroid.xy).reshape(2))


        new_poly = affinity.rotate(self.poly,self.pose[2], origin='centroid', use_radians=True)
        new_poly = affinity.translate(new_poly,self.pose[0],self.pose[1])
        print np.asarray(new_poly.centroid.xy).reshape(2),list(new_poly.exterior.coords)[0] ,'\n', np.asarray(Polygon(new_pts).centroid.xy).reshape(2),list(Polygon(new_pts).exterior.coords)[0]
        
        return Polygon(new_pts)
    
    def eval_pose(self, discrtimestep, v_o, omega):
        displacements = discrtimestep*np.asarray([v_o[0],v_o[1],omega])
        return self.pose + displacements
        # self.pose = self.pose + displacements

    def update_pose(self, pose):
        self.pose = pose
        self.poly_world = affinity.rotate(self.poly,self.pose[2], origin='centroid', use_radians=True)
        self.poly_world = affinity.translate(self.poly_world,self.pose[0],self.pose[1])
        return True
    
    def eval_contact(self,pusher_line,pose,v_push):
        x,y,theta = pose
        new_poly = affinity.rotate(self.poly,pose[2], origin='centroid', use_radians=True)
        new_poly = affinity.translate(new_poly,pose[0],pose[1])
        dummy = pusher_line.intersection(new_poly)
        if not dummy.is_empty:
            intersects = list(dummy.coords)
        else:
            return False, False, False
        new_contact_point = np.array(intersects[0])
        dist = np.linalg.norm(new_contact_point-np.array(pusher_line.coords)[0])
        for pt in intersects[1:]:
            new_dist = np.linalg.norm(pt-np.array(pusher_line.coords)[0])
            if dist > new_dist:
                new_contact_point = pt
                dist = new_dist
        push_vec = v_push/np.linalg.norm(v_push)

        temp_line = LineString([np.array(pusher_line.coords)[0],np.array(pusher_line.coords)[1]+0.1*np.array(push_vec[1],-push_vec[0])])
        temp_intersects = list(temp_line.intersection(new_poly).coords)
        contact_neighbor = np.array(temp_intersects[0])
        dist = np.linalg.norm(contact_neighbor-np.array(temp_line.coords)[0])
        for pt in temp_intersects[1:]:
            new_dist = np.linalg.norm(pt-np.array(temp_line.coords)[0])
            if dist > new_dist:
                contact_neighbor = pt
                dist = new_dist
        tangential = contact_neighbor-new_contact_point
        new_contact_normal = np.array([tangential[1],-tangential[0]])
        if np.dot(new_contact_normal,push_vec)>0:
            new_contact_normal = -new_contact_normal/np.linalg.norm(new_contact_normal)
        else:
            new_contact_normal = new_contact_normal/np.linalg.norm(new_contact_normal)
        return new_contact_point, new_contact_normal, True
    
    def assign_pressure(self,support_points):
        if self.pressure_option is 'uniform':
            pressure_weights = np.ones(len(support_points))
            pressure_weights = pressure_weights/sum(pressure_weights)
            return pressure_weights        
    
    def grid_support_points(self,num_points,option):
        if option is 'polygon':
            object_area = self.poly.area
            
            min_x,min_y,max_x,max_y = self.poly.bounds # enclosing rectangle
            rec_area = (max_x - min_x)*(max_y-min_y)

            grid_size = np.sqrt(object_area/num_points)
            rec_edge_x = np.arange(min_x,max_x,grid_size)
            rec_edge_y = np.arange(min_y,max_y,grid_size)
            for x in rec_edge_x:
                for y in rec_edge_y:
                    grid_center = Point(x,y)
                    if grid_center.within(self.poly):
                        self.support_points.append(np.asarray(grid_center.xy).reshape(2))
            grid = [option,[self.support_points,grid_size]]
            return grid
        # if option == 'rim':
        #     grid = [option]
        #     return 

    def compute_limit_surface_coeffs(self):
        c = 0 # tau_max/f_max
        num_support_points = len(self.support_points)
        for pt in self.support_points:
            c += np.linalg.norm(pt)
        c = c/num_support_points
        return c
