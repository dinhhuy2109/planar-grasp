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
            self.points = np.asarray(self.poly.exterior.coords)
            self.edges = [[self.points[i],self.points[i+1]] for i in range(len(self.points)-1)]
            self.points = self.points[:-1]
            self.centroid = np.asarray(poly.centroid.coords).reshape(2)
            self.edge_normals = [self.compute_normal(edge,self.centroid) for edge in self.edges]
            self.poly_world = self.compute_poly(pose)
            self.support_points = []
            self.num_support_points = num_support_points
            self.grid_support_points(num_support_points,'polygon')
            self.pressure_option = pressure_option
            self.presure_weights = self.assign_pressure(self.support_points)
            self.ls_coeff = self.compute_limit_surface_coeffs()
            self.ls_A = self.compute_limit_surface_A()
            # IPython.embed()
            
        print 'Object init DONE'

    def compute_normal(self,edge,centroid):
        [point1,point2] = edge
        vec1 = point1-point2
        vec2 = point1-centroid
        normal = np.array([vec1[1],-vec1[0]])
        if np.dot(normal,vec2)>0:
            normal = -normal
        return normal/np.linalg.norm(normal)

    def rotate_vector(self, vector, angle):
        R = np.array([[np.cos(angle),-np.sin(angle)],
                      [np.sin(angle), np.cos(angle)]])
        return np.dot(R,vector)
    
    def eval_pose(self, discrtimestep, v_obj):
        x,y,theta= discrtimestep*np.asarray(v_obj)
        T_ = np.array([[np.cos(theta),-np.sin(theta),x],
                      [np.sin(theta), np.cos(theta),y],
                      [0.,0.,1.]])
        return np.dot(self.pose,T_)

    def update_pose(self, pose):
        self.pose = pose
        pts = list(self.poly.exterior.coords)[:-1]
        new_pts = []
        for pt in pts:
            new_pts.append(np.dot(self.pose[:2,:2],pt)+self.pose[:2,2])
        self.poly_world = Polygon(new_pts)
        return True

    def compute_poly(self,pose):
        pts = list(self.poly.exterior.coords)[:-1]
        new_pts = []
        for pt in pts:
            new_pts.append(np.dot(pose[:2,:2],pt)+pose[:2,2])
        return Polygon(new_pts)

    def compute_poly_given_pose(self,poly,pose):
        pts = list(poly.exterior.coords)[:-1]
        new_pts = []
        for pt in pts:
            new_pts.append(np.dot(pose[:2,:2],pt)+pose[:2,2])
        return Polygon(new_pts)

    def eval_contact_with_2shortfingers(self,pusher1,pusher2,obj_pose):
        intersect_line1 = LineString([pusher1.points[0],pusher2.points[0]])
        intersect_line2 = LineString([pusher1.points[1],pusher2.points[1]])
        contact_points = []
        contact_normals = []
        v_pushes = []
        new_poly = self.compute_poly(obj_pose)
        intersect1 = intersect_line1.intersection(new_poly)
        intersect2 = intersect_line2.intersection(new_poly)
        # import IPython
        # IPython.embed()
        if not (intersect1.is_empty or intersect1.is_empty):
            intersects = [np.array(intersect1.coords)]
            intersects.append(np.array(intersect2.coords))
            intersects = np.reshape(intersects,(4,2))
            dists = [np.linalg.norm(np.array(pusher1.points[0] - intersects[0]))]
            dists.append(np.linalg.norm(np.array(pusher2.points[0] - intersects[1])))
            dists.append(np.linalg.norm(np.array(pusher1.points[1] - intersects[2])))
            dists.append(np.linalg.norm(np.array(pusher2.points[1] - intersects[3])))
            minimum = min(dists)
            for i in range(len(dists)):
                if dists[i] == minimum:
                    contact_points.append(np.dot(np.linalg.inv(obj_pose),np.hstack((intersects[i],1.)))[:2]) # in object frame
                    if i%2 == 1:
                        contact_normals.append(-self.edge_normals[0]) # in object frame
                        t_travel=dists[i]/np.linalg.norm(pusher2.v_push)
                        v_pushes.append(np.dot(np.linalg.inv(obj_pose[:2,:2]),pusher2.v_push)) # in object frame
                    else:
                        contact_normals.append(-self.edge_normals[2])
                        t_travel=dists[i]/np.linalg.norm(pusher1.v_push)
                        v_pushes.append(np.dot(np.linalg.inv(obj_pose[:2,:2]),pusher1.v_push))
            IN_CONTACT = True
            return contact_points, contact_normals, v_pushes,  IN_CONTACT, t_travel
        else: # assuming uncertainty is small such that 2 fingers is always near the middle of the object
            raise ValueError('A very specific bad thing happened. Uncertainty is larger than our assumption')
            
                
    
    def eval_contact_with_2pushers(self,pusher1,pusher2,obj_pose,tolerance = 1):
        pusher1_contacts = self.eval_contact_pusher_obj(pusher1,obj_pose,tolerance)
        pusher2_contacts = self.eval_contact_pusher_obj(pusher2,obj_pose,tolerance)
        
        if pusher1_contacts[3] and pusher2_contacts[3]: #both may be in contact
            if pusher1_contacts[-1] - pusher2_contacts[-1] > 1e-4: # pusher 2 in contact first
                return pusher2_contacts
            elif pusher2_contacts[-1] - pusher1_contacts[-1] > 1e-4: # pusher 1 in contact first
                return pusher1_contacts 
            else:
                 return pusher1_contacts[0]+pusher2_contacts[0], pusher1_contacts[1]+pusher2_contacts[1], pusher1_contacts[2]+pusher2_contacts[2], True, pusher1_contacts[-1], '2 pushers'
        elif pusher1_contacts[3]:
            return pusher1_contacts
        elif pusher2_contacts[3]:
            return pusher2_contacts
        return pusher1_contacts # no contact
       
    def eval_contact_pusher_obj(self,pusher,obj_pose,tolerance = 1):
        # inv_pose = np.linalg.inv(obj_pose)
        # inv_pusher = self.compute_poly_given_pose(pusher,inv_pose)
        new_poly = self.compute_poly(obj_pose)
        points = list(new_poly.exterior.coords)[:-1]
        dummy_d = 1000.
        contact_points = []
        contact_normals = []
        v_pushes = []
        dists = []
        t_travel = 0 #from the init pusher to current pose
        for i in range(len(points)):
            pt = points[i]
            projection_line = LineString([pt,pt-pusher.v_push*dummy_d])
            intersect = projection_line.intersection(pusher.poly)
            if not intersect.is_empty:
                dists.append(np.linalg.norm(np.array(pt)-np.array(intersect.coords)[0]))
            else:
                dists.append(dummy_d)
        minimum = min(dists)
        if minimum != dummy_d:
            indices = []
            for i in range(len(dists)):
                if (dists[i]-minimum) < tolerance:
                    indices.append(i)
                    contact_points.append(points[i])
                    contact_normals.append(-pusher.normal)
                    t_travel=dists[i]/np.linalg.norm(pusher.v_push)
                    v_pushes.append(pusher.v_push)
            IN_CONTACT = True
        else:
            IN_CONTACT = False
        return contact_points, contact_normals, v_pushes,  IN_CONTACT, t_travel
    
    def eval_contact(self,pusher_line,pose,v_push):
        new_poly = self.compute_poly(pose)
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

    def compute_limit_surface_A(self):
        A =  np.diag([1.,1.,self.ls_coeff**2])
        return A
