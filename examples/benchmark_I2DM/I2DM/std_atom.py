from std_function import *
class Atom:
    def __init__(self,name:str,
                 charge:int,
                 mass:int,
                 radius:float,
                 position:list,
                 energy:float,
                 velocity:list,
                 alpha:float,
                 beta:float)-> None:
        self.atom_name = name             ##type:str
        self.atom_charge = charge         ##type:int,unit:e
        self.atom_mass = mass             ##type:int,unit:Mu
        self.atom_radius = radius         ##unit:A
        self.atom_position = position     ##unit:A,type:mat
        self.atom_energy = energy
        self.atom_velocity_value = get_velocity_value_from_energy(energy,mass)
        value = dis_bet_points(velocity,[0,0,0])
        if value ==0 or value ==1:
            self.atom_velocity = velocity
        else:
            self.atom_velocity = [x/value for x in velocity]
        self.atom_alpha = alpha
        self.atom_beta = beta
    def get_new_energy(self,energy:float) -> None:
        self.atom_energy = energy
        self.atom_velocity_value = get_velocity_value_from_energy(energy,self.atom_mass)
    def get_new_position(self,position:list) -> None:
        self.atom_position = position
    def get_new_velocity(self,velocity:list) -> None:
        self.atom_velocity = get_direction_of_vector(velocity)
    
    def get_neighbor_atom_orders(self,orders:list) -> None:         
        self.neighbor_atom_orders = orders
    def del_neighbor_atom_orders(self) -> None:
        if self.neighbor_atom_orders:
            del self.neighbor_atom_orders

class TargetAtom(Atom):
    def get_threshold_energy(self,energy:float) -> None:
        self.Td = energy   
    def get_atom_layer(self,layer:int) -> None:
        self.atom_layer = layer
    def get_atom_state(self,state:int) -> None:
        self.atom_state = state   #0-flyout, 1-in solid position，2-in sia state
    
    def get_lattice_order(self,order:int) -> None:             
        self.lattice_order = order  
    def del_lattice_order(self) -> None:
        del self.lattice_order  

    def get_point_order(self,order:int) -> None:
        self.point_order = order
    def del_point_order(self) -> None:
        del self.point_order

    def get_neighbor_point_orders(self,orders:list) -> None:
        self.neighbor_point_orders = orders
    def del_neighbor_point_orders(self) -> None:
        if self.neighbor_point_orders:
            del self.neighbor_point_orders

class ProjectileIon(Atom):
    def get_land_point(self,z_axis) -> None:
        t = abs((self.atom_position[2]-z_axis)/self.atom_velocity[2])
        self.ion_initial_land_point = [x0+t*v0 for x0,v0 in zip(self.atom_position,self.atom_velocity)]