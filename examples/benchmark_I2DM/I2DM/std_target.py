from std_atom import *
class LatticeCell:
    def __init__(self) -> None:
        self.solid_point_in_lattice = []    #solid point order,int
        self.atom_in_lattice = []           #atom order,int
        self.sia_in_lattice = []            #atom order,int
        self.vac_in_lattice = []            #solid point order,int

class SolidPoint:
    def __init__(self,name,position,layer) -> None:
        self.point_name = name
        self.point_position = position
        self.point_layer = layer

class TargetCrystal:
    def __init__(self) -> None:
        self.unit = []                    #TargetAtom
        self.crystal = []                 #TargetAtom
        self.lattice = []                 #LatticeCell
        self.solid_points = []            #SolidPoint
        #self.crystal_copy = []
        #self.lattice_copy = []

        self.displaced_atoms = []         #atom order,after defect recommdation
        self.dis_number = 0
        self.displaced_atoms_total = []
        self.displaced_atom_num = dict()
        self.displaced_atom_data = dict()

        self.vacancies = []               #point order,int
        self.vac_number = 0
        self.primary_vac_number = 0
        self.primary_vac_num = {'SV':{},'DV':{},'MV':{}}
        self.vac_of_atom_type = dict()
        self.vac_num = {'SV':{},'DV':{},'MV':{}}
        self.vac_data = {'SV':{},'DV':{},'MV':{}}

        self.sias = []                    #atom order,int
        self.sia_number = 0
        self.sia_num = {'SSIA':{},'DSIA':{},'MSIA':{}}
        self.sia_data = {'SSIA':{},'DSIA':{},'MSIA':{}}

        self.subs_atoms = []              #atom order,int
        self.subs_number = 0
        self.subs_atom_num = dict()
        self.subs_atom_data = dict()

        self.out_atoms = []               #atom order,int
        self.out_number = 0
        self.out_atom_num = dict()
        self.out_atom_data = dict()

    def add_num_on_displaced_atom_num(self,dis_atom_dict:dict) -> None:
        for key,value in dis_atom_dict.items():
            if key not in self.displaced_atom_num.keys():
                self.displaced_atom_num[key] = value
            else:
                self.displaced_atom_num[key] += value

    def add_num_on_primary_vac_num(self,vac_dict:dict) -> None:
        for key1,value in vac_dict.items():
            for key2 in value.keys():
                if key2 not in self.primary_vac_num[key1].keys():
                    self.primary_vac_num[key1][key2] = value[key2]
                else:
                    self.primary_vac_num[key1][key2] += value[key2]
    def add_num_on_vac_of_atom_type(self,vac_of_atom_type:dict) -> None:
        for key,value in vac_of_atom_type.items():
            if key not in self.vac_of_atom_type.keys():
                self.vac_of_atom_type[key] = value
            else:
                self.vac_of_atom_type[key] += value
    def add_num_on_vac_num(self,vac_dict:dict) -> None:
        for key1,value in vac_dict.items():
            for key2 in value.keys():
                if key2 not in self.vac_num[key1].keys():
                    self.vac_num[key1][key2] = value[key2]
                else:
                    self.vac_num[key1][key2] += value[key2]

    def add_num_on_sia_num(self,sia_dict:dict) -> None:
        for key1 in sia_dict.keys():
            for key2 in sia_dict[key1].keys():
                if key2 not in self.sia_num[key1].keys():
                    self.sia_num[key1][key2] = sia_dict[key1][key2]
                else:
                    self.sia_num[key1][key2] += sia_dict[key1][key2]
    
    def add_num_on_subs_atom_num(self,subs_atom_dict:dict) -> None:
        for key,value in subs_atom_dict.items():
            if key not in self.subs_atom_num.keys():
                self.subs_atom_num[key] = value
            else:
                self.subs_atom_num[key] += value
    
    def add_num_on_out_atom_num(self,out_atom_dict:dict) -> None:
        for key,value in out_atom_dict.items():
            if key not in self.out_atom_num.keys():
                self.out_atom_num[key] = value
            else:
                self.out_atom_num[key] += value
    
    ############################
    ##input target information##
    ##1.structure of target   ##
    ##2.threshold energy      ##
    ############################
    def read_structure_file(self,file_path:str,file_type:int) -> None:
        if file_type == 0:
            file = open(file_path,'r')
            self.system_name = file.readline().strip()

            scale = float(file.readline().strip())
            lattice_x = file.readline().strip().split()            
            lattice_y = file.readline().strip().split()
            lattice_z = file.readline().strip().split()
            vector_x = [float(x)*scale for x in lattice_x]
            vector_y = [float(y)*scale for y in lattice_y]
            vector_z = [float(z)*scale for z in lattice_z]
            self.lattice_vector_matrix = [vector_x,vector_y,vector_z]
            self.lattice_vector_matrix_inverse = inverse_matrix(self.lattice_vector_matrix)

            atom_name_in_unit = file.readline().strip().split()
            atom_num_in_unit = [int(x) for x in file.readline().strip().split()]
            atom_total_num_in_unit = sum(atom_num_in_unit)

            coordinate = file.readline()
            if coordinate[0]=='S' or coordinate[0]=='s':#Selective Dynamics,read the next line
                coordinate = file.readline()
            atom_type_index = 0
            atom_num_judgement = atom_num_in_unit[0]
            atom_order_index = 0
            while atom_order_index < atom_total_num_in_unit:
                if coordinate[0]=='D' or coordinate[0]=='d':#Reduced position
                    reduced_position = [float(x) for x in file.readline().strip().split()[:3]]
                    atom_position = vector_multi_matrix(reduced_position,self.lattice_vector_matrix)
                else:#real position
                    atom_position = [float(x) for x in file.readline().strip().split()[:3]]
                
                if atom_order_index >= atom_num_judgement:
                    atom_type_index += 1
                    atom_num_judgement += atom_num_in_unit[atom_type_index]                    
                atom_name = atom_name_in_unit[atom_type_index]
                index_in_elements = ALL_ELEMENTS_NAME.index(atom_name)
                atom_charge = ALL_ELEMENTS_CHARGE[index_in_elements]
                atom_mass = ALL_ELEMENTS_MASS[index_in_elements]
                atom_radius = ALL_ELEMENTS_RADIUS[index_in_elements]
                atom_alpha = ALL_ELEMENTS_ALPHA[index_in_elements]
                atom_beta = ALL_ELEMENTS_BETA[index_in_elements]
                atom_energy = 0.0
                atom_velocity = [0.0,0.0,0.0]
                new_atom = TargetAtom(atom_name,
                                      atom_charge,
                                      atom_mass,
                                      atom_radius,
                                      atom_position,
                                      atom_energy,
                                      atom_velocity,
                                      atom_alpha,
                                      atom_beta)
                self.unit.append(new_atom)
                atom_order_index += 1
            file.close()
        elif file_type == 1:
            file = open(file_path,'r')
            contents = file.readlines()
            file.close()
            useful_information=[]
            contents = [content.strip().split() for content in contents]
            for content in contents:
                if '_cell_length_a' in content:
                    la_a = float(content[1])
                elif '_cell_length_b' in content:
                    la_b = float(content[1])
                elif '_cell_length_c' in content:
                    la_c = float(content[1])
                elif '_cell_angle_alpha' in content:
                    alpha = float(content[1])*math.pi/180
                elif '_cell_angle_beta' in content:
                    beta = float(content[1])*math.pi/180
                elif '_cell_angle_gamma' in content:
                    gamma = float(content[1])*math.pi/180
                    break
            str_atom_site = []
            for i in range(len(contents)):
                j= len(contents)-i-1
                if len(contents[j])>1:
                    useful_information.insert(0,contents[j])
                else:
                    if contents[j][0].find('_atom_site_')!=-1:
                        str_atom_site.insert(0,contents[j][0])
                    else:
                        break
            index_name = str_atom_site.index('_atom_site_type_symbol')
            index_x = str_atom_site.index('_atom_site_fract_x')
            index_y = str_atom_site.index('_atom_site_fract_y')
            index_z = str_atom_site.index('_atom_site_fract_z')
            vector1 = [la_a,0,0]
            vector2 = [la_b*math.cos(gamma),la_b*math.sin(gamma),0]
            vector3 = [la_c*math.cos(beta),
                       la_c*(math.cos(alpha)-math.cos(beta)*math.cos(gamma))/math.sin(gamma),
                       la_c/math.sin(gamma)*math.sqrt(1+2*math.cos(alpha)*math.cos(beta)*math.cos(gamma)-math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2)]
            self.lattice_vector_matrix = [vector1,vector2,vector3]
            self.lattice_vector_matrix_inverse = self.lattice_matrix**-1
            atom_names = [value[index_name] for value in useful_information]
            atom_positions = [vector_multi_matrix([float(value[index_x]),float(value[index_y]),float(value[index_z])],self.lattice_vector_matrix) for value in useful_information]
            i = 0
            while i < len(atom_names):
                atom_name = atom_names[i]
                atom_position = atom_positions[i]
                index_in_element = ALL_ELEMENTS_NAME.index(atom_name)
                atom_charge =  ALL_ELEMENTS_CHARGE[index_in_element]
                atom_mass = ALL_ELEMENTS_MASS[index_in_element]
                atom_radius = ALL_ELEMENTS_RADIUS[index_in_element]
                atom_alpha = ALL_ELEMENTS_ALPHA[index_in_element]
                atom_beta = ALL_ELEMENTS_BETA[index_in_element]
                atom_energy = 0.0
                atom_velocity = [0.0,0.0,0.0]
                new_atom = TargetAtom(atom_name,
                                      atom_charge,
                                      atom_mass,
                                      atom_radius,
                                      atom_position,
                                      atom_energy,
                                      atom_velocity,
                                      atom_alpha,
                                      atom_beta)
                self.unit.append(new_atom)
                i+=1
            atom_name_statics = dict()
            for name in atom_names:
                if name not in atom_name_statics.keys():
                    atom_name_statics[name]=1
                else:
                    atom_name_statics[name]+=1
            self.system_name= ''
            if len(atom_name_statics)==1:
                for key in atom_name_statics.keys():
                    self.system_name += f'{key}'
            else:
                for key,value in atom_name_statics.items():
                    if atom_name_statics[key]==1:
                        self.system_name += f'{key}'
                    else:
                        self.system_name += f'{key}{value}'
        else:
            raise NameError("File Error!")
                  
    def read_threshold_energy_file(self,file_path:str) -> None:
        file = open(file_path,'r')
        lines = file.readlines()
        file.close()
        i=0
        while i < len(lines):
            if lines[i][0]=='#':
                lines.remove(lines[i])
            else:
                i+=1
        for i in range(len(lines)):
            line = lines[i].strip().split()
            self.unit[i].get_threshold_energy(float(line[1]))
    ###########################
    ##   create the target   ##
    ##1.get target size      ##
    def get_target_size(self,X:int,Y:int,Z:int) -> None:
        self.size_x = X
        self.size_y = Y
        self.size_z = Z
        self.layer_z = dict()
        if len(self.unit)==1:
            self.unit[0].atom_layer(1)
            self.minus_z_axis = self.unit[0].atom_position[2]
            self.layer_z[1] = self.minus_z_axis
            sort = 2
        else:          
            z_aixs = []
            for atom in self.unit:
                z_aixs.append(atom.atom_position[2])
            layer_order = [0]*len(z_aixs)
            z_aixs_copy = [x for x in z_aixs]
            self.minus_z_axis = max(z_aixs_copy)
            sort = 1
            while len(z_aixs_copy)>0:
                MAX =max(z_aixs_copy)
                self.layer_z[sort] = MAX
                i = 0
                while i < len(z_aixs):
                    if z_aixs[i] == MAX:
                        if layer_order[i] ==0:
                            layer_order[i]= sort
                    i+=1
                while z_aixs_copy.count(MAX)>0:
                    z_aixs_copy.remove(MAX)
                sort += 1
            i = 0
            while i < len(self.unit):
                self.unit[i].get_atom_layer(layer_order[i])
                i+=1
        self.layer_num = sort - 1
        self.lattice_num = self.size_x*self.size_y*self.size_z
        self.atom_num = self.lattice_num*len(self.unit)
        la_x = self.lattice_vector_matrix[0]
        la_y = self.lattice_vector_matrix[1]
        la_z = self.lattice_vector_matrix[2]
        v3 = [x+y for x,y in zip(la_x,la_y)]
        new_matrix = matrix_trans([la_x,la_y,v3])
        l_x = max(new_matrix[0])-min([min(new_matrix[0]),0])
        l_y = max(new_matrix[1])-min([min(new_matrix[1]),0])
        self.l_z = max(la_z)-min(la_z)
        for i in range(sort-1):
            self.layer_z[i+1] -= ((self.size_z-1)*self.l_z+ self.minus_z_axis)
        self.range_x = [-l_x*self.size_x/2,l_x*self.size_x/2]
        self.range_y = [-l_y*self.size_y/2,l_y*self.size_y/2]
        self.range_z = [-self.l_z*(self.size_z-1)-self.minus_z_axis,self.l_z-self.minus_z_axis]
        self.sensitive_radius = min(self.size_x*l_x,self.size_y*l_y)*0.45
        
        self.lattice_parameter = max(max(la_x),max(la_y))
        self.source_height = 2*self.range_x[1]*self.lattice_parameter
        self.cut_radius = self.lattice_parameter
        self.unit_area = get_norm_of_vector(vector_cross(self.lattice_vector_matrix[0],self.lattice_vector_matrix[1]))
        self.capture_radius = self.lattice_parameter
        self.N = len(self.unit)/(self.layer_num* self.unit_area)
    
    def get_conditions(self,emitter_parameters:dict, ion_parameters:dict):
        self.conditions = dict()
        self.conditions['SYSTEM'] = self.system_name
        self.conditions['SIZE_X'] = self.size_x
        self.conditions['SIZE_Y'] = self.size_y
        self.conditions['SIZE_Z'] = self.size_z
        self.conditions['ION_TIME'] = ion_parameters['time']
        self.conditions['ION_NAME'] = ion_parameters['name']
        self.conditions['ION_CHARGE'] = ion_parameters['charge']
        self.conditions['ION_MASS'] = ion_parameters['mass']
        self.conditions['ION_ENERGY'] = ion_parameters['energy']
        self.conditions['ION_SANGLE'] = ion_parameters['sangle']
        self.conditions['ION_PANGLE'] = ion_parameters['pangle']
        self.conditions['UNIT_ATOMS'] = len(self.unit)
        self.conditions['INITIAL_ALL_ATOMS'] = self.atom_num
        self.conditions['INITIAL_IRRAD_ATOMS'] = int(len(self.unit)*emitter_parameters['irradiation_area']/self.unit_area)

    def create_target(self) -> None:
        atom_count = 0
        lattice_count = 0
        for k in range(self.size_z):
            for j in range(self.size_y):
                for i in range(self.size_x):
                    lattice_position = vector_multi_matrix([i-self.size_x/2,j-self.size_y/2,k],self.lattice_vector_matrix)
                    lattice_position[2] = lattice_position[2]-(self.size_z-1)*self.l_z-self.minus_z_axis
                    new_lattice = LatticeCell()
                    for atom in self.unit:
                        atom_position = [x+y for x,y in zip(atom.atom_position,lattice_position)]
                        new_atom = TargetAtom(atom.atom_name,
                                              atom.atom_charge,
                                              atom.atom_mass,
                                              atom.atom_radius,
                                              atom_position,
                                              atom.atom_energy,
                                              atom.atom_velocity,
                                              atom.atom_alpha,
                                              atom.atom_beta)
                        #properties for new atom
                        new_atom.get_threshold_energy(atom.Td)
                        new_atom.get_atom_layer(atom.atom_layer+self.layer_num*(self.size_z-k-1))
                        new_atom.get_atom_state(1)
                        new_atom.get_lattice_order(lattice_count)
                        new_atom.get_point_order(atom_count)
                        self.crystal.append(new_atom)
                        #information for new lattice
                        new_lattice.solid_point_in_lattice.append(atom_count)
                        new_lattice.atom_in_lattice.append(atom_count)
                        #information for target
                        new_solid_point = SolidPoint(atom.atom_name,atom_position,atom.atom_layer+self.layer_num*k)
                        self.solid_points.append(new_solid_point)
                        atom_count += 1
                    self.lattice.append(new_lattice)
                    lattice_count += 1
        #self.crystal_copy = copy.deepcopy(self.crystal)
        #self.lattice_copy = copy.deepcopy(self.lattice)
    def refresh_system(self):
        for displaced_atom in self.displaced_atoms_total:
            if self.crystal[displaced_atom].atom_state ==0:
                lattice_order = int(displaced_atom/len(self.unit))
                self.lattice[lattice_order].atom_in_lattice.append(displaced_atom)
                self.lattice[lattice_order].vac_in_lattice.clear()
                self.lattice[lattice_order].sia_in_lattice.clear()
                self.crystal[displaced_atom].get_new_energy(0)
                self.crystal[displaced_atom].get_new_position(self.solid_points[displaced_atom].point_position)
                self.crystal[displaced_atom].get_new_velocity([0,0,0])
                self.crystal[displaced_atom].get_atom_state(1)
                self.crystal[displaced_atom].get_lattice_order(int(displaced_atom/len(self.unit)))
                self.crystal[displaced_atom].get_atom_layer(self.solid_points[displaced_atom].point_layer)
                self.crystal[displaced_atom].get_point_order(displaced_atom)
                self.crystal[displaced_atom].del_neighbor_atom_orders()
            elif self.crystal[displaced_atom].atom_state ==1:
                lattice_before = self.crystal[displaced_atom].lattice_order
                lattice_now = int(displaced_atom/len(self.unit))
                self.lattice[lattice_before].sia_in_lattice.clear()
                self.lattice[lattice_before].vac_in_lattice.clear() 
                self.lattice[lattice_before].atom_in_lattice.remove(displaced_atom)
                self.lattice[lattice_now].sia_in_lattice.clear()
                self.lattice[lattice_now].vac_in_lattice.clear() 
                self.lattice[lattice_now].atom_in_lattice.append(displaced_atom)

                self.crystal[displaced_atom].get_new_energy(0)
                self.crystal[displaced_atom].get_new_position(self.solid_points[displaced_atom].point_position)
                self.crystal[displaced_atom].get_new_velocity([0,0,0])
                self.crystal[displaced_atom].get_lattice_order(lattice_now)
                self.crystal[displaced_atom].get_atom_layer(self.solid_points[displaced_atom].point_layer)
                self.crystal[displaced_atom].get_point_order(displaced_atom)
                self.crystal[displaced_atom].del_neighbor_atom_orders()
            else:
                lattice_before = self.crystal[displaced_atom].lattice_order
                lattice_now = int(displaced_atom/len(self.unit))
                self.lattice[lattice_before].sia_in_lattice.clear()
                self.lattice[lattice_before].vac_in_lattice.clear() 
                self.lattice[lattice_before].atom_in_lattice.remove(displaced_atom)
                self.lattice[lattice_now].sia_in_lattice.clear()
                self.lattice[lattice_now].vac_in_lattice.clear() 
                self.lattice[lattice_now].atom_in_lattice.append(displaced_atom)
                self.crystal[displaced_atom].get_new_energy(0)
                self.crystal[displaced_atom].get_new_position(self.solid_points[displaced_atom].point_position)
                self.crystal[displaced_atom].get_new_velocity([0,0,0])
                self.crystal[displaced_atom].get_lattice_order(lattice_now)
                self.crystal[displaced_atom].get_point_order(displaced_atom)
                self.crystal[displaced_atom].del_neighbor_atom_orders()
                self.crystal[displaced_atom].get_atom_layer(self.solid_points[displaced_atom].point_layer)
                self.crystal[displaced_atom].get_atom_state(1)
        while len(self.crystal)> self.atom_num:
            atom_lattice = self.crystal[len(self.crystal)-1].lattice_order
            if (len(self.crystal)-1) in self.lattice[atom_lattice].atom_in_lattice:
                self.lattice[atom_lattice].atom_in_lattice.remove(len(self.crystal)-1)
            self.lattice[atom_lattice].sia_in_lattice.clear()
            self.lattice[atom_lattice].vac_in_lattice.clear()   
            del self.crystal[-1]
        '''
        self.crystal = copy.deepcopy(self.crystal_copy)
        self.lattice = copy.deepcopy(self.lattice_copy)'''
        self.displaced_atoms_total.clear()
        self.displaced_atoms.clear()
        self.vacancies.clear()
        self.sias.clear()
        self.subs_atoms.clear()
        self.out_atoms.clear()
    #find neighbor items
    def find_neighbor_atoms_of_given_position(self,neighbor_width:int,position:list) -> list:
        whether_in_crystal = self.check_whether_position_in_crystal(position)
        if whether_in_crystal[0]:      
            lattice_order_now = whether_in_crystal[1]
            neighbor_atom_orders = []
            for j in range(-neighbor_width,neighbor_width+1):
                for i in range(-neighbor_width,neighbor_width+1):
                    neighbor_lattice_order = lattice_order_now + i + j*self.size_x
                    if 0 <= neighbor_lattice_order < self.lattice_num:
                        for order in self.lattice[neighbor_lattice_order].atom_in_lattice:
                            distance = dis_bet_points(self.crystal[order].atom_position,position)
                            if distance < self.cut_radius*neighbor_width and distance != 0:
                                neighbor_atom_orders.append(order)
                                
            return [True,neighbor_atom_orders]
        else:
            return [False]
    def find_neighbor_solid_points_of_given_position(self,neighbor_width:int,position:list) -> list:
        whether_in_crystal = self.check_whether_position_in_crystal(position)
        if whether_in_crystal[0]:      
            lattice_order_now = whether_in_crystal[1]
            neighbor_point_orders = []
            for j in range(-neighbor_width,neighbor_width+1):
                for i in range(-neighbor_width,neighbor_width+1):
                    neighbor_lattice_order = lattice_order_now + i + j*self.size_x
                    if 0 <= neighbor_lattice_order < self.lattice_num:
                        for order in self.lattice[neighbor_lattice_order].solid_point_in_lattice:
                            distance = dis_bet_points(self.solid_points[order].point_position,position)
                            if distance < self.cut_radius*neighbor_width:
                                neighbor_point_orders.append(order)
            return [True,neighbor_point_orders]
        else:
            return [False]
    #produce or delete defect
    def check_whether_position_in_crystal(self,position:list) -> list:
        displacement = [self.size_x/2,self.size_y/2,self.size_z/2]
        old_reduced_position = vector_multi_matrix(position,self.lattice_vector_matrix_inverse)
        reduced_position = [x+y for x,y in zip(old_reduced_position,displacement)]
        reduced_position_int = [int(x//1) for x in reduced_position]
        whether_in = True
        if reduced_position_int[0]<0 or reduced_position_int[0] >= self.size_x:
            whether_in = False
        if reduced_position_int[1]<0 or reduced_position_int[1] >= self.size_y:
            whether_in = False
        if reduced_position_int[2]<0 or reduced_position_int[2] >= self.size_z:
            whether_in = False
        if whether_in:
            lattice_order_now = (reduced_position_int[2]*self.size_x*self.size_y + 
                                reduced_position_int[1]*self.size_x + 
                                reduced_position_int[0])
            return [whether_in,lattice_order_now]
        else:
            return [whether_in]
    def whether_can_atom_cascade(self,atom_order:int) -> bool:
        information = self.check_whether_position_in_crystal(self.crystal[atom_order].atom_position)
        whether_in_crystal = information[0]        
        if self.crystal[atom_order].atom_state == 1:
            whether_energy_above_Td = True if (self.crystal[atom_order].atom_energy > self.crystal[atom_order].Td) else False
            #whether_not_displaced = True if (dis_bet_points(self.crystal[atom_order].atom_position,self.solid_points[self.crystal[atom_order].point_order].point_position)
                #<self.crystal[atom_order].atom_radius/2) else False
            if whether_in_crystal:
                lattice_order_before = self.crystal[atom_order].lattice_order
                lattice_order_now = information[1]
                '''known states:
                (1)atom_state = 1
                (2)whether_in_crystal = True
                '''
                if whether_energy_above_Td:
                    #Can cascade
                    self.crystal[atom_order].get_new_energy(self.crystal[atom_order].atom_energy - self.crystal[atom_order].Td)
                    self.sias.append(atom_order)
                    self.vacancies.append(self.crystal[atom_order].point_order)
                    if atom_order in self.subs_atoms:
                        self.subs_atoms.remove(atom_order)
                    self.lattice[lattice_order_before].vac_in_lattice.append(self.crystal[atom_order].point_order)
                    self.lattice[lattice_order_now].sia_in_lattice.append(atom_order)
                    self.crystal[atom_order].get_atom_state(2)
                    self.crystal[atom_order].del_point_order()
                    if lattice_order_now != lattice_order_before:
                        self.lattice[lattice_order_before].atom_in_lattice.remove(atom_order)
                        self.lattice[lattice_order_now].atom_in_lattice.append(atom_order)
                        self.crystal[atom_order].get_lattice_order(lattice_order_now)
                    return True
                else:
                    #cannot cascade
                    self.crystal[atom_order].get_new_energy(0.0)
                    self.crystal[atom_order].get_new_velocity([0,0,0])
                    self.crystal[atom_order].get_new_position(self.solid_points[self.crystal[atom_order].point_order].point_position)
                    return False
            else:
                '''known states:
                (1)atom_state = 1
                (2)whether_in_crystal = False
                (3)whether_in_original_lattice = False
                '''
                lattice_order_before = self.crystal[atom_order].lattice_order
                if whether_energy_above_Td:
                    self.atom_fly_out(atom_order)
                else:
                    self.crystal[atom_order].get_new_energy(0.0)
                    self.crystal[atom_order].get_new_velocity([0,0,0])
                    self.crystal[atom_order].get_new_position(self.solid_points[self.crystal[atom_order].point_order].point_position)
                return False
        elif self.crystal[atom_order].atom_state == 2:
            lattice_order_before = self.crystal[atom_order].lattice_order
            if whether_in_crystal:
                lattice_order_now = information[1]
                if lattice_order_now != lattice_order_before:
                    self.lattice[lattice_order_before].atom_in_lattice.remove(atom_order)
                    self.lattice[lattice_order_before].sia_in_lattice.remove(atom_order)
                    self.lattice[lattice_order_now].atom_in_lattice.append(atom_order)
                    self.lattice[lattice_order_now].sia_in_lattice.append(atom_order)
                    self.crystal[atom_order].get_lattice_order(lattice_order_now)
                return True
            else:
                self.out_atoms.append(atom_order)
                self.sias.remove(atom_order)
                self.lattice[lattice_order_before].atom_in_lattice.remove(atom_order)
                self.lattice[lattice_order_before].sia_in_lattice.remove(atom_order)
                self.crystal[atom_order].get_atom_state(0)
                self.crystal[atom_order].del_lattice_order()
                return False            
    def get_new_atom_from_ion(self,ion:ProjectileIon,lattice_order:int) -> None:
        new_atom = TargetAtom(ion.atom_name,
                              ion.atom_charge,
                              ion.atom_mass,
                              ion.atom_radius,
                              ion.atom_position,
                              0.0,
                              [0.0,0.0,0.0],
                              ion.atom_alpha,
                              ion.atom_beta)
        new_atom.get_threshold_energy(5.0)
        dis = []
        layers = []
        for k in range(self.size_z):
            z = (1-self.size_z+k)*self.l_z-self.minus_z_axis
            for atom in self.unit:
                distance = abs(ion.atom_position[2]-(atom.atom_position[2]+z))
                layer = atom.atom_layer +(self.size_z-k-1)*self.layer_num
                layers.append(layer)
                dis.append(distance)
        MIN_value = min(dis)
        index = dis.index(MIN_value)
        new_atom.get_atom_layer(layers[index])
        new_atom.get_atom_state(2)
        new_atom.get_lattice_order(lattice_order)
        atom_order = len(self.crystal)
        self.crystal.append(new_atom)
        self.sias.append(atom_order)
        self.lattice[lattice_order].atom_in_lattice.append(atom_order)
        self.lattice[lattice_order].sia_in_lattice.append(atom_order)
    def atom_fly_out(self,atom_order:int) -> None:
        lattice_order_before = self.crystal[atom_order].lattice_order
        if self.crystal[atom_order].atom_state == 1:
            self.vacancies.append(self.crystal[atom_order].point_order)
            self.lattice[lattice_order_before].vac_in_lattice.append(self.crystal[atom_order].point_order)
            if atom_order in self.subs_atoms:
                self.subs_atoms.remove(atom_order)
        elif self.crystal[atom_order].atom_state == 2:
            if atom_order in self.sias:
                self.sias.remove(atom_order)
            self.lattice[lattice_order_before].sia_in_lattice.remove(atom_order)
        self.out_atoms.append(atom_order)
        self.lattice[lattice_order_before].atom_in_lattice.remove(atom_order)
        self.crystal[atom_order].get_atom_state(0)
    ###########################
    ###structure relaxation####
    ###########################
    def defect_recombination(self,width:float) -> None:
        i=0
        while i < len(self.sias):
            vacancy_orders = []
            distances = []
            neighbor_point_orders = self.find_neighbor_solid_points_of_given_position(1,self.crystal[self.sias[i]].atom_position)
            for point_order in neighbor_point_orders[1]:
                if point_order in self.vacancies:
                    distance = dis_bet_points(self.solid_points[point_order].point_position,self.crystal[self.sias[i]].atom_position)
                    layer_dis = []
                    for value in self.layer_z.values():
                        layer_dis.append(abs(self.crystal[self.sias[i]].atom_position[2]-value))
                    min_layer_dis = min(layer_dis)
                    diff_layer = abs(layer_dis.index(min_layer_dis) + 1 - self.solid_points[point_order].point_layer)
                    if distance <= self.capture_radius*width and diff_layer <= 1:
                        distances.append(distance)
                        vacancy_orders.append(point_order)
            if len(vacancy_orders) == 0:
                i += 1
            else:
                MIN = min(distances)
                index = distances.index(MIN)
                selected_sia = self.sias[i]
                selected_vac = vacancy_orders[index]
                lattice_order_before = self.crystal[selected_sia].lattice_order
                lattice_order_now = int(selected_vac//len(self.unit))
                self.sias.remove(selected_sia)
                self.vacancies.remove(selected_vac)
                if self.solid_points[selected_vac].point_name != self.crystal[selected_sia].atom_name:
                    self.subs_atoms.append(selected_sia)
                self.lattice[lattice_order_before].atom_in_lattice.remove(selected_sia)
                self.lattice[lattice_order_before].sia_in_lattice.remove(selected_sia)
                self.lattice[lattice_order_now].vac_in_lattice.remove(selected_vac)
                self.lattice[lattice_order_now].atom_in_lattice.append(selected_sia)
                self.crystal[selected_sia].get_atom_state(1)
                self.crystal[selected_sia].get_new_position(self.solid_points[selected_vac].point_position)
                self.crystal[selected_sia].get_lattice_order(lattice_order_now)
                self.crystal[selected_sia].get_point_order(selected_vac)
    ###########################
    ########analyse data#######
    ###########################
    #VACANCY
    def get_vac_of_atom_type(self,point_list:list) -> dict:
        vac_of_atom_type = dict()
        for point_order in point_list:
            name = f'{self.solid_points[point_order].point_name}{self.solid_points[point_order].point_layer}'
            if name in vac_of_atom_type.keys():
                vac_of_atom_type[name]+=1
            else:
                vac_of_atom_type[name]=1
        return vac_of_atom_type
    def devide_vacancy_list(self,point_list:list) -> list:
        vacancy_cluster_list = []
        list_tmp = [x for x in point_list]
        while len(list_tmp)>0:
            cluster = []
            cluster.append(list_tmp[0])
            del list_tmp[0]
            i = 0
            while i<len(cluster):
                j = 0
                while j<len(list_tmp):
                    if ((abs(self.solid_points[list_tmp[j]].point_layer
                             -self.solid_points[cluster[i]].point_layer)<=2)
                         and (dis_bet_points(self.solid_points[list_tmp[j]].point_position,
                                             self.solid_points[cluster[i]].point_position)
                              <self.capture_radius)):
                        cluster.append(list_tmp[j])
                        del list_tmp[j]
                    else:
                        j+=1
                i+=1
            vacancy_cluster_list.append(cluster)
        return vacancy_cluster_list
    def get_vacancy_type_num_dict(self,cluster_list:list) -> dict:
        vac_num = {'SV':{},'DV':{},'MV':{}}
        i = 0
        while i < len(cluster_list):
            if len(cluster_list[i]) ==1:
                name0 = 'SV'
                solid_order = cluster_list[i][0]
                name_i = self.solid_points[solid_order].point_name
                layer_i = self.solid_points[solid_order].point_layer
                name1 = f'V_{name_i}{layer_i}'
            elif len(cluster_list[i]) ==2:
                name0 = 'DV'
                solid_order1 = cluster_list[i][0]
                solid_order2 = cluster_list[i][1]
                name_i1 = self.solid_points[solid_order1].point_name
                name_i2 = self.solid_points[solid_order2].point_name
                layer_i1 = self.solid_points[solid_order1].point_layer
                layer_i2 = self.solid_points[solid_order2].point_layer
                name1 = ((f'V_{name_i1}{layer_i1}+V_{name_i2}{layer_i2}') 
                         if (layer_i1 <= layer_i2) else 
                         (f'V_{name_i2}{layer_i2}+V_{name_i1}{layer_i1}'))
            else:
                name0 = 'MV'
                name1 = 'MV'
            if name1 not in vac_num[name0].keys():
                vac_num[name0][name1] = 1
            else:
                vac_num[name0][name1] += 1
            i+=1
        return vac_num
    def get_vacancy_type_list_dict(self,cluster_list:list) -> dict:
        vac_data = {'SV':{},'DV':{},'MV':{}}
        i = 0
        while i < len(cluster_list):
            if len(cluster_list[i]) ==1:
                name0 = 'SV'
                solid_order = cluster_list[i][0]
                name_i = self.solid_points[solid_order].point_name
                layer_i = self.solid_points[solid_order].point_layer
                name1 = f'V_{name_i}{layer_i}'
            elif len(cluster_list[i]) ==2:
                name0 = 'DV'
                solid_order1 = cluster_list[i][0]
                solid_order2 = cluster_list[i][1]
                name_i1 = self.solid_points[solid_order1].point_name
                name_i2 = self.solid_points[solid_order2].point_name
                layer_i1 = self.solid_points[solid_order1].point_layer
                layer_i2 = self.solid_points[solid_order2].point_layer
                name1 = ((f'V_{name_i1}{layer_i1}+V_{name_i2}{layer_i2}') 
                         if (layer_i1 <= layer_i2) else 
                         (f'V_{name_i2}{layer_i2}+V_{name_i1}{layer_i1}'))
            else:
                name0 = 'MV'
                name1 = 'MV'
            if name1 not in vac_data[name0].keys():
                vac_data[name0][name1] = cluster_list[i]
            else:
                vac_data[name0][name1].append(cluster_list[i])
            i+=1
        return vac_data
    def get_primary_vacancy_data(self,order_list) -> None:
        cluster_list = self.devide_vacancy_list(order_list)
        vac_dict_result = self.get_vacancy_type_num_dict(cluster_list)
        self.add_num_on_primary_vac_num(vac_dict_result)
        self.primary_vac_number += len(order_list)    
    def get_total_displaced_atoms(self,order_list) -> None:
        self.displaced_atoms_total = order_list
    #SIA
    def devide_sia_list(self,sia_list)->list:
        sia_cluster_list = []
        list_tmp = [x for x in sia_list]
        while len(list_tmp)>0:
            cluster = []
            cluster.append(list_tmp[0])
            del list_tmp[0]
            i = 0
            while i<len(cluster):
                j = 0
                while j<len(list_tmp):
                    if ((abs(self.crystal[list_tmp[j]].atom_layer
                             -self.crystal[cluster[i]].atom_layer)<=1)
                         and (dis_bet_points(self.crystal[list_tmp[j]].atom_position,
                                             self.crystal[cluster[i]].atom_position)
                              <self.capture_radius)):
                        cluster.append(list_tmp[j])
                        del list_tmp[j]
                    else:
                        j+=1
                i+=1
            sia_cluster_list.append(cluster)
        return sia_cluster_list
    def get_sia_type_num_dict(self,cluster_list:list) -> dict:
        sia_num = {'SSIA':{},'DSIA':{},'MSIA':{}}
        i = 0
        while i < len(cluster_list):
            if len(cluster_list[i]) ==1:
                name0 = 'SSIA'
                atom_order = cluster_list[i][0]
                name_i = self.crystal[atom_order].atom_name
                layer_i = self.crystal[atom_order].atom_layer
                name1 = f'I_{name_i}{layer_i}'
            elif len(cluster_list[i]) ==2:
                name0 = 'DSIA'
                atom_order1 = cluster_list[i][0]
                atom_order2 = cluster_list[i][1]
                name_i1 = self.crystal[atom_order1].atom_name
                name_i2 = self.crystal[atom_order2].atom_name
                layer_i1 = self.crystal[atom_order1].atom_layer
                layer_i2 = self.crystal[atom_order2].atom_layer
                name1 = ((f'I_{name_i1}{layer_i1}+I_{name_i2}{layer_i2}') 
                         if (layer_i1 <= layer_i2) else 
                         (f'I_{name_i2}{layer_i2}+I_{name_i1}{layer_i1}'))
            else:
                name0 = 'MSIA'
                name1 = 'MSIA'
            if name1 not in sia_num[name0].keys():
                sia_num[name0][name1] = 1
            else:
                sia_num[name0][name1] += 1
            i+=1
        return sia_num
    def get_sia_type_list_dict(self,cluster_list:list) -> dict:
        sia_data = {'SSIA':{},'DSIA':{},'MSIA':{}}
        i = 0
        while i < len(cluster_list):
            if len(cluster_list[i]) ==1:
                name0 = 'SSIA'
                atom_order = cluster_list[i][0]
                name_i = self.crystal[atom_order].atom_name
                layer_i = self.crystal[atom_order].atom_layer
                name1 = f'I_{name_i}{layer_i}'
            elif len(cluster_list[i]) ==2:
                name0 = 'DSIA'
                atom_order1 = cluster_list[i][0]
                atom_order2 = cluster_list[i][1]
                name_i1 = self.crystal[atom_order1].atom_name
                name_i2 = self.crystal[atom_order2].atom_name
                layer_i1 = self.crystal[atom_order1].atom_layer
                layer_i2 = self.crystal[atom_order2].atom_layer
                name1 = ((f'I_{name_i1}{layer_i1}+I_{name_i2}{layer_i2}') 
                         if (layer_i1 <= layer_i2) else 
                         (f'I_{name_i2}{layer_i2}+I_{name_i1}{layer_i1}'))
            else:
                name0 = 'MSIA'
                name1 = 'MSIA'
            if name1 not in sia_data[name0].keys():
                sia_data[name0][name1] = [cluster_list[i]]
            else:
                sia_data[name0][name1].append(cluster_list[i])
            i+=1
        return sia_data
    #dispaced and out atom
    def get_atom_type_num_dict(self,dis_atom_list:list) -> dict:
        atom_num = dict()
        for atom_order in dis_atom_list:
            name0 = f'{self.crystal[atom_order].atom_name}{self.crystal[atom_order].atom_layer}'
            if name0 not in atom_num.keys():
                atom_num[name0] = 1
            else:
                atom_num[name0] += 1
        return atom_num
    def get_atom_type_list_dict(self,dis_atom_list:list) -> dict:
        atom_data = dict()
        for atom_order in dis_atom_list:
            name0 = f'{self.crystal[atom_order].atom_name}{self.crystal[atom_order].atom_layer}'
            if name0 not in atom_data.keys():
                atom_data[name0] = [atom_order]
            else:
                atom_data[name0].append(atom_order)
        return atom_data
    def get_catched_atom_num(self) -> int:
        result = 0
        for atom_order in self.out_atoms:
            if atom_order >= self.atom_num:
                result += 1
        return result
    #subs atom
    def get_subs_atom_type_num_dict(self,subs_atom_list:list) -> dict:
        subs_atom_num = dict()
        for atom_order in subs_atom_list:
            atom_name = self.crystal[atom_order].atom_name
            point_name = self.solid_points[self.crystal[atom_order].point_order].point_name
            name = f'{atom_name}:{point_name}'
            if name not in subs_atom_num.keys():
                subs_atom_num[name] = 1 
            else:
                subs_atom_num[name] += 1
        return subs_atom_num
    def get_subs_atom_type_list_dict(self,subs_atom_list:list) -> dict:
        subs_atom_data = dict()
        for atom_order in subs_atom_list:
            atom_name = self.crystal[atom_order].atom_name
            point_name = self.solid_points[self.crystal[atom_order].point_order].point_name
            name = f'{atom_name}:{point_name}'
            if name not in subs_atom_data.keys():
                subs_atom_data[name] = [atom_order]
            else:
                subs_atom_data[name].append(atom_order)
        return subs_atom_data
    #summary of Defect Information
    def get_defect_data(self,MODE:int) -> None:
        #1.vacancy
        self.add_num_on_vac_of_atom_type(self.get_vac_of_atom_type(self.vacancies))
        self.vac_number += len(self.vacancies)
        vac_cluster = self.devide_vacancy_list(self.vacancies)
        vac_dict_result = self.get_vacancy_type_num_dict(vac_cluster)
        self.add_num_on_vac_num(vac_dict_result)
        #2.sia
        self.sia_number += len(self.sias)
        sia_cluster = self.devide_sia_list(self.sias)
        sia_dict_result = self.get_sia_type_num_dict(sia_cluster)
        self.add_num_on_sia_num(sia_dict_result)
        #3.displaced atom
        self.dis_number += len(self.displaced_atoms)
        dis_atom_dict_result = self.get_atom_type_num_dict(self.displaced_atoms)
        self.add_num_on_displaced_atom_num(dis_atom_dict_result)
        #4.subs atom
        self.subs_number += len(self.subs_atoms)
        subs_atom_dict_result = self.get_subs_atom_type_num_dict(self.subs_atoms)
        self.add_num_on_subs_atom_num(subs_atom_dict_result)
        #5.out atom
        self.out_number += len(self.out_atoms)
        out_atom_dict_result = self.get_atom_type_num_dict(self.out_atoms)
        self.add_num_on_out_atom_num(out_atom_dict_result)
        #6.Continuous mode
        if MODE !=0 :
            self.vac_data = self.get_vacancy_type_list_dict(vac_cluster)
            self.sia_data = self.get_sia_type_list_dict(sia_cluster)
            self.displaced_atom_data = self.get_atom_type_list_dict(self.displaced_atoms)
            self.subs_atom_data = self.get_subs_atom_type_list_dict(self.subs_atoms)
            self.out_atom_data = self.get_atom_type_list_dict(self.out_atoms)
            out_catched_atom_num = self.get_catched_atom_num()
            self.conditions['FINAL_ALL_ATOMS'] = len(self.crystal)-len(self.out_atoms)
            self.conditions['CATCH_ATOMS'] = len(self.crystal)-self.conditions['INITIAL_ALL_ATOMS']-out_catched_atom_num
            
    ###########################
    #output target information#
    ###########################
    def write_structure(self,tag:bool,file_path:str,step:int,MODE:int) -> None:
        if tag and MODE!=0:
            file = open(file_path,'a+')
            file.write(f'ITEM: TIMESTEP\n{step}\n')
            file.write(f'ITEM: NUMBER OF ATOMS\n{len(self.crystal)-len(self.out_atoms)}\n')
            file.write(f'ITEM: BOX BOUNDS pp pp pp\n')
            file.write(f'{self.range_x[0]}\t{self.range_x[1]}\n')
            file.write(f'{self.range_y[0]}\t{self.range_y[1]}\n')
            file.write(f'{self.range_z[0]}\t{self.range_z[1]}\n')
            file.write(f'ITEM: ATOMS id type x y z\n')
            i = 0
            while i < len(self.crystal):
                if self.crystal[i].atom_state!= 0:
                    file.write(f'{i}\t{self.crystal[i].atom_charge}\t')
                    file.write(f'{self.crystal[i].atom_position[0]}\t')
                    file.write(f'{self.crystal[i].atom_position[1]}\t')
                    file.write(f'{self.crystal[i].atom_position[2]}\n')
                i += 1
            file.close()
    def write_summary_defect(self,tag:bool,file_path:str,irradiation_number:int,MODE:int) -> None:
        if tag:
            file = open(file_path,'w')
            file.write('########SYSTEM_PARAMETERS##########\n')
            for key,value in self.conditions.items():
                file.write(f'{key}\t{value}\n')
            irradiated_atom_num = self.conditions['INITIAL_IRRAD_ATOMS']
            file.write('########DISPLACED ATOM DATA#########\n')
            file.write(f'TOTAL_DISPLACED_ATOM\t{self.dis_number}\t{self.dis_number/irradiation_number}')
            if MODE != 0:
                file.write(f'\t{self.dis_number/irradiated_atom_num}')
            for name,num in self.displaced_atom_num.items():
                file.write(f'\n{name}\t{num}\t{num/irradiation_number}')
            if MODE == 0:
                file.write(f'\nATOMS SPUTTERED_BY_ION\t{self.primary_vac_number}\t{self.primary_vac_number/irradiation_number}')
                for pda_name,pda_num in self.primary_vac_num.items():
                    num_of_PDA = 0
                    for value in pda_num.values():
                        num_of_PDA += value
                    file.write(f'\n{pda_name}\t{num_of_PDA}\t{num_of_PDA/irradiation_number}')
                    for key,value in pda_num.items():
                        file.write(f'\n{key}\t{value}\t{value/irradiation_number}')

            file.write('\n###########VACANCY DATA############')
            file.write(f'\nVACANCY\t{self.vac_number}\t{self.vac_number/irradiation_number}')
            if MODE != 0:
                file.write(f'\t{self.vac_number/irradiated_atom_num}')
            file.write('\n')
            for name,type_num in self.vac_of_atom_type.items():
                file.write(f'{name}\t{type_num}\t{type_num/irradiation_number}\n')
            for name,type_num in self.vac_num.items():
                num_of_type = 0
                for value in type_num.values():
                    num_of_type+=value
                file.write(f'{name}\t{num_of_type}\t{num_of_type/irradiation_number}\n')
                for key,value in type_num.items():
                    file.write(f'{key}\t{value}\t{value/irradiation_number}\n')
                            
            file.write('###########SIA DATA############\n')
            file.write(f'SIA\t{self.sia_number}\t{self.sia_number/irradiation_number}')
            if MODE != 0:
                file.write(f'\t{self.sia_number/irradiated_atom_num}')
            file.write('\n')
            for name,type_num in self.sia_num.items():
                num_of_type = 0
                for value in type_num.values():
                    num_of_type+=value
                file.write(f'{name}\t{num_of_type}\t{num_of_type/irradiation_number}\n')
                for key,value in type_num.items():
                    file.write(f'{key}\t{value}\t{value/irradiation_number}\n')
                
            file.write('###########SUBSTANCE ATOM DATA############\n')
            file.write(f'SUBS_ATOM\t{self.subs_number}\t{self.subs_number/irradiation_number}')
            if MODE != 0:
                file.write(f'\t{self.subs_number/irradiated_atom_num}')
            file.write('\n')
            for name,num in self.subs_atom_num.items():
                file.write(f'{name}\t{num}\t{num/irradiation_number}\n')
                
            file.write('###########OUT ATOM DATA############\n')
            file.write(f'OUT_ATOM\t{self.out_number}\t{self.out_number/irradiation_number}')
            if MODE != 0:
                file.write(f'\t{self.out_number/irradiated_atom_num}')
            file.write('\n')
            for name,num in self.out_atom_num.items():
                file.write(f'{name}\t{num}\t{num/irradiation_number}\n')
            file.close()       
    def write_detailed_displaced_atom(self,tag:bool,file_path:str,MODE:int) -> None:
        if tag and MODE!=0:
            file = open(file_path,'w')
            file.write('DISPLACED ATOMS\n')
            file.write(f'TOTAL_NUMBER:{len(self.displaced_atoms)}\n')
            for name,detail in self.displaced_atom_data.items():
                file.write(f'{name}:\n')
                for order in detail:
                    file.write(f'{order}\t')
                file.write('\n')
            file.close()
    def write_detailed_vacancy(self,tag:bool,file_path:str,MODE:int) -> None:
        if tag and MODE!=0:
            file = open(file_path,'w')
            file.write('VACANCIES\n')
            file.write(f'TOTAL NUMBER:{len(self.vacancies)}\n')
            for vac_type,vac_detial in self.vac_data.items():
                file.write(f'####{vac_type}####\n')
                for name,orders in vac_detial.items():
                    file.write(f'{name}:\n')
                    for order in orders:
                        file.write(f'{order}\t')
                    file.write('\n')
            file.close()
    def write_detailed_sia(self,tag:bool,file_path:str,MODE:int) -> None:
        if tag and MODE!=0:
            file = open(file_path,'w')
            file.write('SIAS\n')
            file.write(f'TOTAL NUMBER:{len(self.sias)}\n')
            for sia_type,sia_detial in self.sia_data.items():
                file.write(f'####{sia_type}####\n')
                for name,orders in sia_detial.items():
                    file.write(f'{name}:\n')
                    for order in orders:
                        file.write(f'{order}\t')
                    file.write('\n')
            file.close()
    def write_detailed_subs_atom(self,tag:bool,file_path:str,MODE:int) -> None:
        if tag and MODE!=0:
            file = open(file_path,'w')
            file.write('SUBSTITUTION ATOMS\n')
            file.write(f'TOTAL_NUMBER:{len(self.subs_atoms)}\n')
            for name,detail in self.subs_atom_data.items():
                file.write(f'{name}:\n')
                for order in detail:
                    file.write(f'{order}\t')
                file.write('\n')
            file.close()
    def write_detailed_out_atom(self,tag:bool,file_path:str,MODE:int) -> None:
        if tag and MODE!=0:
            file = open(file_path,'w')
            file.write('FLYING OUT ATOMS\n')
            file.write(f'TOTAL_NUMBER:{len(self.out_atoms)}\n')
            for name,detail in self.out_atom_data.items():
                file.write(f'{name}:\n')
                for order in detail:
                    file.write(f'{order}\t')
                file.write('\n')
            file.close()