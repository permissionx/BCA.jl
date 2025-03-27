from std_io import *
from std_emitter import *
from std_target import *
import numpy as np 
class MySimulation:
    def __init__(self) -> None:
        self.io = InOut()
        self.source = Emitter(self.io.in_emitter)
        self.target = TargetCrystal()
        self.initialization()
        self.irradiation()
        self.output_finally()
    ##########################
    ##system initialization###
    ##########################
    def initialization(self) -> None:
        self.target.read_structure_file(self.io.in_structure,self.io.input_tags['STRUC_TYPE'])
        self.target.read_threshold_energy_file(self.io.in_energy)
        self.target.get_target_size(self.io.input_tags['SIZE_A'],self.io.input_tags['SIZE_B'],self.io.input_tags['SIZE_C'])
        self.target.create_target()
        self.target.write_structure(self.io.input_tags['OUTSTRUC'],self.io.out_structure,0,self.io.input_tags['MODE'])
        #self.target.write_structure(1,self.io.out_structure,0,1)
        self.source.get_emitter_parameters(self.target.sensitive_radius)
        if self.io.input_tags['MODE'] == 0:
            self.source.ion_parameters['time'] = self.io.input_tags['INCITIME']
        self.target.get_conditions(self.source.emitter_parameters,self.source.ion_parameters)
        self.current_num = 0
        self.irradiation_percentage = 0
        self.outputnum = self.source.ion_parameters['time']/100.0
        self.primary_neighbor_width = 1 if int(math.sqrt(self.source.ion_parameters['charge'])/1.4) < 1 else int(math.sqrt(self.source.ion_parameters['charge'])/1.4)
        #self.primary_neighbor_width = 1
    ##########################
    ##irradiation simulation##
    ##########################
    def irradiation(self) -> None:
        while self.current_num < self.source.ion_parameters['time']:
            if self.current_num % 100 == 0:
                print("Irradiation number:", self.current_num)
            #output process
            if self.current_num/self.outputnum >= self.irradiation_percentage:
                self.output_process()
            while self.current_num/self.outputnum >= self.irradiation_percentage:
                self.irradiation_percentage += 1
            #start irradiation
            #step1:produce ion and check its neighbor atoms
            self.current_num += 1
            self.source.create_projectile()
            self.source.incident_ion.atom_energy = 1000
            self.source.incident_ion.atom_velocity = [0.,0.,-1.]
            self.source.incident_ion.atom_position = [15.386099491907107, 19.127032099321863, 5.0]
            self.target.crystal = self.target.crystal[:4]
            #self.target.atom_num = 4
            self.target.crystal[0].atom_position = [16.061705, 18.143275510554176, 0.0]
            self.target.crystal[1].atom_position = [16.76004, 19.35282721125779, 0.0]
            self.target.crystal[2].atom_position = [14.665035, 18.143275510554176, 0.0]
            self.target.crystal[3].atom_position = [13.966700000000001, 19.35282721125779, 0.0]
            self.source.incident_ion.get_land_point(0)
            #find_ion_neighbors = self.target.find_neighbor_atoms_of_given_position(self.primary_neighbor_width,self.source.incident_ion.ion_initial_land_point)           
            find_ion_neighbors = [True,[0,1,2,3]]
            if find_ion_neighbors[0]:
                self.source.incident_ion.get_neighbor_atom_orders(find_ion_neighbors[1])
            else:
                continue
            #step2:primary collision
            ion_index = -1
            exempt_list = []
            cascade_atom =[]
            primary = self.find_partner(self.source.incident_ion,exempt_list)
            #print(f"\ninitial:")
            #print(self.source.incident_ion.atom_energy)
            #print(self.source.incident_ion.atom_velocity)
                #primary = [True,[order_t],[p],[foot]] or [False]
            while primary[0]:
                print([self.target.crystal[id].atom_position for id in primary[1]])
                print(self.source.incident_ion.atom_energy)
                print(self.source.incident_ion.atom_position)
                self.collision(self.source.incident_ion,primary[1],primary[2],primary[3])
                print(self.source.incident_ion.atom_velocity)
                print([self.target.crystal[id].atom_velocity for id in primary[1]])
                print([self.target.crystal[id].atom_energy for id in primary[1]])
                exit()
                #print(f"primary:")
                #print(f"ion:{self.source.incident_ion.atom_energy}\t{self.source.incident_ion.atom_velocity}")
                #print(f"atom{self.target.crystal[primary[1][0]].atom_layer}:{self.target.crystal[primary[1][0]].atom_energy}\t{self.target.crystal[primary[1][0]].atom_velocity}")

                exempt_list += primary[1]
                for atom_order in primary[1]:
                    if self.target.whether_can_atom_cascade(atom_order):
                        cascade_atom.append(atom_order)
                        self.target.displaced_atoms.append(atom_order)
                whether_ion_in_crystal = self.target.check_whether_position_in_crystal(self.source.incident_ion.atom_position)
                if whether_ion_in_crystal[0]:
                    if self.source.incident_ion.atom_energy <= self.io.input_tags['E_CUT']:
                        lattice_now = whether_ion_in_crystal[1]#capture new atom
                        self.target.get_new_atom_from_ion(self.source.incident_ion,lattice_now)
                        ion_index = len(self.target.crystal)-1
                        break
                    else:
                        find_ion_neighbors = self.target.find_neighbor_atoms_of_given_position(self.primary_neighbor_width,self.source.incident_ion.atom_position)
                        self.source.incident_ion.get_neighbor_atom_orders(find_ion_neighbors[1])
                else:
                    break
                primary = self.find_partner(self.source.incident_ion,exempt_list)
            if self.io.input_tags['MODE'] == 0:
                self.target.get_primary_vacancy_data(cascade_atom)
                if len(cascade_atom)==0:
                    self.target.get_defect_data(0)
                    self.target.refresh_system()
                    continue
            if len(cascade_atom)==0:              
                continue
            #step3.cascade process
            i = 0
            while i<len(cascade_atom)-1:
                j = i+1
                while j<len(cascade_atom):
                    if self.target.crystal[cascade_atom[i]].atom_energy < self.target.crystal[cascade_atom[j]].atom_energy:
                        tmp = cascade_atom[i]
                        cascade_atom[i] = cascade_atom[j]
                        cascade_atom[j] = tmp
                    j+=1
                i+=1
            i = 0
            while i < len(cascade_atom):
                exempt_list = [order for order in cascade_atom]
                if ion_index != -1:
                    exempt_list.append(ion_index)
                width = 1 if self.target.crystal[cascade_atom[i]].atom_charge <= 10 else 2
                self.target.crystal[cascade_atom[i]].get_neighbor_atom_orders(
                    self.target.find_neighbor_atoms_of_given_position(
                        width,self.target.crystal[cascade_atom[i]].atom_position)[1])
                cascade = self.find_partner(self.target.crystal[cascade_atom[i]],exempt_list)
                while cascade[0]:
                    self.collision(self.target.crystal[cascade_atom[i]],cascade[1],cascade[2],cascade[3])
                    exempt_list += cascade[1]
                    for partner_order in cascade[1]:
                        if self.target.whether_can_atom_cascade(partner_order):
                            cascade_atom.append(partner_order)
                            self.target.displaced_atoms.append(partner_order)
                    if self.target.whether_can_atom_cascade(cascade_atom[i]):
                        if self.target.crystal[cascade_atom[i]].atom_energy <= self.io.input_tags['E_CUT']:
                            break
                        else:
                            self.target.crystal[cascade_atom[i]].get_neighbor_atom_orders(
                                self.target.find_neighbor_atoms_of_given_position(
                                    width,self.target.crystal[cascade_atom[i]].atom_position)[1])
                            cascade = self.find_partner(self.target.crystal[cascade_atom[i]],exempt_list)
                    else:
                        break
                if cascade[0]==False:
                    self.target.atom_fly_out(cascade_atom[i])
                i += 1
            if self.io.input_tags['MODE'] == 0:
                self.target.displaced_atoms_total = cascade_atom
            #STEP4: defect recombination
            self.target.defect_recombination(self.io.input_tags['CAP_R'])
            #STEP5: get defect data
            if self.io.input_tags['MODE'] == 0:
                self.target.get_defect_data(0)
                self.target.refresh_system()
        if self.io.input_tags['MODE'] != 0:
            self.target.get_defect_data(1)
    ##########################
    ########BCA model#########
    ##########################
    def find_partner(self,atom_i:Atom, exempt_list:list) -> list:
        useful_atom_orders = []
        impact_parameters = []
        foots = []
        names = []
        distance_between_x0_and_foot = []
        for atom_order in atom_i.neighbor_atom_orders:
            if (atom_order in exempt_list) or (atom_order in self.target.out_atoms):
                continue
            information = dis_bet_line_point(atom_i.atom_position,
                                             atom_i.atom_velocity,
                                             self.target.crystal[atom_order].atom_position)
            '''test_0 = impact_parameter_extimate(atom_i.atom_charge,
                                               self.target.crystal[atom_order].atom_charge,
                                                atom_i.atom_mass,
                                                self.target.crystal[atom_order].atom_mass,
                                                atom_i.atom_energy,
                                                self.target.crystal[atom_order].Td)
            if test_0 ==0:
                radius_cut = 0
            else:
                radius_cut = min(test_0,self.target.crystal[atom_order].atom_radius+self.source.incident_ion.atom_radius)'''
            #radius_cut = self.target.crystal[atom_order].atom_radius+self.source.incident_ion.atom_radius
            radius_cut = self.target.cut_radius
            if (information[0] > 0 and 
                information[1] < radius_cut and
                information[2] > atom_i.atom_radius/2):
                name = f'{self.target.crystal[atom_order].atom_name}-{self.target.crystal[atom_order].atom_layer}'
                names.append(name)
                useful_atom_orders.append(atom_order)
                impact_parameters.append(information[1])
                distance_between_x0_and_foot.append(information[2])
                foots.append(information[3])
        if len(useful_atom_orders) == 0:
            return [False]
        elif len(useful_atom_orders) == 1:
            return [True,useful_atom_orders,impact_parameters,foots,distance_between_x0_and_foot]
        else:
            min_d = min(distance_between_x0_and_foot)
            index = distance_between_x0_and_foot.index(min_d)
            result_orders = []
            result_ps = []
            result_foots = []
            for i in range(len(impact_parameters)):
                if (distance_between_x0_and_foot[i] - min_d)<self.target.crystal[useful_atom_orders[index]].atom_radius/2:
                    result_orders.append(useful_atom_orders[i])
                    result_ps.append(impact_parameters[i])
                    result_foots.append(foots[i])
            return [True,result_orders,result_ps,result_foots]
    def collision(self,atom_i:Atom,atoms_t:list,impact_ps:list,foots:list) -> None:
        '''
        atoms_t: list of atom order
        '''
        num_t = len(atoms_t)
        Ts = []
        Qs = []
        new_vectors_for_i = []
        new_vectors_for_t = []
        x1s = []
        x2s = []
        #1.calculate all informations
        i=0
        while i < num_t:
            CON = self.target.crystal[atoms_t[i]].atom_mass/atom_i.atom_mass
            Em = 0.5*atom_i.atom_charge**(4/3)*atom_i.atom_mass*CONSTANT_MU*CONSTANT_VB**2/CONSTANT_EV
            E0 = atom_i.atom_energy
            Er = atom_i.atom_energy*CON/(1+CON)
            Tmax = 4*Er/(1+CON)
            impact_p = impact_ps[i]
            R = R_min(atom_i.atom_charge,
                      self.target.crystal[atoms_t[i]].atom_charge,
                      Er,
                      impact_p)
            theta = get_theta(atom_i.atom_charge,
                              self.target.crystal[atoms_t[i]].atom_charge,
                              Er,
                              impact_p,
                              R,
                              self.io.input_tags['INTEG_START_ACCUR'],
                              self.io.input_tags['INTEG_FAC'])
            tao = get_tao(atom_i.atom_charge,
                          self.target.crystal[atoms_t[i]].atom_charge,
                          Er,
                          impact_p,
                          R,
                          self.io.input_tags['INTEG_START_ACCUR'],
                          self.io.input_tags['INTEG_FAC'])
            Q = electron_energy_loss(atom_i.atom_charge,
                                     self.target.crystal[atoms_t[i]].atom_charge,
                                     atom_i.atom_mass,
                                     self.target.crystal[atoms_t[i]].atom_mass,
                                     atom_i.atom_energy,
                                     impact_p,
                                     self.target.lattice_parameter,
                                     self.target.N,
                                     atom_i.atom_alpha,
                                     atom_i.atom_beta)
            if E0 < Em:
                T_e = Tmax*math.sin(theta/2)**2
                psi = math.atan(CON*math.sin(theta)/(1+CON*math.cos(theta)))
                phi = math.atan(math.sin(theta)/(1-math.cos(theta)))
            else:
                f = math.sqrt(1-Q/Er)
                psi = math.atan(CON*f*math.sin(theta)/(1+CON*f*math.cos(theta)))   #for atom_i
                phi = math.atan(f*math.sin(theta)/(1-f*math.cos(theta)))           #for atom_t
                T_e = (f*math.sin(theta/2)**2+(1-f)**2/4)*Tmax    #energy transferred to atom_t
            initial_velocity = atom_i.atom_velocity
            vector_r = [x-y for x,y in zip(self.target.crystal[atoms_t[i]].atom_position,atom_i.atom_position)]#position r
            spiale = vector_cross(vector_r,initial_velocity)
            x1 = (2*tao+(CON-1)*impact_p*math.tan(theta/2))/(1+CON)
            x2 = impact_p*math.tan(theta/2)-x1
            new_vector_for_i = get_direction_of_vector(vector_rotation(initial_velocity,spiale,psi))
            new_vector_for_t = get_direction_of_vector(vector_rotation(initial_velocity,spiale,-phi))
            Qs.append(Q)
            Ts.append(T_e)
            x1s.append(x1)
            #print(x1)
            x2s.append(x2)
            #print(x2)
            new_vectors_for_i.append(new_vector_for_i)
            new_vectors_for_t.append(new_vector_for_t)
            i+=1
            #print(f"Er: {Er}\ntheta: {theta}\ntau: {tao}\n tanphi: {math.tan(phi)}\ntan_psi: {math.tan(psi)}\nte: {T_e}\nx1: {x1}\nx2: {x2}\nQ: {Q}")  #debug 
            #print("\n")
        #2.change atoms states
        Ttot = sum(Ts)
        beta = num_t*atom_i.atom_energy/(num_t*atom_i.atom_energy+(num_t-1)*Ttot)
        new_Ts = [energy*beta for energy in Ts]
        Qtot_final = sum(Qs) * beta
        Ttot_final = sum(new_Ts)
        initial_velocity = atom_i.atom_velocity
        initial_energy = atom_i.atom_energy
          #for atom_i
        T_final = atom_i.atom_energy - Qtot_final - Ttot_final 
        if T_final <=0:
            T_final = 0
        atom_i.get_new_energy(T_final)#new energy for atom_i

        new_vector_for_i = [0,0,0]
        if  T_final !=0:
            new_vector_for_i[0] += math.sqrt(2*atom_i.atom_mass * initial_energy)*initial_velocity[0]
            new_vector_for_i[1] += math.sqrt(2*atom_i.atom_mass * initial_energy)*initial_velocity[1]
            new_vector_for_i[2] += math.sqrt(2*atom_i.atom_mass * initial_energy)*initial_velocity[2]
            if num_t==1:
                new_vector_for_i = new_vectors_for_i[0]
            else:
                i = 0
                while i < num_t:
                    new_vector_for_i[0] += -new_vectors_for_t[i][0]*math.sqrt(2*new_Ts[i]*self.target.crystal[atoms_t[i]].atom_mass)
                    new_vector_for_i[1] += -new_vectors_for_t[i][1]*math.sqrt(2*new_Ts[i]*self.target.crystal[atoms_t[i]].atom_mass)
                    new_vector_for_i[2] += -new_vectors_for_t[i][2]*math.sqrt(2*new_Ts[i]*self.target.crystal[atoms_t[i]].atom_mass)
                    i+=1
        atom_i.get_new_velocity(new_vector_for_i)#new velocity direction for atom_i

        new_x1 = sum([x*beta for x in x1s])
        new_foot = [0,0,0]
        for vector in foots:
            new_foot[0] += vector[0]
            new_foot[1] += vector[1]
            new_foot[2] += vector[2]
        new_foot = [x/num_t for x in new_foot]
        new_position_for_i = [x-new_x1*y for x,y in zip(new_foot,initial_velocity)] 
        atom_i.get_new_position(new_position_for_i)#new position for atom_i

          #for atom_t
        i = 0
        while i < num_t:
            self.target.crystal[atoms_t[i]].get_new_energy(new_Ts[i])
            self.target.crystal[atoms_t[i]].get_new_velocity(new_vectors_for_t[i])
            initial_position = self.target.crystal[atoms_t[i]].atom_position
            new_position_for_t = [x+ beta*x2s[i]*y for x,y in zip(initial_position,initial_velocity)]
            self.target.crystal[atoms_t[i]].get_new_position(new_position_for_t)
            i+=1
        
        initial_p = np.sqrt(2 * atom_i.atom_mass * initial_energy) * np.array(initial_velocity)
        after_p = np.sqrt(2 * atom_i.atom_mass * atom_i.atom_energy) * np.array(atom_i.atom_velocity)
        while i < num_t:
            target = self.target.crystal[atoms_t[i]]
            after_p += np.sqrt(2 * target.atom_mass * target.atom_energy) * np.array(target.atom_velocity)
        with open("p_debug.csv", 'a') as file:
            file.write(f"{num_t},{initial_p[0]},{initial_p[1]},{initial_p[2]},{after_p[0]},{after_p[1]},{after_p[2]}\n")
    ##########################
    ####output information####
    ##########################
    def output_process(self) -> None:
        file = open(self.io.out_process,'a+')
        if self.current_num==0:
            file.write(f'TOTAL\tCURRENT\tPERCENTAGE\n')
        file.write(f"{self.source.ion_parameters['time']}\t{self.current_num}\t{self.irradiation_percentage}%\n")
        #print(f"{self.source.ion_parameters['time']}\t{self.current_num}\t{self.irradiation_percentage}%")
        file.close()
        self.target.write_structure(self.io.input_tags['OUTSTRUC'],self.io.out_structure,self.irradiation_percentage,self.io.input_tags['MODE'])
    def output_finally(self) -> None:
        self.output_process()
        self.target.write_structure(self.io.input_tags['OUTSTRUC'],self.io.out_structure,100,self.io.input_tags['MODE'])
        self.target.write_summary_defect(self.io.input_tags['OUTDEFECT'],self.io.out_defect,self.source.ion_parameters['time'],self.io.input_tags['MODE'])
        self.target.write_detailed_displaced_atom(self.io.input_tags['OUTDISATOM'],self.io.out_displaced_atom,self.io.input_tags['MODE'])
        self.target.write_detailed_vacancy(self.io.input_tags['OUTVACANCY'],self.io.out_vacancy,self.io.input_tags['MODE'])
        self.target.write_detailed_sia(self.io.input_tags['OUTSIA'],self.io.out_sia,self.io.input_tags['MODE'])
        self.target.write_detailed_subs_atom(self.io.input_tags['OUTSUBSATOM'],self.io.out_subs_atom,self.io.input_tags['MODE'])
        self.target.write_detailed_out_atom(self.io.input_tags['OUTFLYATOM'],self.io.out_out_atom,self.io.input_tags['MODE'])
        