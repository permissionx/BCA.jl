from std_atom import *
import random
class Emitter:
    def __init__(self,file_path) -> None:
        self.read_parameters(file_path)
        self.get_ion_parameters()

    def read_parameters(self,file_path) -> None:
        file = open(file_path,'r',encoding= 'gbk')
        contents = file.readlines()
        file.close()
        parameters = ['SHAPE','NAME','CHARGE','MASS','ENERGY','FLUX','SANGLE','PANGLE']
        self.parameters = dict()
        for content in contents:
            if content[0]!='#':
                content1 = [x.strip() for x in content.strip().split('=')]
                if content1[0] in parameters:
                    parameters.remove(content1[0])
                else:
                    raise NameError(f'The name of tag {content1[0]} is not correct, or this tag is repetitive.')
                content2 = content1[1].split()[0]
                self.parameters[content1[0]]= content2

    def get_ion_parameters(self)-> None:
        self.ion_parameters = dict()
        self.ion_parameters['name'] = self.parameters['NAME']
        self.ion_parameters['charge']=int(self.parameters['CHARGE'])
        self.ion_parameters['mass']=int(self.parameters['MASS'])
        index = ALL_ELEMENTS_NAME.index(self.parameters['NAME'])
        self.ion_parameters['radius'] = ALL_ELEMENTS_RADIUS[index]
        self.ion_parameters['energy'] = float(self.parameters['ENERGY'])
        self.ion_parameters['alpha'] = ALL_ELEMENTS_ALPHA[index]
        self.ion_parameters['beta'] = ALL_ELEMENTS_BETA[index]
        self.ion_parameters['sangle'] = float(self.parameters['SANGLE'])*math.pi/180
        self.ion_parameters['pangle'] = float(self.parameters['PANGLE'])*math.pi/180
        self.ion_parameters['R'] = multi_matrices([rotation_matrix('y',-self.ion_parameters['sangle']),rotation_matrix('z',self.ion_parameters['pangle'])])
        self.ion_parameters['direction'] = vector_multi_matrix([0,0,-1],self.ion_parameters['R'])
        
    def get_emitter_parameters(self,R0) -> None:
        self.emitter_parameters = dict()
        self.emitter_parameters['height'] = 1000    #A
        if self.parameters['SHAPE']=='C' or self.parameters['SHAPE']=='.C.':
            self.emitter_parameters['shape'] = 'C'
            self.emitter_parameters['radius'] = R0*math.cos(self.ion_parameters['sangle'])
            self.emitter_parameters['area'] = self.emitter_parameters['radius']**2*math.pi      #A2
            self.emitter_parameters['irradiation_area'] = self.emitter_parameters['area']/math.cos(self.ion_parameters['sangle'])
        elif self.parameters['SHAPE']=='S' or self.parameters['SHAPE']=='.S.':
            self.emitter_parameters['shape'] = 'S'
            self.emitter_parameters['length'] = R0*math.cos(self.ion_parameters['sangle'])
            self.emitter_parameters['area'] = self.emitter_parameters['length']**2
            self.emitter_parameters['irradiation_area'] = self.emitter_parameters['area']/math.cos(self.ion_parameters['sangle'])
        elif self.parameters['SHAPE']=='T' or self.parameters['SHAPE']=='.T.':
            self.emitter_parameters['shape'] = 'T'
            self.emitter_parameters['length'] = R0*math.cos(self.ion_parameters['sangle'])
            self.emitter_parameters['area'] = self.emitter_parameters['length']**2*math.sqrt(3)/4
            self.emitter_parameters['irradiation_area'] = self.emitter_parameters['area']/math.cos(self.ion_parameters['sangle'])
        else:
            raise ValueError(f'The shape of ion emitter is out of this code.')
        self.ion_parameters['time'] = int(float(self.parameters['FLUX'])*self.emitter_parameters['area']*1e-4)
    
    def random_position(self) -> list:
        if self.emitter_parameters['shape'] == 'C':
            x = (2*random.random()-1) * self.emitter_parameters['radius']
            y = (2*random.random()-1) * self.emitter_parameters['radius']
            if x**2+y**2 <= self.emitter_parameters['radius']**2:
                return vector_multi_matrix([x,y,self.emitter_parameters['height']],self.ion_parameters['R'])
            else:
                return self.random_position()
        elif self.emitter_parameters['shape'] == 'S':
            x = (random.random()-0.5) * self.emitter_parameters['length']
            y = (random.random()-0.5) * self.emitter_parameters['length']
            return vector_multi_matrix([x,y,self.emitter_parameters['height']],self.ion_parameters['R'])
        else:
            x = (random.random()-0.5) * self.emitter_parameters['length']
            y = random.random()* self.emitter_parameters['length']
            if y < math.sqrt(3)*(abs(x)+0.5*self.emitter_parameters['length']):
                return vector_multi_matrix([x,y,self.emitter_parameters['height']],self.ion_parameters['R'])
            else:
                return self.random_position()

    def gaussian_distribution(self) -> float:
        Energy=random.gauss(self.ion_parameters['energy'],0.1*self.ion_parameters['energy'])
        if 0.7*self.ion_parameters['energy']<Energy<1.3*self.ion_parameters['energy']:
            return Energy
        else:
            return self.gaussian_distribution()
    
    def create_projectile(self) -> None:
        self.incident_ion = ProjectileIon(self.ion_parameters['name'],
                                          self.ion_parameters['charge'],
                                          self.ion_parameters['mass'],
                                          self.ion_parameters['radius'],
                                          self.random_position(),
                                          self.ion_parameters['energy'],
                                          self.ion_parameters['direction'],
                                          self.ion_parameters['alpha'],
                                          self.ion_parameters['beta'])
        #self.incident_ion.atom_position = [0,1.36,100]
        #self.incident_ion.atom_energy = 100000000