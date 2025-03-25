import os
import datetime
class InOut:
    def __init__(self) -> None:
        self.file_path = os.getcwd()
        self.get_input_file_path()
        self.get_output_file_path()

    def get_input_file_path(self) -> None:
        '''
           Get INPUT files pathes.
        '''
        self.in_path = self.file_path + '/INPUT'
        self.in_structure = self.in_path + f'/InStructure.dat'   
        self.in_tag = self.in_path + '/InTag.dat'
        self.in_energy = self.in_path + '/InEnergy.dat'                           
        self.in_emitter = self.in_path +'/InEmitter.dat' 
        self.read_tag_file()

    def get_output_file_path(self) -> None:
        '''
           Get OUTPUT files pathes.
        '''
        self.out_path = self.file_path +'/OUTPUT'
        time_str = str(datetime.datetime.now())
        time_now = time_str[:19].replace(':','-')
        time_now = time_now.replace(' ','@')
        self.out_path += f'/{time_now}'
        isExistOutput = os.path.exists(self.out_path)
        if not isExistOutput:
            os.makedirs(self.out_path)
        self.out_structure = self.out_path + '/OutStructure.dat'
        self.out_process = self.out_path + '/StdOut'
        self.out_displaced_atom = self.out_path + '/OutDisplacedAtoms.dat'
        self.out_defect = self.out_path + '/OutDefects.dat'
        self.out_vacancy = self.out_path + '/OutVacancies.dat'
        self.out_sia = self.out_path + '/OutSias.dat'
        self.out_subs_atom = self.out_path + '/OutSubsAtoms.dat'
        self.out_out_atom = self.out_path + '/OutFlyingAtoms.dat'
            
    def read_tag_file(self) -> None:
        '''
           Read output tags
        '''
        tags_1 = ['MODE','INCITIME']
        tags_2 = ['OUTSTRUC','OUTDEFECT','OUTDISATOM','OUTVACANCY','OUTSIA','OUTSUBSATOM','OUTFLYATOM']
        tags_3 = ['STRUC_TYPE','SIZE_A','SIZE_B','SIZE_C']
        tags = ['OUTSTRUC','OUTDEFECT','OUTDISATOM','OUTVACANCY','OUTSIA','OUTSUBSATOM','OUTFLYATOM',
                'STRUC_TYPE','SIZE_A','SIZE_B','SIZE_C',
                'CAP_R','E_CUT','INTEG_START_ACCUR','INTEG_FAC',
                'MODE','INCITIME']
        yes = ['T','t','.T.','.t.','True','true','TRUE']
        file = open(self.in_tag,'r')
        contents = file.readlines()
        file.close()
        self.input_tags = dict()
        for content in contents:
            if content[0]!='#':
                content1 = [x.strip() for x in content.strip().split('=')]
                if content1[0] in tags:
                    content2 = content1[1].split()[0]
                    if content1[0] in tags_1 or content1[0] in tags_3:
                        self.input_tags[content1[0]] = int(content2)
                    elif content1[0] in tags_2:
                        self.input_tags[content1[0]] = True if content2 in yes else False
                    else:
                        self.input_tags[content1[0]] = float(content2)
                    tags.remove(content1[0])
                else:
                    raise NameError(f'The name of tag {content1[0]} is not correct, or this tag is repetitive.')
        if self.input_tags['MODE']==1 and 'INCITIME' in tags:
            tags.remove('INCITIME')
        self.input_tags['INTEG_START_ACCUR'] = self.input_tags['INTEG_START_ACCUR']*1e-15
        if len(tags)!=0:
            raise ValueError(f'Missing settings of {tags}')