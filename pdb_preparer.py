#!/nfs/u1/local/anaconda3/bin/python
from optparse import OptionParser
import re
import numpy as np
import pandas as pd
import copy
from itertools import groupby
import operator

class optParser():
    def __init__(self,fakeArgs):
        usage = 'This script is used to generate a protein pdb file or a protein pdb file and small molecule pdb file for automd.sh or gengmxfep.sh, which requires the input pdb file processed by Protein Preparation Wizard module of Schrodinger Maestro.\n \nPresented by phloroglucinol.'
        parser = OptionParser(usage)
        parser.add_option('-f', '--file', dest = 'file', help = 'The name of your pdb file that will be recognized by this script. Default: protein.pdb', default = 'protein.pdb')
        parser.add_option('-l', '--iflig', dest = 'iflig', help = 'Whether a small molecule ligand exists in the pdb file. Default: 1', default = 1)

        if fakeArgs:
            self.option, self.args = parser.parse_args(fakeArgs)
        else:
            self.option, self.args = parser.parse_args()


class ATOM():
    def __init__(self, tpl):
        self.record_name = tpl[0]
        self.atom_num = tpl[1]
        self.atom_name = tpl[2]
        self.altLoc = tpl[3]
        self.res_name = tpl[4]
        self.chainID = tpl[5]
        self.res_seq = tpl[6]
        self.iCode = tpl[7]
        self.x_coor = tpl[8]
        self.y_coor = tpl[9]
        self.z_coor = tpl[10]
        self.occupancy = tpl[11]
        self.tempFactor = tpl[12]
        self.element = tpl[13]
        self.charge = tpl[14]
        if self.charge == 1:
            self.chg_str = '1+'
        elif self.charge == -1:
            self.chg_str ='1-'
        else:
            self.chg_str = ''
        xyz=[self.x_coor, self.y_coor, self.z_coor]
        self.xyz = np.array(xyz)
        self.terminal = False

    def __str__(self):
        return "Atom: {} in {}_{} with atom serial number {} and element {}".format(self.atom_name, self.res_name, self.res_seq, self.atom_num, self.element)
    
    def __repr__(self):
        return self.__str__()
    
    def write_atm_line(self):
        l_str='{: <6s}{: >5d} {: ^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{: >6.2f}{: >6.2f}          {: >2s}{:2s}'.format(
            self.record_name,
            self.atom_num,
            self.atom_name,
            self.altLoc,
            self.res_name,
            self.chainID,
            self.res_seq,
            self.iCode,
            self.x_coor,
            self.y_coor,
            self.z_coor,
            self.occupancy,
            self.tempFactor,
            self.element,
            self.chg_str,
            )
        if not self.terminal:
            l_str += '\n'
        else:
            l_str += '\nTER\n'
        return l_str
    
class RESIDUE():
    def __init__(self, seq, lst):
        self.seq = seq
        self.atm_in_res = lst
        self.name = lst[0].res_name
        ori_atm_in_res = self.atm_in_res
        if self.name == 'MOL':
            self.convert_MOL_atom()
        elif self.name == 'PTR':
            self.name = 'Y2P'
            self.convert_PTR2TYR()
        elif self.name == 'ZNA':
            self.convert_RECORD()
            self.atm_in_res[0].atom_name = 'Zn' # could be changed to Zn according to the prepi file you use
            self.atm_in_res[-1].terminal = True
        elif self.name == 'MG':
            self.convert_RECORD()
            self.atm_in_res[-1].terminal = True
        elif self.name == 'TPO':
            self.convert_RECORD()
        elif self.name == 'SO4':
            self.convert_RECORD()
            self.atm_in_res[-1].terminal = True
        elif self.name == 'HOH':
            self.atm_in_res[-1].terminal = True
        elif self.name == 'CYX':
            for atm in ori_atm_in_res:
                if atm.atom_name == 'HG':
                    self.atm_in_res.remove(atm)
        elif self.name == 'ASP':
            for atm in ori_atm_in_res:
                if atm.atom_name in ['HD2','HA2' ]:
                    self.atm_in_res.remove(atm)
                elif 'HX' in atm.atom_name:
                    self.atm_in_res.remove(atm)
        elif self.name == 'GLU':
            for atm in ori_atm_in_res:
                if atm.atom_name in ['HE2','HA2']:
                    self.atm_in_res.remove(atm)
        elif self.name in ['HIE']:
            for atm in ori_atm_in_res:
                if atm.atom_name == 'HD1':
                    self.atm_in_res.remove(atm)

        self.res_df = pd.DataFrame()
        for atm in ori_atm_in_res:
            # remove the conformation ######################################
            if atm.altLoc == '' or atm.altLoc == 'A':
                pass
            else:
                self.atm_in_res.remove(atm)
            ################################################################
            # atom_name_lst.append(atm.atom_name)
            if atm.atom_name == 'C':
                self.C_atom_xyz = atm.xyz
            elif atm.atom_name == 'N':
                self.N_atom_xyz = atm.xyz
            

    @property
    def net_chg(self):
        net_charge = 0
        for i in self.atm_in_res:
            net_charge += i.charge          
        return net_charge
    
    def __str__(self):
        return "Residue: {}_{}, and the atm_in_res is {}".format(self.name, self.seq, self.atm_in_res)
    
    def __repr__(self):
        return self.__str__()
    
    def convert_PTR2TYR(self):
        new_atm_in_res = copy.deepcopy(self.atm_in_res)
        for i in new_atm_in_res:
            if i.res_name == 'PTR':
                i.record_name = 'ATOM'
                i.res_name = 'Y2P'
                if i.atom_name == 'OH':
                    i.atom_name = 'OG'
        self.atm_in_res = new_atm_in_res
    
    def convert_MOL_atom(self):
        _atm_in_res = copy.deepcopy(self.atm_in_res)
        sorted_atm_in_res = sorted(_atm_in_res, key=operator.attrgetter('element'))
        grouped_atoms = [list(grp_result) for key, grp_result
                         in groupby(sorted_atm_in_res, key=operator.attrgetter('element'))]

        new_atm_in_res = []
        for every_group_atoms in grouped_atoms:
            for idx, atm in enumerate(every_group_atoms):
                atm.atom_name = f'{atm.element}{idx+1}'
                new_atm_in_res.append(atm)
        self.atm_in_res = _atm_in_res

    def convert_RECORD(self):
        new_atm_in_res = copy.deepcopy(self.atm_in_res)
        for i in new_atm_in_res:
            i.record_name = 'ATOM'
        self.atm_in_res = new_atm_in_res

    def generate_res_df(self):
        column_lst = ['atom_num', 'atom_name', 'res_num', 'res_name',
                      'chainID', 'xyz', 'occupancy', 'tempFactor', 'charge']
        index_lst = [i.atom_name for i in self.atm_in_res]
        self.res_df = pd.concat([pd.DataFrame({'atom_num': [i.atom_num],
                                               'atom_name': [i.atom_name],
                                               'res_num': [i.res_seq],
                                               'res_name': [i.res_name],
                                               'chainID': [i.chainID],
                                               'xyz': [i.xyz],
                                               'occupancy': [i.occupancy],
                                               'tempFactor': [i.tempFactor],
                                               'charge': [i.charge]
                                               },columns=column_lst) for i in self.atm_in_res], ignore_index=True)
        self.res_df.index = index_lst
        return self.res_df
    
    def write_res_line(self):
        l_str = ''
        for i in self.atm_in_res:
            l_str += i.write_atm_line()
        return l_str
    
    def rm_hydrogen(self):
        h_pattern = re.compile(r'^H.*')
        new_atm_in_res = self.atm_in_res.copy()
        for i in new_atm_in_res:
            if h_pattern.match(i.atom_name):
                self.atm_in_res.remove(i)
    
    def change_name(self, new_name):
        self.name = new_name
        for i in self.atm_in_res:
            i.res_name = new_name


class PDB_PREPARER():
    '''
    This script need a pdb files, 
    which is handled by Schrodinger-Maestro software
    (including adding hydrogens, 
    filling missing atoms, 
    optimize hydrogen bond assignment).
    '''
    def __init__(self, file_name, iflig):
        self.file_name = file_name
        self.iflig = iflig       
        self.atom_lst = []
        self.res_seq_lst = []
        self.residue_lst = []
        self.cys_res_lst = []
        self.cyx_pairs = []
        self.stretch_dct = {}
        self.read_pdb(self.file_name)
        self.common_residue_name = ["ALA","ARG","ASH","ASN","ASP","CYM","CYS",
                                    "CYX","GLH","GLN","GLU","GLY","HID","HIE",
                                    "HIP","HIS","ILE","LEU","LYN","LYS","MET",
                                    "PHE","PRO","SER","THR","TRP","TYR","VAL",
                                    "PTR","Y2P","HSE","HSD","HSP"]
        
    def read_pdb(self, pdbin):
        
        file_in = open(pdbin)
        lines = file_in.readlines()
        self.parse_pdb_lines(lines)
    
    def parse_pdb_lines(self, pdb_lines):
        
        for line in pdb_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                
                record_name = line[:6].strip()
                atom_num = int(line[6:11])
                atom_name = line[11:16].strip()
                altLoc = line[16].strip() # to identify different conformation
                
                if line.startswith('HETATM'):
                    if line[17:20].strip() == 'HOH':
                        res_name = 'HOH'
                    elif line[17:20].strip() == 'PTR':
                        res_name = 'PTR'
                    elif line[17:20].strip() == 'ZN':
                        res_name = 'ZNA'
                    elif line[17:20].strip() == 'MG':
                        res_name = 'MG'
                    elif line[17:20].strip() == 'SO4':
                        res_name = 'SO4'
                    elif line[17:20].strip() == 'TPO':
                        res_name = 'TPO'
                    elif line[17:20].strip() == 'CA':
                        res_name = 'CA'
                    else:
                        res_name = 'MOL'
                else:
                    res_name = line[17:20].strip()
                
                chainID = line[21].strip()
                res_seq = int(line[22:26])
                iCode = line[26]
                x_coor = float(line[30:38])
                y_coor = float(line[38:46])
                z_coor = float(line[46:54])
                occupancy = line[54:60].strip()
                tempFactor = line[60:66].strip()
                element = line[76:78].strip()
                charge = line[78:80].strip()
                if charge == "":
                    chg = 0
                elif charge == "1+":
                    chg = 1
                elif charge == "1-":
                    chg = int(-1)
                elif charge == "2+":
                    chg = 2
                elif charge == "2-":
                    chg = int(-2)
                
                if occupancy == "":
                    occupancy = 0.0
                else:
                    occupancy = float(occupancy)

                if tempFactor == "":
                    tempFactor = 0.0
                else:
                    tempFactor = float(tempFactor)

                atom_tuple = (record_name, 
                              atom_num, 
                              atom_name,
                              altLoc,
                              res_name,
                              chainID,
                              res_seq,
                              iCode,
                              x_coor,
                              y_coor,
                              z_coor,
                              occupancy,
                              tempFactor,
                              element,
                              chg,
                             )
                atm_obj = ATOM(atom_tuple)
                self.res_seq_lst.append(res_seq)
                self.atom_lst.append(atm_obj)
        
        #Generate a list(one_res_atms) of atoms belonging to the same residue (according to residue name, sequence and chainID).
        #Then generate the residue_lst consist of all residues(object type) in protein.
        one_res_atms = []
        seq = self.res_seq_lst[0] #will update in the following cycle
        name = self.atom_lst[0].res_name #will update in the following cycle
        chain = self.atom_lst[0].chainID #will update in the following cycle
        for atm in self.atom_lst:         
            if atm.res_seq == seq and atm.res_name == name and atm.chainID == chain:
                one_res_atms.append(atm)
            else:
                oneres = RESIDUE(seq, one_res_atms)
                # print(f'___________{oneres}_____________\n')
                self.residue_lst.append(oneres)
                one_res_atms = [atm] #reinitialize {one_res_atms} with the first atom in the next residue.
            seq = atm.res_seq
            name = atm.res_name
            chain = atm.chainID
        oneres = RESIDUE(seq, one_res_atms) #generate RESIDUE object for the last res.
        self.residue_lst.append(oneres)

        self.cys_res_v_sg_k = {}
        for res in self.residue_lst:
            if res.name == 'CYS' or res.name == 'CYX':
                for atm in res.atm_in_res:
                    if atm.atom_name == 'SG':
                        self.cys_res_v_sg_k[atm] = res
        self.cyx_pairs = []
        processed_cys = set()
        for sg_atm, cys in self.cys_res_v_sg_k.items():
            for sg_atm_, cys_ in self.cys_res_v_sg_k.items():
                if cys in processed_cys or cys_ in processed_cys:
                    continue
                if sg_atm == sg_atm_:
                    continue
                sgs_distance = np.around(np.sqrt(np.sum((sg_atm.xyz - sg_atm_.xyz) ** 2)),2)
                if sgs_distance <= 2.06: # DISTANCE BETWEEN CYX.SG AND CYX.SG IS 2.05 A.
                    self.cyx_pairs.append((cys, cys_))
                    processed_cys.add(cys)
                    processed_cys.add(cys_)
                    for res in self.residue_lst:
                        if res.atm_in_res == cys.atm_in_res:
                            res.change_name('CYX')
                        elif res.atm_in_res == cys_.atm_in_res:
                            res.change_name('CYX')
        self.log_sbond_info('model')


            
    def log_sbond_info(self, rec_obj_name='model'):
        '''
        Loging the sbond information as tleap bond style, the default protein obj name is 'model'
        '''
        with open('sbond.lst_bypdb_preparer', 'w') as f:
            for sbond_pair in self.cyx_pairs:
                print(f'bond {rec_obj_name}.{sbond_pair[0].seq}.SG bond {rec_obj_name}.{sbond_pair[1].seq}.SG', file=f)
            

    def generate_stretch(self):
        '''
        This function will generate self.stretch_dct (a dictionary), each of whose values is 
        one stretch (sequential residues). The keys of self.stretch_dct is int that start from zero.
        
        Example:
        >>> fakeArgs="-f 5WA5_handled.pdb -l 1" #only keep this for test purpose
        >>> opts=optParser(fakeArgs.strip().split()) #only keep this for test purpose
        >>> ts = PDB_PREPARER(opts.option.file, bool(opts.option.iflig))
        >>> ts.generate_stretch()
        >>> ts.stretch_dct
        {0: [Residue: GLU_1,
        .......
          Residue: GLU_326],
         1: [Residue: MOL_330]}
        '''
        stretch_num = 0
        needadd=False
        one_stretch = []
        total_residue_num = len(self.residue_lst)
        for res in self.residue_lst:
            if needadd:
                stretch_num+=1
                one_stretch=[]
                one_stretch.append(res)
                self.stretch_dct[stretch_num] = one_stretch
                next_res_idx = self.residue_lst.index(res)+1
                if next_res_idx < total_residue_num:
                    NC_distan = self.cal_NC_distance(res, self.residue_lst[next_res_idx])
                    # if res.seq+1 == self.residue_lst[next_res_idx].seq:
                    if NC_distan < 1.4:
                        needadd=False
                    else:
                        needadd=True
                else:
                    pass                    
            else:
                one_stretch.append(res)
                self.stretch_dct[stretch_num]=one_stretch
                next_res_idx = self.residue_lst.index(res)+1
                if next_res_idx < total_residue_num:
                    NC_distan = self.cal_NC_distance(res, self.residue_lst[next_res_idx])
                    # if res.seq+1 == self.residue_lst[next_res_idx].seq:
                    if NC_distan < 1.4:
                        needadd=False
                    else:
                        needadd=True
                else:
                    pass
#         return self.stretch_dct

    def cal_NC_distance(self, residue_obj1, residue_obj2):
        '''
        This function will calculate the distance from atom C in the first given residue to atom N 
        in the next residue, which is used to determine whether these two residues are connected.
        The data type of return value is float or int. If two residues are in different strenches, the return value will be 2.
        
        Example:
        >>> fakeArgs="-f 5WA5_handled.pdb -l 1" #only keep this for test purpose
        >>> opts=optParser(fakeArgs.strip().split()) #only keep this for test purpose
        >>> ts = PDB_PREPARER(opts.option.file, bool(opts.option.iflig))
        >>> residue_1 = ts.residue_lst[0]
        >>> residue_2 = ts.residue_lst[1]
        >>> ts.cal_NC_distance(residue_1, residue_2)
        1.33
        '''
        try:
            xyzC = residue_obj1.C_atom_xyz
            xyzN = residue_obj2.N_atom_xyz
            NC_distance = np.sqrt(np.sum((xyzN - xyzC) ** 2))
            # print(f'Atom_C_xyz is {xyzC}, the next residue\'s atom_N_xyz is {xyzN}, distance is {NC_distance}')
        except:
            # print(f'<{residue_obj2.name} {residue_obj2.seq}> dosen\'t have N_atom.')
            NC_distance = 2.0
        return NC_distance

    def rm_stretch_fl_hydrogen(self):
        '''
        This function will remove hydrogen of the first residue and last residue of 
        one stretch (sequential residues) except water and ligands.
        
        Example:
        
        >>> fakeArgs="-f 5WA5_handled.pdb -l 1" #only keep this for test purpose
        >>> opts=optParser(fakeArgs.strip().split()) #only keep this for test purpose
        >>> ts = PDB_PREPARER(opts.option.file, bool(opts.option.iflig))
        >>> ts.generate_stretch()
        >>> ts.rm_stretch_fl_hydrogen()
        >>> ts.stretch_dct[0][0].atm_in_res
        >>> ts.stretch_dct[0][-1].atm_in_res
        '''
        for stretch in self.stretch_dct.values():
            if stretch[0].name != 'MOL' and stretch[0].name != 'HOH':
                stretch[0].rm_hydrogen()
                if stretch[-1].name != 'MOL' and stretch[-1].name != 'HOH':
                    stretch[-1].rm_hydrogen()

    def get_MOL_net_chg(self):
        '''
        This function will return the net_charge of the ligand in the pdb file. 
        By default, the number of ligand in the pdb file is one.
        The data type of return value is int.
        
        Example:
        >>> fakeArgs="-f 5WA5_handled.pdb -l 1" #only keep this for test purpose
        >>> opts=optParser(fakeArgs.strip().split()) #only keep this for test purpose
        >>> ts = PDB_PREPARER(opts.option.file, bool(opts.option.iflig))
        >>> ts.get_MOL_net_chg()
        1
        '''
        for i in self.residue_lst:
            if i.name == 'MOL':
                return i.net_chg

    def correct_HIS_name(self, stretch_dict):
        '''
        This function will correct the HIS name according to the numbers of HD and HE within the HIS,
        two HD and two HE: HIP,
        two HD and one HE: HID,
        one HD and two HE: HIE,
        
        Example:
        >>> fakeArgs="-f 5WA5_handled.pdb -l 1" #only keep this for test purpose
        >>> opts=optParser(fakeArgs.strip().split()) #only keep this for test purpose
        >>> ts = PDB_PREPARER(opts.option.file, bool(opts.option.iflig))
        >>> ts.generate_stretch()
        >>> self.correct_HIS_name(self.stretch_dct)
        '''
        dct = copy.deepcopy(stretch_dict)
        for stretch in dct.values():
            for res in stretch:
                if res.name == 'HIS':
                    res_df = res.generate_res_df()
                    HD_lst = []
                    HD_pattern = re.compile(r'^HD.*')
                    HE_lst = []
                    HE_pattern = re.compile(r'^HE.*')
                    for i in res_df['atom_name']:
                        if HD_pattern.search(i):
                            HD_lst.append(i)
                        elif HE_pattern.search(i):
                            HE_lst.append(i)
                    if len(HD_lst) == 1 and len(HE_lst) == 2:
                        res.change_name('HIE')
                    elif len(HD_lst) == 2 and len(HE_lst) == 1:
                        res.change_name('HID')
                    elif len(HD_lst) == 2 and len(HE_lst) == 2:
                        res.change_name('HIP')
        self.stretch_dct = dct
#         return self.stretch_dct

    
    def prepare_4_amberleap(self):
        '''
        This function will generate the self.stretch_dct, whose the first and the last residues hydrogen is removed 
        and HIS residues names are corrected.
        
        Example:
        >>> fakeArgs="-f 5WA5_handled.pdb -l 1" #only keep this for test purpose
        >>> opts=optParser(fakeArgs.strip().split()) #only keep this for test purpose
        >>> ts = PDB_PREPARER(opts.option.file, bool(opts.option.iflig))
        >>> ts.prepare_4_amberleap()
        '''
        self.generate_stretch()
        self.rm_stretch_fl_hydrogen()
        self.correct_HIS_name(self.stretch_dct)
    
    def write_pdb(self, stretch_dict, protein_pdb_name='protein.pdb', ligand_pdb_name='mol.pdb'):
        '''
        This function need a self.stretch_dct, it will generate a protein.pdb and a mol.pdb (if self.iflig=1)

        Example:
        >>> fakeArgs="-f 5WA5_handled.pdb -l 1" #only keep this for test purpose
        >>> opts=optParser(fakeArgs.strip().split()) #only keep this for test purpose
        >>> ts = PDB_PREPARER(opts.option.file, bool(opts.option.iflig))
        >>> ts.write_pdb(ts.stretch_dct, protein_pdb_name='protein.pdb', ligand_pdb_name='mol.pdb')
        '''
        dct = copy.deepcopy(stretch_dict)
        if self.iflig:
            mol_line = []
            rec_line = []
            for stretch in dct.values():
                for res in stretch:
                    if res.name == 'MOL':
                        mol_line.append(res.write_res_line())
                    else:
                        rec_line.append(res.write_res_line())
                # if stretch[0].name != 'MOL':
                if stretch[0].name in self.common_residue_name:
                    rec_line.append('TER\n')
                
            with open(ligand_pdb_name, 'w', encoding='utf-8') as ligfile:
                for i in mol_line:
                    ligfile.write(i)
            with open(protein_pdb_name, 'w', encoding='utf-8') as recfile:
                for i in rec_line:
                    recfile.write(i)
        else:
            rec_line = []
            for stretch in dct.values():
                for res in stretch:
                    if res.name != 'MOL':
                        rec_line.append(res.write_res_line())
                if stretch[0].name in self.common_residue_name:
                    rec_line.append('TER\n')
            with open(protein_pdb_name, 'w', encoding='utf-8') as recfile:
                for i in rec_line:
                    recfile.write(i)
if __name__ == "__main__":
    # fakeArgs="-f protein_prev_1.pdb -l 1" #only keep this for test purpose
    # opts=optParser(fakeArgs.strip().split()) #only keep this for test purpose
    opts=optParser('')
    # print(int(opts.option.iflig))
    # print(bool(opts.option.iflig))
    ts = PDB_PREPARER(opts.option.file, int(opts.option.iflig))
    ts.prepare_4_amberleap()
    ts.write_pdb(ts.stretch_dct,)
    if int(opts.option.iflig):
        print(ts.get_MOL_net_chg())
    else:
        pass
    
