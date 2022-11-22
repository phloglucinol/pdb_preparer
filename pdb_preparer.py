#!/nfs/u1/local/amber20/miniconda/bin/python
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
        if self.name == 'MOL':
            self.convert_MOL_atom()
        elif self.name == 'PTR':
            self.name = 'Y2P'
            self.convert_PTR2TYR()
        elif self.name == 'ZNA':
            self.convert_RECORD()
            self.atm_in_res[0].atom_name = 'Zn'
            self.atm_in_res[-1].terminal = True
        elif self.name == 'MG':
            self.convert_RECORD()
            self.atm_in_res[-1].terminal = True
        elif self.name == 'SO4':
            self.convert_RECORD()
            self.atm_in_res[-1].terminal = True
        elif self.name == 'HOH':
            self.atm_in_res[-1].terminal = True
        self.res_df = pd.DataFrame()

    @property
    def net_chg(self):
        net_charge = 0
        for i in self.atm_in_res:
            net_charge += i.charge          
        return net_charge
    
    def __str__(self):
        return "Residue: {}_{}".format(self.name, self.seq)
    
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
                altLoc = line[16].strip()
                
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
        
        order_duplicate = []
        for i in self.res_seq_lst:
            if not i in order_duplicate:
                order_duplicate.append(i)
        self.res_seq_lst=order_duplicate

        for seq in self.res_seq_lst:
            one_res_atms = []
            for atm in self.atom_lst:          
                if atm.res_seq == seq:
                    one_res_atms.append(atm)
            oneres = RESIDUE(seq, one_res_atms)
            self.residue_lst.append(oneres)

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
                    if res.seq+1 == self.residue_lst[next_res_idx].seq:
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
                    if res.seq+1 == self.residue_lst[next_res_idx].seq:
                        needadd=False
                    else:
                        needadd=True
                else:
                    pass
#         return self.stretch_dct
    
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
    
