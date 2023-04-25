# pdb_preparer
2022/11/21  Add RESIDUE.convert_MOL_atom() function for rename of the atom_name of the atom with the same element.
2022/11/21 	Add the code for convert ZNA's atom_name to 'Zn' in RESIDUE()
2022/11/21  Add the code for recognize the SO4 residue.
2022/11/22  Now can do the correct TER writing.
2023/04/21  Now generate RESIDUE() object by recognizing residue names, residue squence and chainID at the same time. 
            Now can divide strench by calculating distance from C atom in the current residue to N atom in the next residue. 
