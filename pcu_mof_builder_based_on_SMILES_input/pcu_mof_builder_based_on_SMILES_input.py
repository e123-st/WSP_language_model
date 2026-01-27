from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pormake as pm
import os

def MOFbuilder_error_fixed(database_path:str,MOF_output_path:str,O_edge=None,N_edge=None):
     
    node_name = 'node'
    
    cifs = os.listdir(MOF_output_path)
    if O_edge != None :
       for i in range(len(O_edge)):
           O_linker = 'O_edge_'+str(O_edge[i])
           for j in range(len(cifs)):
               O_name=cifs[j].strip('.cif').split('-')[1]
               if O_name == O_linker :
                  N_linker = cifs[j].strip('.cif').split('-')[2]                  
                  MOF_output_name=cifs[j].strip('.cif')
                  MOF_Build(database_file_path=database_path,
                            output_file_path=MOF_output_path,
                            output_name=MOF_output_name,
                            node=node_name,
                            edge1=O_linker,
                            edge2=N_linker)
            
    if N_edge != None :
       for i in range(len(N_edge)):
           N_linker = 'N_edge_'+str(N_edge[i])
           for j in range(len(cifs)):
               N_name=cifs[j].strip('.cif').split('-')[2]               
               if N_name == N_linker :                 
                  O_linker = cifs[j].strip('.cif').split('-')[1]                  
                  MOF_output_name=cifs[j].strip('.cif')
                  MOF_Build(database_file_path=database_path,
                            output_file_path=MOF_output_path,
                            output_name=MOF_output_name,
                            node=node_name,
                            edge1=O_linker,
                            edge2=N_linker)        

    return print('--------Finish!--------')


def MOFbuilder(mofid_path:str,database_path:str,MOF_output_path:str,N_bond_length:float):
    
    node_name = 'node'
    #node_path = database_path+'/node' 
    O_linker_path = database_path+'/O_edge'
    N_linker_path = database_path+'/N_edge'
    if os.path.exists(O_linker_path) == False :
       os.mkdir(O_linker_path) 
    if os.path.exists(N_linker_path) == False :
       os.mkdir(N_linker_path) 
    
    if os.path.exists(MOF_output_path) == False :
       os.mkdir(MOF_output_path)
              
    O_linker = []
    N_linker = []
    aromatic_n_1 = []
    MOFid = []
    txtfile = open(mofid_path,'r').readlines()
    for i  in range(len(txtfile)):
        if '&&' not in txtfile[i]:
           mofid = txtfile[i].strip('\n').split(' ')[0]
        if '&&' in txtfile[i]:
           mofid = txtfile[i].strip('\n').split('&&')[0]   
        MOFid.append(mofid)
        mofid = mofid.split('.')
        for j in range(1,len(mofid)):
            if '[O-]' in mofid[j] :
               O_linker.append(mofid[j])
            else:
                N_linker.append(mofid[j])
    
    O_linker = np.unique(O_linker)      
    N_linker = np.unique(N_linker)  
    
    O_linker_database_output_path = database_path+'/O_linker.txt'
    N_linker_database_output_path = database_path+'/N_linker.txt'
    aromatic_n_1_output_path = database_path+'/N_linker_potential_error.txt'

    output=open(O_linker_database_output_path,'w') 
    for i in range(len(O_linker)):
        output.write(str(i+1)+':  '+O_linker[i]+'\n')
    output.close()
    
    output=open(N_linker_database_output_path,'w') 
    for i in range(len(N_linker)):
        output.write(str(i+1)+':  '+N_linker[i]+'\n')
    output.close()
    
    for i in range(len(O_linker)):
        SMILES = O_linker[i]
        xyz_name = 'O_edge_' + str(i+1)
        output_path = O_linker_path        
        SMILES_TO_XYZ_FOR_PORMAKE(SMILES,xyz_name,output_path,N_bond_length)
        
        
    for i in range(len(N_linker)):
        SMILES = N_linker[i]
        xyz_name = 'N_edge_' + str(i+1)
        output_path = N_linker_path        
        XYZ=SMILES_TO_XYZ_FOR_PORMAKE(SMILES,xyz_name,output_path,N_bond_length)
        if len(XYZ) != 0 :
           aromatic_n_1.append(XYZ[0])
        
    output=open(aromatic_n_1_output_path,'w') 
    for i in range(len(aromatic_n_1)):
        output.write(aromatic_n_1[i]+'\n')
    output.close()    
        
    for i in range(len(MOFid)):
        mofid = MOFid[i].split('.')
        O = 0
        for j in range(1,len(mofid)):
            linker = mofid[j]   
            if O == 0 :
               for k in range(len(O_linker)):                
                   if linker == O_linker[k] :
                      O_edge = 'O_edge_' + str(k+1)
                      O = O+1
                      break
                   if k == len(O_linker)-1 and linker != O_linker[k] :
                      for l in range(len(N_linker)):
                          if linker == N_linker[l] :
                             N_edge = 'N_edge_' + str(l+1)
                             break
            if O == 1 :
               for k in range(len(N_linker)):                
                   if linker == N_linker[k] :
                      N_edge = 'N_edge_' + str(k+1)
                      break                      
        
        MOF_output_name=str(i+1)+'-'+O_edge+'-'+N_edge
        MOF_Build(database_file_path=database_path,
                  output_file_path=MOF_output_path,
                  output_name=MOF_output_name,
                  node=node_name,
                  edge1=O_edge,
                  edge2=N_edge)     
        
    return print('--------Finish!--------')


def MOF_Build(database_file_path:str,
              output_file_path:str,
              output_name:str,
              node:str,
              edge1:str,#[O-]
              edge2:str):
        
    #软件自带数据库
    database = pm.Database()
    builder=pm.Builder()
    
    node_path = database_file_path+'/node/'+node+'.xyz'
    edge1_path = database_file_path+'/O_edge/'+edge1+'.xyz'
    edge2_path = database_file_path+'/N_edge/'+edge2+'.xyz'

    node = pm.BuildingBlock(node_path)
    edge1 = pm.BuildingBlock(edge1_path)
    edge2 = pm.BuildingBlock(edge2_path)
    topology = database.get_topo("pcu")

    node_bbs = {0:node}
    edge_bbs = {(0,0):edge1}

    bbs = builder.make_bbs_by_type(topology=topology, node_bbs=node_bbs, edge_bbs=edge_bbs)

    bbs[3]=edge2.copy()

    MOF = builder.build(topology=topology,bbs=bbs)

    cif_save_path = output_file_path+'/'+output_name+'.cif'
    
    return MOF.write_cif(cif_save_path)


def read_bond(mol):
    
    bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == 1.0 :
           bonds.append(str(bond.GetBeginAtomIdx())+' '+str(bond.GetEndAtomIdx())+' S')   
        if bond.GetBondTypeAsDouble() == 1.5 :
           bonds.append(str(bond.GetBeginAtomIdx())+' '+str(bond.GetEndAtomIdx())+' A')
        if bond.GetBondTypeAsDouble() == 2.0 :
           bonds.append(str(bond.GetBeginAtomIdx())+' '+str(bond.GetEndAtomIdx())+' D') 
        if bond.GetBondTypeAsDouble() == 3.0 :
           bonds.append(str(bond.GetBeginAtomIdx())+' '+str(bond.GetEndAtomIdx())+' T') 
    
    return bonds


def SMILES_TO_XYZ_FOR_PORMAKE(smiles,xyz_name,output_path,N_bond_length=1.0):
    
    output_path = output_path+'/'+xyz_name+'.xyz'
    aromatic_n_1 = []
    
    SMILES = smiles
    
    if '[C]' in SMILES:
       SMILES = SMILES.replace('[C]','C')
    
    if '[O-]C(=O)' not in SMILES and 'C(=O)[O-]' not in SMILES :
        
       mol = Chem.MolFromSmiles(SMILES)
       mol = Chem.AddHs(mol)
       AllChem.EmbedMolecule(mol)
       AllChem.UFFOptimizeMolecule(mol)       

       bonds = read_bond(mol) 
       
       nh_atoms = []       
       neighbor_h_atoms = []
       n_no_h_atoms = []     
       aromatic_n = []
       aromatic_ring_n = []
       ring_info = mol.GetRingInfo()
       atom_rings = ring_info.AtomRings()
       #n_aromatic_ring = []
       for atom in mol.GetAtoms():
           if atom.GetIsAromatic() and atom.GetAtomicNum() == 7:#芳香N
              aromatic_n.append(atom.GetIdx())
           if atom.GetAtomicNum() == 7 :                
              if atom.IsInRingSize(5) == True or atom.IsInRingSize(6) == True :#在5/6元环
                 for ring in atom_rings :
                     for atom_idx in ring :
                         if atom_idx == atom.GetIdx() :
                            R = ring
                            #n_aromatic_ring.append(R)
                            aromatic_ring = True
                            for idx in R :
                                if len(mol.GetAtomWithIdx(idx).GetNeighbors()) > 3 : #平面 最多3个原子
                                   aromatic_ring = False 
                                   break
                            if aromatic_ring == True :
                               aromatic_ring_n.append(atom.GetIdx()) 

       #n_aromatic_ring = np.unique(n_aromatic_ring)
                                    
       if len(aromatic_n) >= 2 and len(aromatic_n) == len(aromatic_ring_n) :  #rdkit识别的芳香N和芳香环数上n一致，潜在配位N已确定，芳香N不会连O
          for atom in mol.GetAtoms():       
            # 检查该原子是否连接了H原子
              has_h = False
              has_O = False
              #nh = 0 
              if atom.GetAtomicNum() == 7 :
                 for i in range(len(aromatic_n)):
                     if atom.GetIdx() == aromatic_n[i] :
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 1:
                               has_h = True   
                               #nh = nh+1               
                        if has_h == True : #最多只能连1个H                       
                           nh_atoms.append(atom.GetIdx())
                           neighbor_h_atoms.append(neighbor.GetIdx())
                        if has_h == False :  
                           n_no_h_atoms.append(atom.GetIdx())
       
       if len(aromatic_n) < 2 or len(aromatic_n) != len(aromatic_ring_n):               
          for atom in mol.GetAtoms():       
            # 检查该原子是否连接了H原子
              has_c = False
              has_h = False
              has_O = False
              nh = 0 
              nc = 0
              if atom.GetAtomicNum() == 7 :
                 for neighbor in atom.GetNeighbors():
                     if neighbor.GetAtomicNum() == 1:
                        has_h = True   
                        nh = nh+1  
                     if neighbor.GetAtomicNum() == 6:  
                        has_c = True 
                        nc = nc+1  
                 if has_h == True and nh == 1 : #最多只能连1个H 通常没有氢 1-2个C，3个C表明在内部不在边缘，不可配位，不直接连O
                    for neighbor in atom.GetNeighbors(): 
                        if neighbor.GetAtomicNum() == 8:
                           has_O = True   
                           break
                    if has_O == False and has_c == True and nc <= 2 : 
                       nh_atoms.append(atom.GetIdx())
                       neighbor_h_atoms.append(neighbor.GetIdx())
                 if has_h == False :  
                    for neighbor in atom.GetNeighbors(): 
                        if neighbor.GetAtomicNum() == 8:
                           has_O = True   
                           break 
                    if has_O == False and nc <= 2 :  
                       n_no_h_atoms.append(atom.GetIdx())           
                                                                                                   
       mol = pybel.readstring("smi", SMILES)
       mol.make3D()
       mol.write("xyz", output_path, overwrite=True)
                     
       X=[]
       N_atoms = []
       molecules=open(output_path,'r').readlines()   
       for i in range(1,len(molecules)):
           coordinate = molecules[i].split(' ')[0]
           if coordinate[0] == 'N':
              for j in range(len(nh_atoms)):
                  if i-2 == nh_atoms[j] :
                     N_atoms.append(i-2) #N原子的索引
                     break
                  
              for j in range(len(n_no_h_atoms)):
                  if i-2 == n_no_h_atoms[j] :
                     N_atoms.append(i-2) #N原子的索引
                     break    
                                   
       for i in range(len(N_atoms)):
           coordinate = molecules[N_atoms[i]+2].split(' ')
           for j in range(1,len(coordinate)):
               if coordinate[j] != '':            
                  x=coordinate[j]
                  X.append(eval(x))
                  break

       number_N = []
       right_N_x = np.max(X)
       left_N_x = np.min(X)
              
       if len(aromatic_n) != len(aromatic_ring_n) :          
          aromatic_x = []          
          n_N_number = smiles.count('n') + smiles.count('N')
          if n_N_number > 2 :
             aromatic_n_1.append(smiles)              
             for i in range(len(aromatic_ring_n)):
                 coordinate = molecules[aromatic_ring_n[i]+2].split(' ')
                 for j in range(1,len(coordinate)):
                     if coordinate[j] != '':            
                        x=eval(coordinate[j])
                        aromatic_x.append(x)
                        break      
             right_N_aromatic_x = np.max(aromatic_x)
             left_N_aromatic_x = np.min(aromatic_x)   
             if right_N_x != right_N_aromatic_x :
                right_N_x = right_N_aromatic_x 
             if left_N_x != left_N_aromatic_x :
                left_N_x = left_N_aromatic_x
                                                                         
       for i in range(len(X)):
           if X[i] == right_N_x:   
              right_N_number=N_atoms[i]
              number_N.append(N_atoms[i]) 
              coordinate = molecules[N_atoms[i]+2].split(' ')              
              for j in range(1,len(coordinate)):
                  if coordinate[j] != '':            
                     x=format(eval(coordinate[j])+N_bond_length,'.5f')
                     for k in range(j+1,len(coordinate)):
                         if coordinate[k] != '': 
                            y=coordinate[k]
                            break
                     break   
              z=coordinate[-1]
              X_write = 'X'+'         '+str(x)+'         '+y+'         '+z
              molecules.append(X_write)
              break   
                 
       for i in range(len(X)):
           if X[i] == left_N_x:
              left_N_number=N_atoms[i] 
              number_N.append(N_atoms[i]) 
              coordinate = molecules[N_atoms[i]+2].split(' ')
              for j in range(1,len(coordinate)):
                  if coordinate[j] != '':            
                     x=format(eval(coordinate[j])-N_bond_length,'.5f')
                     for k in range(j+1,len(coordinate)):
                         if coordinate[k] != '': 
                            y=coordinate[k]
                            break
                     break   
              z=coordinate[-1]
              X_write = 'X'+'         '+str(x)+'         '+y+'         '+z
              molecules.append(X_write)
              break  

       neighbor_coord_nh = [] 
       for i in range(len(nh_atoms)):
           if right_N_number == nh_atoms[i] or left_N_number == nh_atoms[i]:
              molecules[neighbor_h_atoms[i]+2]=-1
              neighbor_coord_nh.append(neighbor_h_atoms[i]) 
       
       if len(neighbor_coord_nh) == 0 :
          for i in range(len(bonds)) : 
              molecules.append(bonds[i]+'\n')
       
       if len(neighbor_coord_nh) != 0 : 
          neighbor_coord_nh.sort() 
          for i in range(len(bonds)) :
              bond = bonds[i].split(' ')           
              atom_number_0 = bond[0]
              atom_number_1 = bond[1] 
              
              for j in range(len(neighbor_coord_nh)):                                
                  if eval(atom_number_0) == neighbor_coord_nh[j] or eval(atom_number_1) == neighbor_coord_nh[j] :
                     break
                  if j == len(neighbor_coord_nh) - 1 and eval(atom_number_0) != neighbor_coord_nh[j] and eval(atom_number_1) != neighbor_coord_nh[j] :
                     if eval(atom_number_0) < neighbor_coord_nh[0]:
                        m_0 = str(eval(bond[0])) 
                     if eval(atom_number_0) > neighbor_coord_nh[-1]:
                        m_0 = str(eval(bond[0])-len(neighbor_coord_nh))
                     if eval(atom_number_0) > neighbor_coord_nh[0] and eval(atom_number_0) < neighbor_coord_nh[-1]:
                        for k in range(len(neighbor_coord_nh)-1):
                            if eval(atom_number_0) > neighbor_coord_nh[k] and eval(atom_number_0) < neighbor_coord_nh[k+1]:
                               m_0 = str(eval(bond[0])-k-1) 
                               break 
                        
                     if eval(atom_number_1) < neighbor_coord_nh[0]:
                        m_1 = str(eval(bond[1])) 
                     if eval(atom_number_1) > neighbor_coord_nh[-1]:
                        m_1 = str(eval(bond[1])-len(neighbor_coord_nh))
                     if eval(atom_number_1) > neighbor_coord_nh[0] and eval(atom_number_1) < neighbor_coord_nh[-1]:
                        for k in range(len(neighbor_coord_nh)-1):
                            if eval(atom_number_1) > neighbor_coord_nh[k] and eval(atom_number_1) < neighbor_coord_nh[k+1]:
                               m_1 = str(eval(bond[1])-k-1) 
                               break   
                        
                     molecules.append(m_0+' '+m_1+' '+bond[2]+'\n')
                                                 
       len_molecule = len(molecules)-len(bonds)-2   
                                        
       new_molecule=[]
       new_molecule.append(str(eval(molecules[0])+2-len(neighbor_coord_nh)))
       new_molecule.append('\n')
       for i in range(1,len(molecules)):
           if molecules[i] != -1:
              new_molecule.append(molecules[i])
       
       new_molecule.append(str(number_N[0])+' '+ str(len_molecule-2)+' S')  
       new_molecule.append('\n')
       new_molecule.append(str(number_N[1])+' '+ str(len_molecule-1)+' S')       

    if '[O-]C(=O)' in SMILES or 'C(=O)[O-]' in SMILES :  
                          
       #mol = Chem.MolFromSmiles(SMILES.replace('[O-]','[OH]'))
       mol = Chem.MolFromSmiles(SMILES)
       mol = Chem.AddHs(mol)
       AllChem.EmbedMolecule(mol)
       AllChem.UFFOptimizeMolecule(mol)

       bonds = read_bond(mol) 
       
       C_atoms = []       
       O_atoms = [] 
       for atom in mol.GetAtoms():
           if atom.GetAtomicNum() == 6:
              O = 0
              O_neighbor_number = 0
              for neighbor in atom.GetNeighbors():
                  if neighbor.GetAtomicNum() == 8:
                     O = O + 1 
                     O_neighbor_number = O_neighbor_number + len(neighbor.GetNeighbors())
                     
              if O == 2 and O_neighbor_number == 2 :       
                 C_atoms.append(atom.GetIdx())
                 for neighbor in atom.GetNeighbors():
                     if neighbor.GetAtomicNum() == 8:
                        O_atoms.append(neighbor.GetIdx())   
                        for O_neighbor in neighbor.GetNeighbors():
                            if O_neighbor.GetAtomicNum() == 1:
                               O_atoms.append(O_neighbor.GetIdx())                        
                               break    
                                               
                     
       output_mol = pybel.readstring("smi", SMILES)
       output_mol.make3D()
       output_mol.write("xyz", output_path, overwrite=True) 
              
       new_molecule=[]
       molecules=open(output_path,'r').readlines()
       new_molecule.append(str(eval(molecules[0])))   
       new_molecule.append('\n')
       new_molecule.append('\n')  
       
       N_COO = 0
       for i in range(2,len(molecules)):
           element = molecules[i].split(' ')[0]
           if element != 'C' and element != 'O' :
              new_molecule.append(molecules[i])
              
           if element == 'C' :
              for j in range(len(C_atoms)):                  
                  if i == C_atoms[j]+2 :
                     coordinate = molecules[i].split(' ')                                 
                     for k in range(1,len(coordinate)):
                         if coordinate[k] != ''  :            
                            x=coordinate[k]
                            for l in range(k+1,len(coordinate)):
                                if coordinate[l] != '': 
                                   y=coordinate[l]
                                   break                        
                            break 
                     z=coordinate[-1]
                     X_write = 'X'+'         '+x+'         '+y+'         '+z
                     new_molecule.append(X_write) 
                     N_COO = N_COO + 1
                     break
                  if j == len(C_atoms)-1 and i != C_atoms[j]+2 :
                     new_molecule.append(molecules[i]) 
                     
           if element == 'O' :
              for j in range(len(O_atoms)):
                  if i == O_atoms[j]+2 :                       
                     break 
                  if j == len(O_atoms)-1 and i != O_atoms[j]+2 :
                     new_molecule.append(molecules[i]) 
                               
                                                                            
       O_atoms.sort()
       #print(O_atoms)
       for i in range(len(bonds)) :
           bond = bonds[i].split(' ')           
           atom_number_0 = bond[0]
           atom_number_1 = bond[1] 
           
           for j in range(len(O_atoms)):
               if eval(atom_number_0) == O_atoms[j] or eval(atom_number_1) == O_atoms[j] :
                  break
               if j == len(O_atoms) - 1  and eval(atom_number_0) != O_atoms[j] and eval(atom_number_1) != O_atoms[j] :
                  if eval(atom_number_0) < O_atoms[0]:
                     m_0 = str(eval(bond[0])) 
                  if eval(atom_number_0) > O_atoms[-1]:
                     m_0 = str(eval(bond[0])-len(O_atoms))
                  if eval(atom_number_0) > O_atoms[0] and eval(atom_number_0) < O_atoms[-1]:
                     for k in range(len(O_atoms)-1):
                         if eval(atom_number_0) > O_atoms[k] and eval(atom_number_0) < O_atoms[k+1]:
                            m_0 = str(eval(bond[0])-k-1) 
                  if eval(atom_number_1) < O_atoms[0]:
                     m_1 = str(eval(bond[1])) 
                  if eval(atom_number_1) > O_atoms[-1]:
                     m_1 = str(eval(bond[1])-len(O_atoms))
                  if eval(atom_number_1) > O_atoms[0] and eval(atom_number_1) < O_atoms[-1]:
                     for k in range(len(O_atoms)-1):
                         if eval(atom_number_1) > O_atoms[k] and eval(atom_number_1) < O_atoms[k+1]:
                            m_1 = str(eval(bond[1])-k-1)           
                  new_molecule.append(m_0+' '+m_1+' '+bond[2]+'\n')   
                             
       new_molecule[0] = str(eval(molecules[0])-N_COO*2)   
                                                        
    output=open(output_path,'w')
    for i in range(len(new_molecule)):
        output.write(str(new_molecule[i]))
    output.close()
    
    return aromatic_n_1


#MOFbuilder(mofid_path=r"C:\Users\E\Desktop\mofid.txt",
#           database_path=r"C:\Users\E\Desktop\database",
#           MOF_output_path=r"C:\Users\E\Desktop\MOFbuild",
#           N_bond_length=1.0)
        
#MOFbuilder_error_fixed(database_path=r"C:\Users\E\Desktop\database",
#                       MOF_output_path=r"C:\Users\E\Desktop\MOFbuild",
#                       O_edge=None,
#                       N_edge=[1])
#SMILES_TO_XYZ_FOR_PORMAKE(SMILES='[O-]C(=O)C=CC(=O)[O-]',xyz_name='edge1',output_path=r'C:\Users\E\Desktop\cif\database')
#SMILES_TO_XYZ_FOR_PORMAKE(SMILES='n1ccc(cc1)C1=N[N]C(=N[N]1)c1ccncc1',xyz_name='edge2',output_path=r'C:\Users\E\Desktop\cif\database')

#MOF_Build(database_file_path=r'C:\Users\E\Desktop\cif\database',
#              output_file_path=r'C:\Users\E\Desktop\cif\database',
#              output_name='test',
#              node='node1',
#              edge1='edge1',
#              edge2='edge2')

