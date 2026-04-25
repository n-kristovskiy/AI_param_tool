from collections import deque, defaultdict
from rdkit import Chem

def rdkit_pdb_modification_3(rdkit_mol, resname="MOD", resid=1, segid="A",
                             force_numeric=False, C_terminal=False):
    """
    Универсальная модификация атомных имён для аминокислот и модификаций.
    Правила:
    - Backbone: N, CA, C, O, OXT
    - HA назначается альфа-протону
    - Протоны аминогруппы: H1 H2 H3
    - Боковые атомы используют греческие уровни (B,G,D,E,Z,H,T,I,K,L,M,N)
    - Протоны наследуют имя родителя
    - Numeric fallback если длина имени >4
    """

    greek_levels = {1:"B",2:"G",3:"D",4:"E",5:"Z",6:"H",7:"T",8:"I",
                    9:"K",10:"L",11:"M",12:"N"}

    atom_names = {}
    used_names = set()
    visited = set()
    numeric_mode = force_numeric
    heavy_counter = 0

    # Найти backbone
    CA_idx, backbone = find_backbone_match(rdkit_mol, C_terminal)
    backbone_names = ["N","CA","C","O"]
    if len(backbone)==5:
        backbone_names=["N","CA","C","OC1","OC2"]

    for idx,name in zip(backbone,backbone_names):
        atom_names[idx] = name
        used_names.add(name)
        if rdkit_mol.GetAtomWithIdx(idx).GetSymbol() != "H":
            heavy_counter += 1
    N_idx = backbone[0]

    # HA
    for nbr in rdkit_mol.GetAtomWithIdx(CA_idx).GetNeighbors():
        if nbr.GetSymbol() == "H":
            atom_names[nbr.GetIdx()] = "HA"
            used_names.add("HA")

    # H аминогруппы
    amine_counter = 1
    for nbr in rdkit_mol.GetAtomWithIdx(N_idx).GetNeighbors():
        if nbr.GetSymbol() == "H":
            name = f"H{amine_counter}"
            atom_names[nbr.GetIdx()] = name
            used_names.add(name)
            amine_counter += 1

    # BFS боковой цепи
    queue = deque([(CA_idx, 0)])
    visited.add(CA_idx)
    branch_counter = defaultdict(int)

    while queue:
        current_idx, depth = queue.popleft()
        atom = rdkit_mol.GetAtomWithIdx(current_idx)
        heavy_neighbors = []
        hydrogens = []

        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in backbone or nbr.GetIdx() in visited:
                continue
            if nbr.GetSymbol() == "H":
                hydrogens.append(nbr)
            else:
                heavy_neighbors.append(nbr)

        for i, nbr in enumerate(heavy_neighbors, start=1):
            idx = nbr.GetIdx()
            visited.add(idx)
            branch_counter[current_idx] += 1
            branch_id = branch_counter[current_idx]

            element = nbr.GetSymbol()
            name = None

            # Greek mode
            if not numeric_mode:
                if depth+1 in greek_levels:
                    greek = greek_levels[depth+1]
                    candidate = f"{element}{greek}" if branch_id==1 else f"{element}{greek}{branch_id}"
                    if len(candidate) <= 4 and candidate not in used_names:
                        name = candidate
                    else:
                        numeric_mode = True
                else:
                    numeric_mode = True

            # Numeric mode
            if numeric_mode:
                heavy_counter += 1
                name = f"{element}{heavy_counter}"
                if len(name) > 4:
                    raise ValueError("Atom index overflow in numeric mode")

            atom_names[idx] = name
            used_names.add(name)
            queue.append((idx, depth+1))

            # Присвоение водородов сразу после heavy атома
            h_counter = 1
            for h in nbr.GetNeighbors():
                if h.GetSymbol() != "H" or h.GetIdx() in atom_names:
                    continue
                if numeric_mode:
                    parent_index = "".join(c for c in name if c.isdigit())
                    hname = f"H{parent_index}{h_counter}"
                else:
                    suffix = name[1:]
                    hname = f"H{suffix}" if len(heavy_neighbors)==1 else f"H{suffix}{h_counter}"
                if len(hname) > 4:
                    raise ValueError("Hydrogen name overflow")
                atom_names[h.GetIdx()] = hname
                used_names.add(hname)
                h_counter += 1

        # Водороды текущего атома (например, CA)
        if current_idx in atom_names and hydrogens:
            parent_name = atom_names[current_idx]
            h_counter = 1
            for h in hydrogens:
                if h.GetIdx() in atom_names:
                    continue
                if numeric_mode:
                    parent_index = "".join(c for c in parent_name if c.isdigit())
                    hname = f"H{parent_index}{h_counter}"
                else:
                    suffix = parent_name[1:]
                    hname = f"H{suffix}" if len(hydrogens)==1 else f"H{suffix}{h_counter}"
                if len(hname) > 4:
                    raise ValueError("Hydrogen name overflow")
                atom_names[h.GetIdx()] = hname
                used_names.add(hname)
                h_counter += 1

    # Запись в PDB
    for atom in rdkit_mol.GetAtoms():
        idx = atom.GetIdx()
        if idx not in atom_names:
            raise ValueError(f"Atom {idx} was not assigned a name")
        name = atom_names[idx]
        pdb_name = name[:4].ljust(4)

        info = Chem.AtomPDBResidueInfo()
        info.SetName(pdb_name)
        info.SetResidueName(resname)
        info.SetResidueNumber(resid)
        info.SetChainId(segid)
        atom.SetProp("AtomName", name.strip())
        atom.SetMonomerInfo(info)

    check_duplicate_atom_names(rdkit_mol)
    return rdkit_mol