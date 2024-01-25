from chemfig_struct import *

def GetCoordRange(all_atoms):
    min_x = math.inf
    min_y = math.inf
    max_x = -math.inf
    max_y = -math.inf
    for atom in all_atoms:
        min_x = min(atom.pos_x, min_x)
        min_y = min(atom.pos_y, min_y)
        max_x = max(atom.pos_x, max_x)
        max_y = max(atom.pos_y, max_y)
    return min_x, min_y, max_x, max_y

def NormFlatChemfig(inAtom:Atom): 
    all_atoms = GetAllAtoms(inAtom)
    out_text = []
    all_atoms = sorted(all_atoms, key=lambda x: x.pos_x)
    last_bond = None
    for atom_id, atom in enumerate(all_atoms):
        all_bonds = atom.in_bonds + atom.out_bonds
        if last_bond is not None:
            if last_bond not in all_bonds:
                return inAtom
            else:
                all_bonds.remove(last_bond)
        if len(all_bonds) <= 0 and atom_id < len(all_atoms) - 1:
            return inAtom
        if len(all_bonds) > 1:
            return inAtom
        if atom_id < len(all_atoms) - 1:
            if len(all_bonds) != 1:
                return inAtom
            new_bond = all_bonds[0]
            if new_bond.begin_atom == atom:
                if new_bond.end_atom != all_atoms[atom_id+1]:
                    return inAtom
                delta = math.fabs(new_bond.m_angle - 0) % 360
                if delta > 180:
                    delta = 360 - delta
                if delta > 15:
                    return inAtom
            else:
                if new_bond.begin_atom != all_atoms[atom_id+1]:
                    return inAtom
                delta = math.fabs(new_bond.m_angle - 180) % 360
                if delta > 180:
                    delta = 360 - delta
                if delta > 15:
                    return inAtom
        else:
            if len(all_bonds) != 0:
                return inAtom
            new_bond = None
        # output text
        if last_bond is not None and last_bond.m_length > 0.25:
            if last_bond.m_type not in bond_in_out_dict:
                return inAtom
            type_str = bond_in_out_dict[last_bond.m_type]
            out_text.append(type_str)
        out_text += atom.m_text
        last_bond = new_bond
    return out_text

def GetAllAtoms(atom:Atom):
    all_atoms = set()
    atom_stack = [atom]
    while len(atom_stack) > 0:
        curAtom = atom_stack.pop()
        if curAtom in all_atoms:
            continue
        if curAtom is None:
            continue
        all_atoms.add(curAtom)
        for bond in curAtom.out_bonds:
            atom_stack.append(bond.end_atom)
        for bond in curAtom.in_bonds:
            atom_stack.append(bond.begin_atom)
    all_atoms = list(all_atoms)
    all_atoms = sorted(all_atoms, key = lambda x:x.name)
    return list(all_atoms)

def SimulateCoord(rootAtom: Atom, scale=1):
    atom_stack = [(rootAtom, None)]
    all_atoms = []
    none_bonds = set()
    while len(atom_stack):
        curAtom, lastBond = atom_stack.pop(0)
        if curAtom is None:
            continue
        if curAtom in all_atoms:
            continue
        if lastBond is None:
            curAtom.pos_x = 0
            curAtom.pos_y = 0
        else:
            if lastBond.m_angle is None:
                assert lastBond.m_length is None
                none_bonds.add(lastBond)
                if lastBond.begin_atom.pos_x is not None and lastBond.end_atom.pos_x is not None:
                    # pdb.set_trace()
                    lastBond.m_angle = math.atan2(-lastBond.end_atom.pos_y + lastBond.begin_atom.pos_y, lastBond.end_atom.pos_x - lastBond.begin_atom.pos_x) * 180.0 / math.pi
                    lastBond.m_length = math.sqrt(math.pow(-lastBond.end_atom.pos_y + lastBond.begin_atom.pos_y, 2) + math.pow(lastBond.end_atom.pos_x - lastBond.begin_atom.pos_x, 2))
                else:
                    atom_stack.append((curAtom, lastBond))
                    continue
            elif lastBond.end_atom == curAtom:
                lastAtom = lastBond.begin_atom
                try:
                    curAtom.pos_x = lastAtom.pos_x + lastBond.m_length * scale * math.cos(lastBond.m_angle * math.pi / 180)
                    curAtom.pos_y = lastAtom.pos_y - lastBond.m_length * scale * math.sin(lastBond.m_angle * math.pi / 180)
                except BaseException as e:
                    pdb.set_trace()
            elif lastBond.begin_atom == curAtom:
                lastAtom = lastBond.end_atom
                curAtom.pos_x = lastAtom.pos_x - lastBond.m_length * scale * math.cos(lastBond.m_angle * math.pi / 180)
                curAtom.pos_y = lastAtom.pos_y + lastBond.m_length * scale * math.sin(lastBond.m_angle * math.pi / 180)

        all_atoms.append(curAtom)

        for bond in curAtom.out_bonds:
            if bond.end_atom is not None:
                atom_stack.append((bond.end_atom, bond))
            else:
                end_atom = Atom("")
                bond.end_atom = end_atom
                end_atom.in_bonds.append(bond)
                atom_stack.append((bond.end_atom, bond))
        for bond in curAtom.in_bonds:
            if bond.begin_atom is not None:
                atom_stack.append((bond.begin_atom, bond))
            else:
                begin_atom = Atom("")
                bond.begin_atom = begin_atom
                begin_atom.out_bonds.append(bond)
                atom_stack.append((bond.begin_atom, bond))
    for bond in none_bonds:
        if bond.begin_atom.pos_x is not None and bond.end_atom.pos_x is not None:
            bond.m_angle = math.atan2(-bond.end_atom.pos_y + bond.begin_atom.pos_y, bond.end_atom.pos_x - bond.begin_atom.pos_x) * 180.0 / math.pi
            bond.m_length = math.sqrt(math.pow(-bond.end_atom.pos_y + bond.begin_atom.pos_y, 2) + math.pow(bond.end_atom.pos_x - bond.begin_atom.pos_x, 2))

    return all_atoms

# connect circle atom to all side atoms
def NormCircleAtom(circleAtom: Atom, base=1e10, all_connect=0):
    if len(circleAtom.out_bonds) +  len(circleAtom.in_bonds) != 1:
        return

    if len(circleAtom.out_bonds) == 1: #circle -> begin
        begin_atom = circleAtom.out_bonds[0].end_atom
        origin_length = circleAtom.out_bonds[0].m_length
        begin_atom.in_bonds.remove(circleAtom.out_bonds[0])
        circleAtom.out_bonds.clear()
    else: #circle <- begin
        begin_atom = circleAtom.in_bonds[0].begin_atom
        origin_length = circleAtom.in_bonds[0].m_length
        begin_atom.out_bonds.remove(circleAtom.in_bonds[0])
        circleAtom.in_bonds.clear()

    ref_angle = math.atan2(-begin_atom.pos_y + circleAtom.pos_y, begin_atom.pos_x - circleAtom.pos_x) * 180.0 / math.pi
    ref_length = math.sqrt(math.pow(begin_atom.pos_y - circleAtom.pos_y, 2) + math.pow(begin_atom.pos_x - circleAtom.pos_x, 2))
    visited = set()
    path = []
    
    def dfs_find(curAtom, lastBond=None, visited=[]):
        if curAtom in visited:
            return [visited]
        cur_visited = visited + [curAtom]
        cur_angle = math.atan2(-curAtom.pos_y + circleAtom.pos_y, curAtom.pos_x - circleAtom.pos_x) * 180.0 / math.pi
        ret = []
        cand_bonds = []
        for bond in curAtom.out_bonds + curAtom.in_bonds:
            next_atom = bond.end_atom if bond.end_atom != curAtom else bond.begin_atom
            if len(visited)>0 and next_atom == visited[-1]:
                continue
            if next_atom.m_text == "\\circle":
                continue
            next_angle = math.atan2(-next_atom.pos_y + circleAtom.pos_y, next_atom.pos_x - circleAtom.pos_x) * 180.0 / math.pi
            next_length = math.sqrt(math.pow(next_atom.pos_y - circleAtom.pos_y, 2) + math.pow(next_atom.pos_x - circleAtom.pos_x, 2))
            delta_length = math.fabs(next_length - ref_length)
            cand_bonds.append((bond, next_atom, next_angle, next_length, delta_length))
        cand_bonds = sorted(cand_bonds, key=lambda x:x[-1])
        for bond, next_atom, next_angle, next_length, delta_length in cand_bonds:
            if (next_angle-cur_angle) % 360 >= 180:
                continue
            if math.fabs(next_length - ref_length) > ref_length * 0.5:
                continue
            ret += dfs_find(next_atom, bond, cur_visited)
        return ret
    
    det_results = dfs_find(begin_atom, None, [])
    if len(det_results) <= 0:
        det_results = [[begin_atom]]

    ring_atoms = det_results[0]
    ring_atoms = sorted(ring_atoms, key = lambda x:x.pos_x*base - x.pos_y)
    selected_atom = ring_atoms[0]

    for selected_atom in ring_atoms:
        new_bond = Bond("-:")
        new_bond.begin_atom = selected_atom
        new_bond.end_atom = circleAtom
        new_bond.m_angle = math.atan2( - new_bond.end_atom.pos_y + new_bond.begin_atom.pos_y, new_bond.end_atom.pos_x - new_bond.begin_atom.pos_x) * 180.0 / math.pi
        new_bond.m_length = origin_length
        selected_atom.out_bonds.append(new_bond)
        circleAtom.in_bonds.append(new_bond)
        if not all_connect:
            break

    return

def NormAllCircleAtom(rootAtom:Atom, all_connect=0):
    all_atoms = GetAllAtoms(rootAtom)
    min_x, min_y, max_x, max_y = GetCoordRange(all_atoms)
    max_value = max(max_y-min_x, max_y-min_y)
    for atom in all_atoms:
        if atom.m_text == "\\circle":
            NormCircleAtom(atom, base=max_value, all_connect=all_connect)
    return all_atoms
