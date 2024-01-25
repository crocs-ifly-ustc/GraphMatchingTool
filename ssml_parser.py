import utils
import re
import networkx as nx
from chemfig_struct import *
import argparse


match_virtual = {'}':'{', ')':'(', 'branch)':'branch('}

def add_atom(atom_dict, node_tag, name=""):
    if node_tag not in atom_dict:
        atom_dict[node_tag] = Atom(name)
    else:
        if atom_dict[node_tag].m_text != "":
            if name != "":
                print("conflict atom old = {} new = {}".format(atom_dict[node_tag].m_text, name))
                # pdb.set_trace()
                # raise ValueError("conflict atom old = {} new = {}".format(atom_dict[node_tag].m_text, name))
        else:
            atom_dict[node_tag].m_text = name


def get_atom_group(item_list):
    '''scan atom group'''
    lengths = [len(item) for item in item_list if isinstance(item, str)]

    start_index = None
    end_index = None

    new_item_list = []
    for i, item in enumerate(item_list):
        item_type = judge_str_item_type(item)
        if item_type == "atom":
            if start_index is None:
                start_index = i
                new_item_list.append(item)
            else:
                new_item_list[-1] += item
        else:
            if start_index is not None:
                start_index = None
            new_item_list.append(item)
    
    return new_item_list

def judge_str_item_type(item: str):
    if '?' in item:
        begin_conn_pattern = re.compile(r'\?\[[a-zA-Z]+\]')
        begin_result = begin_conn_pattern.findall(item)
        if len(begin_result) > 0:
            return 'reconn_begin'
        else:
            return 'reconn_end'

    for bond in bond_types:
        if item.startswith(bond+"["):
            return 'bond_atom'
    
    atom_pattern = re.compile(r'[a-zA-Z]+|\\circle')
    atom_result = atom_pattern.findall(item)
    if len(atom_result) > 0:
        if atom_result[0] == item:
            return 'atom'
    
    if item == 'branch(':
        return 'branch_begin'
    if item == 'branch)':
        return 'branch_end'

    return 'atom'

def attr_obtain(str):
    item_type = judge_str_item_type(str)
    if item_type == 'bond_atom':
        # bond_type bond_angle
        bond_type_pattern = re.compile(r'.*\[')
        bond_angle_pattern = re.compile(r'[\d\.]+')
        bond_type = bond_type_pattern.findall(str)[0][:-1]

        attr_list = bond_angle_pattern.findall(str)
        bond_angle = float(attr_list[0])
        bond_length = float(attr_list[1]) if len(attr_list) >= 2 else 1
        return bond_type, bond_angle, bond_length
    elif item_type == 'atom':
        # atom name
        return str
    elif item_type == 'reconn_end':
        # reconn bond type
        bond_type_pattern = re.compile(r'\{(.+)\}')
        bond_type = bond_type_pattern.findall(str)[0]
        assert bond_type is not None, "bond_type is null."
        return bond_type
    else:
        print("Error in attr_obtain, type not defined.")
        sys.exit()

def parse_ssml(inputStr, is_debug=False):
    # traverse every atom, v
    item_list = inputStr.split()
    item_list = get_atom_group(item_list)
    if is_debug:
        print(item_list)
    virtual_stack = [] # bracket stack, used to match brackets

    cur_atom = None
    cur_bond = None

    reconn_begin_atom_dict = {} # start of reconnect
    branch_stack = [None for i in range(1000)] # stack for backtrack
    stack_level = 0
    is_reconn = False # reconnection mark
    is_branch_end = False
    cur_reconn_tag = ''

    node_tag = 0
    branch_begin_tag = 0

    atom_dict = {}
    
    for ssml_item in item_list:
        ssml_item_type = judge_str_item_type(ssml_item)
        if is_debug:
            print("cur: ", ssml_item, ssml_item_type, stack_level)
        if ssml_item_type == 'atom':
            add_atom(atom_dict, node_tag, ssml_item)
        elif ssml_item_type == 'bond_atom':
            add_atom(atom_dict, node_tag + 1)

            bond_type, bond_angle, bond_length = attr_obtain(ssml_item)
            bond = Bond(bond_type)
            bond.m_angle = bond_angle
            bond.m_length = bond_length
            bond.end_atom = atom_dict[node_tag + 1]
            atom_dict[node_tag + 1].in_bonds.append(bond)
           
            if is_branch_end:
                branch_begin_tag = branch_stack[tgt_stack_level]
                is_branch_end = False
                bond.begin_atom = atom_dict[branch_begin_tag]
                atom_dict[branch_begin_tag].out_bonds.append(bond)
            else:
                add_atom(atom_dict, node_tag)
                bond.begin_atom = atom_dict[node_tag]
                atom_dict[node_tag].out_bonds.append(bond)
            node_tag += 1

        elif ssml_item_type == 'reconn_begin':
            is_reconn = True
            tag = re.match("\?\[([a-zA-Z]+)\]",ssml_item).group(1)
            cur_reconn_tag = tag
            reconn_begin_atom_dict[tag] = node_tag

        elif ssml_item_type == 'reconn_end':
            cur_reconn_tag, bond_type = re.match("\?\[([a-zA-Z]+)[,]\{([^\{\}]+)\}\]",ssml_item).groups()
            reconn_atom = reconn_begin_atom_dict[cur_reconn_tag] # get start of reconnection

            bond = Bond(bond_type)
            bond.begin_atom = atom_dict[reconn_atom]
            bond.end_atom = atom_dict[node_tag]
            atom_dict[reconn_atom].out_bonds.append(bond)
            atom_dict[node_tag].in_bonds.append(bond)
            bond.m_angle = None
            bond.m_length = None
            
        elif ssml_item_type == 'branch_begin':
            if not is_branch_end:
                branch_stack[stack_level] = node_tag # begin branch, push stack
            else:
                pass
            stack_level += 1
            
        elif ssml_item_type == 'branch_end':
            stack_level -= 1
            tgt_stack_level = stack_level
            is_branch_end = True

        elif ssml_item_type == 'virtual':
            if len(virtual_stack) == 0:
                virtual_stack.append(ssml_item) # put into stack
            else:
                cur_virtual = virtual_stack[-1] # stack top
                if ssml_item not in match_virtual.keys(): # left bracket, put info stack
                    virtual_stack.append(ssml_item)
                else: # right bracket, match
                    virtual_stack.pop() # pop from stack
    return atom_dict[0] if len(atom_dict) > 0 else None


def main(args):
    import tqdm
    import image_render
    import viz
    import shutil
    import cv2
    
    with open(args.input, "r", encoding="utf-8") as fin:
        lines = fin.readlines()
    

    for ind, line in enumerate(tqdm.tqdm(lines)):
        rel_img_path, ssml_lab = line.split("\t")
        img_path = os.path.join(args.img_prefix, rel_img_path)
        raw_img = cv2.imread(img_path)
        # shutil.copy(img_path, "raw.jpg")
        new_text, replace_dict, text = utils.replace_chemfig(ssml_lab)
        viz_img_arr = []
        for key, chemfig_str in replace_dict.items():
            root_atom = parse_ssml(" ".join(chemfig_str.split(" ")[2:-1]))
            viz_img = image_render.rend(root_atom, scale=100)
            viz_img_arr.append(viz_img)
        viz_img_all = viz.hstack_images(viz_img_arr)
        out_img = viz.vstack_images([raw_img, viz_img_all])
        
        output_img_path = os.path.join(args.output, rel_img_path)
        os.makedirs(os.path.dirname(output_img_path), exist_ok=True)
        cv2.imwrite(output_img_path, out_img)
        
        # pdb.set_trace()

    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser("")
    parser.add_argument("-input", type=str, default="test_decode_result.txt")
    parser.add_argument("-img_prefix", type=str, default="./data/mini_valid/mini_valid")
    parser.add_argument("-output", type=str, default="./viz_test_decode_result/")
    parser.add_argument("-num_workers", type=int, default=32)
    args = parser.parse_args()
    main(args)