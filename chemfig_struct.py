import os, sys
import re
import pdb
import numpy as np
import math
import warnings
import copy

bond_types = ["-", "=", "~", ">", "<", ">:", "<:", ">|", "<|", "-:", "=_", "=^", "~/"]
bond_types = sorted(bond_types, key=lambda x: -len(x))

bond_in_out_dict = {"-": "-", "=": "=", "~": "\equiv", "-:": ""}

directed_bond_types = {
    ">": "<",
    "<": ">",
    ">:": "<:",
    "<:": ">:",
    ">|": "<|",
    "<|": ">|",
}


class Atom(object):
    index = 0

    def __init__(self, text=""):
        self.name = "Atom_{}".format(Atom.index)
        Atom.index += 1
        self.m_text = text
        self.pos_x = None
        self.pos_y = None
        self.in_bonds = []
        self.out_bonds = []
        pass

    def __repr__(self):
        out_str = "<Atom"
        out_str += " name={}".format(self.name)
        out_str += " text={}".format(self.m_text)
        out_str += " x={}".format(self.pos_x)
        out_str += " y={}".format(self.pos_y)
        out_str += ">"
        return out_str

    # def norm_chem_above_below(self, text_arr):
    #     rootNode = transcription.parse_tree(text_arr)
    #     node_stack = [rootNode]
    #     cursor = 0
    #     while len(node_stack) > 0 and cursor < len(node_stack):
    #         curNode = node_stack[cursor]
    #         cursor += 1
    #         for word_id, word in enumerate(curNode.words):
    #             if type(word) is Node:
    #                 node_stack.append(word)
    #             elif word == "\\chemabove":
    #                 curNode.words[word_id] = "\\Chemabove"
    #                 pass
    #             elif word == "\\Chembelow" or word == "\\chembelow":
    #                 curNode.words[word_id] = "\\Chemabove"
    #                 # swap
    #                 first = curNode.words[word_id+1]
    #                 second = curNode.words[word_id+2]
    #                 curNode.words[word_id+1] = second
    #                 curNode.words[word_id+2] = first
    #     rootNode.fix_op_nest()
    #     return rootNode.flatten()

    def normed_text(self):
        return self.m_text
        # if self.m_text is not None and self.m_text != "":
        #     text_arr = self.m_text.split(" ")
        #     text_arr = self.norm_chem_above_below(text_arr)
        #     return " ".join(text_arr)
        # else:
        #     return self.m_text

class Bond(object):
    index = 0
    __default__ = {"m_angle": 0, "m_length": 1, "m_start": 0, "m_end": 0}

    def __init__(self, b_type="-"):
        self.name = "Bond_{}".format(Bond.index)
        Bond.index += 1
        b_type = b_type.replace("_", "").replace("^", "")

        self.m_type = b_type
        self.m_angle = None
        self.m_length = 1
        self.begin_atom = None
        self.end_atom = None

        self.ring_ids = {}


    def __repr__(self):
        out_str = "<Bond"
        out_str += " name={}".format(self.name)
        out_str += " type={}".format(self.m_type)
        out_str += " angle={}".format(self.m_angle)
        out_str += " length={}".format(self.m_length)
        out_str += " begin_atom={}".format(self.begin_atom.name if self.begin_atom else "none")
        out_str += " end_atom={}".format(self.end_atom.name if self.end_atom else "none")
        out_str += ">"
        return out_str