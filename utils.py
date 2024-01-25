from collections import OrderedDict
import re
from chemfig_struct import *
import Levenshtein

pair_dict = {"}":"{", "]":"["}
def replace_chemfig(text):
    replace_dict = OrderedDict()
    ind = 0
    new_text = ""
    while True:
        pos = text.find("\\chemfig")
        if pos == -1:
            break
        cur_pos = pos + 8
        cur_left_pair = None
        cur_left_pos = None
        curLevel = 0
        range_cnt = {"[":0, "{":0}
        while cur_pos < len(text):
            ch = text[cur_pos]
            if ch == "[" or ch == "{":
                if cur_left_pair is None:
                    cur_left_pair = ch
                    cur_left_pos = cur_pos
                    curLevel = 1
                elif cur_left_pair == ch:
                    curLevel += 1
            elif ch == "}" or ch == "]":
                if cur_left_pair == pair_dict[ch]:
                    curLevel -= 1
                    if curLevel == 0:
                        range_cnt[cur_left_pair] += 1
                        if range_cnt["["] > 1:
                            raise ValueError("multiple attr range")
                        if range_cnt["{"] >= 1:
                            # pdb.set_trace()
                            break
                        cur_left_pair = None
                        cur_left_pos = None
                    
            elif cur_left_pair is None:
                if ch != " ":
                    warnings.warn("format err, input = {}".format(text))
                    break
                    # raise ValueError("format err, input = {}".format(text))
            else:
                pass
            cur_pos += 1
        if cur_left_pos is None:
            text = text[pos+8:]
            continue
        beginPos = cur_left_pos
        endPos = cur_pos + 1
        rep_key = "\\chem{}".format(chr(ord('a') + ind))
        ind += 1
        replace_dict[rep_key] = "\\chemfig "+text[beginPos:endPos]
        text = text[0:pos] + " " + rep_key + " " + text[endPos:]
        new_text += replace_dict[rep_key] + " "

        pos = cur_pos + 1
    return new_text, replace_dict, text

def cal_edit_ops(str1, str2):
    char_idx_dict = dict()
    for item in str1:
        if item not in char_idx_dict:
            char_idx_dict[item] = chr(len(char_idx_dict))
    for item in str2:
        if item not in char_idx_dict:
            char_idx_dict[item] = chr(len(char_idx_dict))
    str1 = ''.join([char_idx_dict[item] for item in str1])
    str2 = ''.join([char_idx_dict[item] for item in str2])
    ops = Levenshtein.editops(str1, str2)
    return ops