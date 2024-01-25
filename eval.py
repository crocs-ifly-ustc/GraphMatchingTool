import os, sys
import pdb
import argparse
import Levenshtein
import tqdm
import shutil
import cv2
from multiprocessing import Process, synchronize, Lock, Manager, Pool
import multiprocessing
from six.moves import queue
import warnings
from collections import OrderedDict
import utils
from chemfig_struct import *
import ssml_parser
import chemfig_struct
import chemfig_ops
import graph_cmp
from viz_struct import *



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

def count_ops(ops):
    insert_nums = sum([1 for op_name, *_ in ops if op_name == 'delete'])
    substitute_nums = sum([1 for op_name, *_ in ops if op_name == 'replace'])
    delete_nums = sum([1 for op_name, *_ in ops if op_name == 'insert'])
    assert delete_nums + substitute_nums + insert_nums == len(ops)
    return delete_nums, substitute_nums, insert_nums

def compare_struct(lab_rep_dict, rec_rep_dict):
    ret = 0
    extra_cnt = 0
    lab_keys = set(list(lab_rep_dict.keys()))
    lab_cnt = len(lab_keys)

    all_keys = []
    all_keys += list(lab_rep_dict.keys())
    all_keys += list(rec_rep_dict.keys())
    all_keys = list(set(all_keys))

    for key in all_keys:
        if key in rec_rep_dict and key not in lab_rep_dict:
            extra_cnt += 1
        if key not in rec_rep_dict:
            rec_rep_dict[key] = {"text": "\\smear", "root": None}
        if key not in lab_rep_dict:
            lab_rep_dict[key] = {"text": "\\smear", "root": None}


    correct_cnt = 0
    for key in lab_rep_dict:
        if key not in rec_rep_dict:
            ret = 1
            continue
        lab_text = lab_rep_dict[key]["text"]
        lab_atom = lab_rep_dict[key]["root"]
        rec_text = rec_rep_dict[key]["text"]
        rec_atom = rec_rep_dict[key]["root"]

        rec_rep_dict[key]["lab_atom"] = lab_atom
        rec_rep_dict[key]["rec_atom"] = rec_atom
        rec_rep_dict[key]["lab_text"] = lab_text
        rec_rep_dict[key]["rec_text"] = rec_text

        if lab_atom and rec_atom:
            res_graph = graph_cmp.compare_graph(lab_atom, rec_atom)
            rec_rep_dict[key]["res_graph"] = res_graph
            rec_rep_dict[key]["res_cmp"] = res_graph
            if res_graph == 1:
                ret = 1
            else:
                if key in lab_keys:
                    correct_cnt += 1
        else:
            res_text = 0 if lab_text.replace(" ", "") == rec_text.replace(" ", "") else 1
            rec_rep_dict[key]["res_text"] = res_text
            rec_rep_dict[key]["res_cmp"] = res_text
            if res_text == 1:
                ret = 1
            else:
                if key in lab_keys:
                    correct_cnt += 1
    return ret, lab_cnt, correct_cnt, extra_cnt

def parse_rep_dict(rep_dict: OrderedDict, allow_exception=True):
    for key, value in rep_dict.items():
        begin = value.find("{")
        end = value.rfind("}")
        if begin == -1 or end == -1:
            rep_dict[key] = {"text": "\\chemfig { " + value + " }", "root": None}
            raise ValueError("can not find bracket pair in chemfig domain, {}".format(vaue))
        else:
            in_value = value[begin + 1:end].replace("\\\\", "").replace("\\enter", "").replace("\\space", " ")
            try:
                root_atom = ssml_parser.parse_ssml(in_value)
            except BaseException as e:
                if allow_exception:
                    print(e)
                    root_atom = None
                else:
                    raise e
            rep_dict[key] = {"text": "\\chemfig { " + in_value + " }", "root": root_atom}

def remove_fake_chemfig(rep_dict: OrderedDict, global_text=""):
    new_items = []
    new_global_text = global_text
    for key, value in rep_dict.items():
        text = value["text"]
        atom = value["root"]
        if atom is not None:
            chemfig_ops.SimulateCoord(atom)
            ret = chemfig_ops.NormFlatChemfig(atom)
        else:
            ret = None
        if isinstance(ret, Atom) or ret is None:
            new_items.append((key, value))
        else:
            new_text = " ".join(ret)
            new_global_text = new_global_text.replace(key, new_text)
    ind = 0
    key_map = {}
    new_rep_dict = OrderedDict()

    new_text_chemfig = []
    for key, value in new_items:
        new_key = "\\chem{}".format(chr(ord("a") + ind))
        ind += 1
        key_map[key] = new_key
        new_rep_dict[new_key] = value
        new_text_chemfig.append(value["text"])
    # pdb.set_trace()
    new_text_chemfig = " ".join(new_text_chemfig) 

    text_arr = new_global_text.split(" ")
    for i, text in enumerate(text_arr):
        if text in key_map:
            text_arr[i] = key_map[text]
    new_global_text = " ".join(text_arr)

    # new_text_chemfig = " ".join( list(new_rep_dict.values()) )
    return new_text_chemfig, new_rep_dict, new_global_text

def process_chemfig_str(inputStr, allow_exception = False):
    chem_trans, rep_dict, rep_trans = utils.replace_chemfig(inputStr)
    parse_rep_dict(rep_dict, allow_exception)
    chem_trans, rep_dict, rep_trans = remove_fake_chemfig(rep_dict, rep_trans)
    return chem_trans, rep_dict, rep_trans

def do_single_task(result, result_id, args, m_metrics, m_metrics_lock, 
            records_queue, records_queue_lock, 
            m_result_lines, m_result_lines_lock, 
            m_exception_lines, m_exception_lines_lock
    ):
    # reset id
    chemfig_struct.Atom.index = 0
    chemfig_struct.Bond.index = 0


    lab_line = result["lab"]
    lab = lab_line.split("\t")[1].strip()

    rec_line = result["rec"]
    rec = rec_line.split("\t")[1].strip()

    lab_arr = list(filter(None, lab.split(" ")))
    rec_arr = list(filter(None, rec.split(" ")))

    # base string compare
    base_ops = cal_edit_ops(rec_arr, lab_arr)

    # struct compare
    rec_chemfig, rec_rep_dict, rec_global_text = process_chemfig_str(rec, allow_exception=True)
    lab_chemfig, lab_rep_dict, lab_global_text  = process_chemfig_str(lab, allow_exception=False)

    all_ops = compare_struct(lab_rep_dict, rec_rep_dict)

    cmp_res, single_n, single_acc, single_extra = all_ops
    cmp_rec_global_text = rec_global_text
    cmp_lab_global_text = lab_global_text

    
    # update metrics
    ref_cmp_res = None
    with m_metrics_lock:
        #update base
        cur_metric = m_metrics["base"]
        cur_metric["sent_n"] += 1
        if len(base_ops) == 0:
            cur_metric["sent_correct"] += 1
        _d, _s, _i = count_ops(base_ops)
        cur_metric["d"] += _d
        cur_metric["s"] += _s
        cur_metric["i"] += _i
        cur_metric["n"] += len(lab_arr)
        m_metrics["base"] = cur_metric
        if args.ref_metric == "base":
            ref_cmp_res = 0 if len(base_ops) == 0 else 1

        #update struct
        cur_metric = m_metrics["struct"]
        cur_metric["sent_n"] += 1
        if cmp_res == 0:
            cur_metric["sent_correct"] += 1
        _d, _s, _i = (0, 0, 0)
        cur_metric["d"] += _d
        cur_metric["s"] += _s
        cur_metric["i"] += _i
        cur_metric["n"] += 1
        m_metrics["struct"] = cur_metric
        if args.ref_metric == "struct":
            ref_cmp_res = cmp_res

        # update struct line
        cur_metric = m_metrics["struct.line"]
        cur_metric["sent_n"] += 1
        cmp_line = 1
        if cmp_res == 0 and cmp_rec_global_text.replace(" ", "") == cmp_lab_global_text.replace(" ", ""):
            cur_metric["sent_correct"] += 1
            cmp_line = 0
        cur_metric["n"] += 1
        m_metrics["struct.line"] = cur_metric
        if args.ref_metric == "struct.line":
            ref_cmp_res = cmp_line



        

    if m_result_lines_lock is not None:
        # img_path, struct, struct.line, lab, rec
        cur_out_line = "{}\t{}\t{}\t{}\t{}\n".format(result["key"], cmp_res, cmp_line, lab, rec)
        with m_result_lines_lock:
            m_result_lines.put(cur_out_line)

    # do viz struct
    tgt_rec_rep_dict = rec_rep_dict
    if args.viz > 0:
        ori_trans = None
        if len(tgt_rec_rep_dict) == 0:
            tgt_rec_rep_dict = {"\\chema": {"res_text": cmp_res, "res_cmp": cmp_res, "lab_text": "\\smear", "rec_text": "\\smear"}}
        img_path = os.path.join(args.img_prefix, result["key"]+".jpg")
        viz_img = viz_struct_res(img_path, tgt_rec_rep_dict, ori_trans)
        img_name = os.path.splitext(os.path.basename(img_path))[0]

        output_path = os.path.join(args.viz_output, "res[{}]/{}.jpg".format(cmp_res, img_name))
        outDir = os.path.dirname(output_path)
        if not os.path.exists(outDir):
            os.makedirs(outDir, exist_ok=True)
        cv2.imwrite(output_path, viz_img)
    return

def try_do_single_task(result, result_id, args, m_metrics, m_metrics_lock, 
                        records_queue, records_queue_lock, 
                        m_result_lines, m_result_lines_lock, 
                        m_exception_lines, m_exception_lines_lock
                        ):
    try:
        do_single_task(result, result_id, args, m_metrics, m_metrics_lock, 
            records_queue, records_queue_lock, 
            m_result_lines, m_result_lines_lock, 
            m_exception_lines, m_exception_lines_lock
        )
    except BaseException as e:
        err_string = "[exception] result_id={}\terror={}\tline={}".format(result_id, e, str(result))
        print(err_string)
        with m_exception_lines_lock:
            m_exception_lines.put(err_string+"\n")

    if records_queue_lock is not None:
        with records_queue_lock:
            records_queue.put(1)

def main(args):
    rec_dict = {}
    with open(args.rec, "r", encoding="utf-8") as fin:
        lines = fin.readlines()
    for line in lines:
        key = line.split("\t")[0]
        key = os.path.splitext(os.path.basename(key))[0]
        rec_dict[key] = line
    
    lab_dict = {}
    with open(args.lab, "r", encoding="utf-8") as fin:
        lines = fin.readlines()
    for line in lines:
        key = line.split("\t")[0]
        key = os.path.splitext(os.path.basename(key))[0]
        lab_dict[key] = line

    results = []
    for key in lab_dict:
        result = {}
        result["key"] = key
        result["lab"] = lab_dict[key]
        result["rec"] = rec_dict[key] if key in rec_dict else ""
        results.append(result)

    # init output viz dir
    viz_dir = args.viz_output
    if viz_dir is not None:
        if not os.path.exists(viz_dir):
            if args.viz > 0:
                os.makedirs(viz_dir, exist_ok=True)

    # init metrics
    manager = Manager()
    base_item = {"sent_n": 0, "sent_correct": 0, "d": 0, "s": 0, "i": 0, "n": 0}
    base_item = manager.dict(base_item)
    metrics = {"base": base_item.copy()}
    metrics["struct"] = base_item.copy()
    metrics["struct.line"] = base_item.copy()

    #gen tasks
    records_queue = manager.Queue()
    records_queue_lock = manager.Lock()
    m_metrics = manager.dict(metrics)
    m_metrics_lock = manager.Lock()

    m_exception_lines = manager.Queue()
    m_exception_lines_lock = manager.Lock()

    m_result_lines = manager.Queue()
    m_result_lines_lock = manager.Lock()

    all_tasks = []
    result_id = -1
    if args.num_workers <= 0:
        for result in tqdm.tqdm(results):
            result_id += 1
            cur_task = (
                result, result_id, args, m_metrics, m_metrics_lock, records_queue, records_queue_lock, m_result_lines, m_result_lines_lock, m_exception_lines, m_exception_lines_lock,
            )
            if result["key"].find("valid_00824") == -1:
                continue
            do_single_task(*cur_task)
    else:
        for result in tqdm.tqdm(results):
            result_id += 1
            cur_task = (
                result, result_id, args, m_metrics, m_metrics_lock, records_queue, records_queue_lock, m_result_lines, m_result_lines_lock, m_exception_lines, m_exception_lines_lock,
            )
            all_tasks.append(cur_task)
            pass

        def print_error(error):
            print("error:", error)

        poolSize = args.num_workers
        pool = Pool(poolSize)
        pool.starmap_async(try_do_single_task, all_tasks, error_callback=print_error)
        pool.close()
        tq = tqdm.tqdm(total=len(all_tasks))
        count = 0
        print("begin")
        #try:
        while count < len(all_tasks):
            try:
                c = records_queue.get_nowait()
            except queue.Empty:
                continue
            count += 1
            tq.update(1)
        pool.join()

    # print metric
    for metric_name in m_metrics:
        metric = m_metrics[metric_name]
        sent_correct = metric["sent_correct"]
        sent_n = metric["sent_n"]
        d = metric["d"]
        s = metric["s"]
        i = metric["i"]
        n = metric["n"]
        if sent_n == 0:
            sent_n = 1e-10
        if n == 0:
            n = 1e-10
        sent_acc = float(sent_correct) / (sent_n)
        word_acc = float(n - d - s - i) / (n)
        word_cor = float(n - d - s) / (n)
        print("------ metric {} ------".format(metric_name))
        print("wacc={:.4f} % wcor={:.4f} % d={} s={} i={} n={}".format(word_acc * 100.0, word_cor * 100.0, d, s, i, n))
        print("sent acc = {:.4f} %( {}/{} )".format(sent_acc * 100.0, sent_correct, sent_n))

    if not m_exception_lines.empty():
        prefix, ext = os.path.splitext(args.rec)
        exception_path = "{}.exception.{}{}".format(prefix, args.ref_metric, ext)
        try:
            with open(exception_path, "w") as fout:
                while not m_exception_lines.empty():
                    line = m_exception_lines.get()
                    fout.write(line)
        except BaseException as e:
            pass
    
    with m_result_lines_lock:
        with open(args.output, "w") as fout:
            while not m_result_lines.empty():
                line = m_result_lines.get()
                fout.write(line)

    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser("graph matching tool")
    parser.add_argument("-rec", type=str, default="./test_decode_result.txt", help="output result by your model")
    parser.add_argument("-lab", type=str, default="./data/mini_valid/mini_valid_ssml_normed.txt", help="label file")
    parser.add_argument("-output", type=str, default="./result.txt", help="output file")
    parser.add_argument("-viz", type=int, default=0, help="if bigger than zero, will do visualization")
    parser.add_argument("-img_prefix", type=str, default="data/mini_valid/mini_valid", help="img root path")
    parser.add_argument("-ref_metric", type=str, default="struct", help="when viz is on, will divide result into true or false according to this metric")
    parser.add_argument("-viz_output", type=str, default=None, help="visualization output dir")
    parser.add_argument("-num_workers", type=int, default=64, help="for multi process")
    args = parser.parse_args()
    main(args)
