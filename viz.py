# -*- coding: utf-8 -*- 
import sys
import argparse
import json
import pdb
import os
import sys
import time
import numpy as np
import logging
import random
import re
import cv2
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
import math
import pdb

font_path = "./simsun.ttf"

#img is opencv img, textLines
def PrintImageWithInfo(srcImg, textLines):
    _h,_w,_c = srcImg.shape
    textLineCnt = len(textLines)
    maxTextLength = 0
    outH = _h
    for i in range(textLineCnt):
        if maxTextLength < len(textLines[i]):
            maxTextLength = len(textLines[i])

    labSize = int(_w/maxTextLength)
    if labSize < 15:
        labSize = 15
    src_font = ImageFont.truetype(font_path, size=labSize)
    labWidth = labSize*maxTextLength
    
    margin = 1
    outH = _h+(labSize + margin)*textLineCnt
    outW = max([_w, labWidth])
    outImgPIL = Image.new("RGB", (outW, outH), (255, 255, 255))
    draw = ImageDraw.Draw(outImgPIL)
    realW = 0
    for j in range(textLineCnt):
        curLine = textLines[j]
        if type(curLine) is not str:
            curLine = curLine.decode("utf-8")
        lw, lh = src_font.getsize(curLine)
        if realW<lw:
            realW=lw
        draw.text((0, _h+(margin+labSize)*j+margin), curLine, (0, 0, 0), font=src_font)
    outImg=np.frombuffer(outImgPIL.tobytes(),dtype=np.uint8).reshape((outH,outW,-1)).copy()
    outW2 = max([_w, realW])
    outImg = outImg[:,0:outW2,:]
    outImg[0:_h, 0:_w, :] = srcImg
    return outImg

def PrintImageWithInfoColor(srcImg, textLines, colorConfig):
    _h,_w,_c = srcImg.shape
    textLineCnt = len(textLines)
    maxTextLength = 0
    outH = _h
    for i in range(textLineCnt):
        if maxTextLength < len(textLines[i]):
            maxTextLength = len(textLines[i])

    labSize = int(_w/maxTextLength)
    if labSize < 15:
        labSize = 15
    src_font = ImageFont.truetype('simsun.ttc', size=labSize)
    labWidth = labSize*maxTextLength
    
    margin = 1
    outH = _h+(labSize + margin)*textLineCnt
    outW = max([_w, labWidth])
    outImgPIL = Image.new("RGB", (outW, outH), (255, 255, 255))
    draw = ImageDraw.Draw(outImgPIL)
    realW = 0
    for j in range(textLineCnt):
        curLine = textLines[j]
        if type(curLine) is not unicode:
            curLine = curLine.decode("utf-8")
        lw, lh = src_font.getsize(curLine)
        if realW<lw:
            realW=lw
        if colorConfig[j] is None:
            draw.text((0, _h+(margin+labSize)*j+margin), curLine, (0, 0, 0), font=src_font)
        else:
            curX = 0
            for idx in range(len(curLine)):
                cw, ch = src_font.getsize(curLine[idx])
                #pdb.set_trace()
                draw.text((curX, _h+(margin+labSize)*j+margin), curLine[idx], colorConfig[j][idx], font=src_font)
                curX += cw
    outImg=np.frombuffer(outImgPIL.tobytes(),dtype=np.uint8).reshape((outH,outW,-1)).copy()
    outW2 = max([_w, realW])
    outImg = outImg[:,0:outW2,:]
    outImg[0:_h, 0:_w, :] = srcImg
    return outImg


# rend

def adjust_font_and_width(in_text, fixed_height, ref_width, min_font_size=15, line_space=2):
    # min_font_size = 15
    line_space = 2
    cur_text_width = ref_width
    cur_lab_font_size = 25
    if cur_lab_font_size < min_font_size:
        cur_lab_font_size = min_font_size
    font = ImageFont.truetype(font_path, cur_lab_font_size, encoding="utf-8")
    while True:
        success = True
        cur_x = 0
        cur_y = 0
        cur_h = 0
        for ch in in_text:
            w, h = font.getsize(ch)
            if cur_x + w >= cur_text_width:
                cur_x = 0
                cur_y += (cur_h + line_space)
                cur_h = 0
            cur_x += w
            cur_h = max(cur_h, h)
            if cur_y + cur_h >= fixed_height:
                success = False
                break
        if success is True:
            break
        else:
            if cur_lab_font_size > min_font_size:
                cur_lab_font_size -= 1
                font = ImageFont.truetype(font_path, cur_lab_font_size, encoding="utf-8")
            else:
                cur_text_width += 100
    return cur_lab_font_size, cur_text_width

def drawText(text, height, width, font_size, line_space=2, align = "center"):
    font = ImageFont.truetype(font_path, font_size, encoding="utf-8")
    max_w = 0
    max_h = 0
    cur_x = 0
    cur_y = 0
    cur_h = 0
    img_np = np.ones((height, width, 3), dtype=np.uint8)*255
    #img = Image.new('RGB', (width, height), "#FFFFFF")
    img = Image.fromarray(img_np)
    imgDraw = ImageDraw.Draw(img)
    for ch in text:
        w, h = font.getsize(ch)
        if cur_x + w >= width:
            cur_x = 0
            cur_y += (cur_h + line_space)
            cur_h = 0
        imgDraw.text((cur_x, cur_y), ch, fill=(0,0,0), font=font)
        cur_x += w
        cur_h = max(cur_h, h)
        max_w = max([max_w, cur_x])
        max_h = max([max_h, cur_y + cur_h])
    # pdb.set_trace()
    img_np = np.array(img)
    img_out = np.ones((height, width, 3), dtype=np.uint8)*255
    #img_np[:max_h, :max_w, :]
    if align == "center":
        start_x = (width - max_w) // 2
        start_y = (height - max_h) // 2
    elif align == "left":
        start_x = 0
        start_y = 0
    else:
        raise NotImplementedError("unsupport align = {}".format(align))
    img_out[start_y:start_y+max_h, start_x:start_x+max_w, :] = img_np[:max_h, :max_w, :]
    return img_out

def drawTextOnImage(img_np, text, x, y, font_size, color):
    font = ImageFont.truetype(font_path, font_size, encoding="utf-8")
    img = Image.fromarray(img_np)
    imgDraw = ImageDraw.Draw(img)
    imgDraw.text((x, y), text, fill=color, font=font)
    out_img_np = np.array(img)
    return out_img_np
    
def adjust_font_and_height(in_text, fixed_width, ref_height, min_font_size=15, line_space=2):
    # min_font_size = 15
    line_space = 2
    cur_text_height = ref_height
    cur_lab_font_size = 25
    if cur_lab_font_size < min_font_size:
        cur_lab_font_size = min_font_size
    font = ImageFont.truetype(font_path, cur_lab_font_size, encoding="utf-8")
    while True:
        success = True
        cur_x = 0
        cur_y = 0
        cur_h = 0
        max_height = 0
        for ch in in_text:
            w, h = font.getsize(ch)
            if cur_x + w >= fixed_width:
                cur_x = 0
                cur_y += (cur_h + line_space)
                cur_h = 0
            cur_x += w
            cur_h = max(cur_h, h)
            if cur_y + cur_h >= cur_text_height:
                success = False
                break
            else:
                max_height = max([max_height, cur_y + cur_h])
        if success is True:
            break
        else:
            if cur_lab_font_size > min_font_size:
                cur_lab_font_size -= 1
                font = ImageFont.truetype(font_path, cur_lab_font_size, encoding="utf-8")
            else:
                cur_text_height += (cur_lab_font_size + line_space)
    # pdb.set_trace()
    return cur_lab_font_size, max_height

def vstack_images(imgs, margin=10):
    outH = 0
    outW = 0
    outC = 3
    for img in imgs:
        _h, _w, _c = img.shape
        outH += _h + margin
        outW = max([outW, _w])
    outH -= margin

    outImg = np.ones((outH, outW, outC), dtype=np.uint8)*255

    curY = 0
    for img in imgs:
        _h, _w, _c = img.shape
        outImg[curY:curY+_h, :_w, :] = img
        curY += _h + margin
    return outImg

def hstack_images(imgs, margin=10):
    outH = 0
    outW = 0
    outC = 3
    for img in imgs:
        _h, _w, _c = img.shape
        #outH += _h + margin
        outH = max([outH, _h])
        # outW = max([outW, _w])
        outW += _w + margin
    outW -= margin

    outImg = np.ones((outH, outW, outC), dtype=np.uint8)*255
    curX = 0
    for img in imgs:
        _h, _w, _c = img.shape
        #outImg[curY:curY+_h, :_w, :] = img
        outImg[:_h, curX:curX+_w, :] = img
        curX += _w + margin
    return outImg