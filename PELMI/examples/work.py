from re import search
from PIL import Image
import numpy as np
import generalizedreedsolo
from numpy import fromfile, array, uint8
import random
from iamge_matric import image_quantity_test
from inhance_image_all import inhance_image



def error_prone(sequence):
    temp = ""
    for i in range(len(sequence)):
        temp += sequence[i]
    for i in range(3):
        str1000 = temp[:10000]
        errorSum = str1000.count('ATA') + str1000.count('TAT') + str1000.count('GAC') + str1000.count(
        'CAC') + str1000.count('GTC') + str1000.count('GTG') + str1000.count('GCG') + str1000.count(
        'AAA') + str1000.count('TTT') + str1000.count('CCC') + str1000.count('CGC') + str1000.count('GGG')
        sum = errorSum/10000
        print(sum,"modify 的含量")

# 图像转为二进制数据    输出：二进制字符串  str
def image_to_bit(img1):
    img = Image.open(img1)
    red_channel, green_channel, blue_channel = img.split()
    # 转换为Numpy数组
    img_np = np.array(red_channel)
    img_np1 = np.array(green_channel)
    img_np2 = np.array(blue_channel)
    # 转换为一维数组
    quantized_coef_1d = img_np.flatten()
    quantized_coef_1d1 = img_np1.flatten()
    quantized_coef_1d2 = img_np2.flatten()
    # 一维数组转换为二进制数据
    compressed_data = ''.join([format(int(x), '08b') for x in quantized_coef_1d])
    compressed_data1 = ''.join([format(int(x), '08b') for x in quantized_coef_1d1])
    compressed_data2 = ''.join([format(int(x), '08b') for x in quantized_coef_1d2])
    return compressed_data, compressed_data1, compressed_data2

def bits_to_list(binary_string, segment_length):
    matrix = [int(bit) for bit in binary_string]  # 将字符串转为list数据类型
    if len(matrix) % segment_length != 0:
        matrix += [0] * (segment_length - len(matrix) % segment_length)
    # 将list转为二维list
    matrix = array(matrix)
    matrix = matrix.reshape(int(len(matrix) / segment_length), segment_length)
    return matrix.tolist()

# 添加索引部分
def connect(index, bit_segment, index_binary_length):
    bin_index = list(map(int, list(str(bin(index))[2:].zfill(index_binary_length))))
    one_list = bin_index + bit_segment
    return one_list
def connect_all(bit_segments, index_binary_length):
    connected_bit_segments = []
    for row in range(len(bit_segments)):
        connected_bit_segments.append(connect(row, bit_segments[row], index_binary_length))
    return connected_bit_segments

# 数据分割   返回的数据类型 ['010101', '010101']
def bits_to_divide(bit_segments):
    bit_segments_divide_all = []
    bit_segments_divide1_all = []
    for i in range(len(bit_segments)):
        bit_segments_divide = []
        bit_segments_divide1 = []
        data = bit_segments[i]
        # 分块
        blocks = [data[j:j + 8] for j in range(0, len(data), 8)]
        # 将每块的前四位和后四位分别合成一个列表
        result = [[block[0:4], block[4:]] for block in blocks]
        for j in range(len(result)):
            bit_segments_divide.extend(result[j][0])
            bit_segments_divide1.extend(result[j][1])

        temp_bit_divide =  ''.join(map(str, bit_segments_divide))
        temp_bit_divide1 = ''.join(map(str, bit_segments_divide1))

        bit_segments_divide_all.append(temp_bit_divide)
        bit_segments_divide1_all.append(temp_bit_divide1)

    return bit_segments_divide_all, bit_segments_divide1_all


def binary_to_dna(binary_data):
    dna_data = ''
    mapping = {
        '00': 'A',
        '01': 'T',
        '10': 'C',
        '11': 'G'
    }
    for i in range(0,len(binary_data), 2):
        dna_data += mapping[binary_data[i: i +2]]
    return dna_data

def binary_to_dna_right(binary_data):
    dna_right_all = []
    mapping = {
        '00': 'A',
        '01': 'T',
        '10': 'C',
        '11': 'G'
    }
    mapping1 = {
        '00': 'T',
        '01': 'A',
        '10': 'G',
        '11': 'C'
    }
    mapping2 = {
        '00': 'G',
        '01': 'C',
        '10': 'A',
        '11': 'T'
    }
    mapping3 = {
        '00': 'C',
        '01': 'G',
        '10': 'T',
        '11': 'A'
    }

    dna_data = ""
    for i in range(0,len(binary_data), 2):
        dna_data += mapping[binary_data[i: i +2]]
    dna_right_all.append(dna_data)
    dna_data = ""
    for i in range(0,len(binary_data), 2):
        dna_data += mapping1[binary_data[i: i +2]]
    dna_right_all.append(dna_data)
    dna_data = ""
    for i in range(0,len(binary_data), 2):
        dna_data += mapping2[binary_data[i: i +2]]
    dna_right_all.append(dna_data)
    dna_data = ""
    for i in range(0,len(binary_data), 2):
        dna_data += mapping3[binary_data[i: i +2]]
    dna_right_all.append(dna_data)
    return dna_right_all

# DNA数据块进行重组
def dna_block_res(dna_left, dna_right_all):
    temp_blovk_res = []
    for i in dna_right_all:
        temp = dna_left[:2] + i[:2] + dna_left[2:] + i[2:]
        temp_blovk_res.append(temp)
    return temp_blovk_res

def homopolymer(sequence, max_homopolymer):
    homopolymers = "A{%d,}|C{%d,}|G{%d,}|T{%d,}" % tuple([1 + max_homopolymer] * 4)
    return False if search(homopolymers, sequence) else True

def gc_content(sequence, max_content):
    return (1 - max_content) <= (float(sequence.count("C") + sequence.count("G")) / float(len(sequence))) <= max_content

def gc_content_count(sequence, max_content):
    content = float(sequence.count("C") + sequence.count("G")) / float(len(sequence))
    return abs(0.5 - content)

def modify_content(sequence_str):
    sum = 0
    patterns = ['TAT', 'CGC', 'GTC', 'GTG', 'GAC', 'CAC', 'GCG', 'AGA', 'ATA', 'AAA', 'CCC', 'TTT', 'GGG']
    for pattern in patterns:
        count = sequence_str.count(pattern)
        sum += count
    return sum

# 返回编码数据块和  规则（0---3）
def check(data_left, data_right, dna_divide):
    dna_left = binary_to_dna(data_left)     # 这是个字符串的数据类型
    dna_right_all = binary_to_dna_right(data_right)     # 这是个列表的数据类型
    temp_blovk_res_all = dna_block_res(dna_left, dna_right_all)  # 将数据块进行重组
    # 约束的逻辑，首先是 ----数据块的均聚物小于4
    # 第二：GC含量尽可能的靠近  0.52
    # 不期望基序尽可能的少
    max_homopolymer = 4
    max_content = 0.57
    rule = []
    rule_satgchomo = []
    temp_filter = []
    # 将满足均聚物约束的筛选出来
    for i in range(len(temp_blovk_res_all)):
        temp_block = temp_blovk_res_all[i]
        if max_homopolymer and not homopolymer(temp_block, max_homopolymer):
            continue
        rule.append(i)
        temp_filter.append(temp_block)
    # 将整条链进行整合
    dna_divide_all = []
    dna_divide_all_satgchomo = []   # 替身的作用
    for i in range(len(temp_filter)):
        dna_divide_temp = dna_divide
        dna_divide_temp += temp_filter[i]
        dna_divide_all.append(dna_divide_temp)
    # 对其GC含量进行约束
    gc_sat_sum = 0
    gc_count = []
    for i in range(len(dna_divide_all)):
        for j in range(len(dna_divide_all)):
             if gc_content(dna_divide_all[j], max_content):
                 gc_sat_sum += 1
        if gc_sat_sum == 0:
            for t in range(len(dna_divide_all)):
                temp_gc_content = gc_content_count(dna_divide_all[t], max_content)
                gc_count.append(temp_gc_content)
            min_index = gc_count.index(min(gc_count))
            return dna_divide_all[min_index][-8:], rule[min_index]
        else:
            if gc_content(dna_divide_all[i], max_content):
                dna_divide_all_satgchomo.append(dna_divide_all[i])
                rule_satgchomo.append(rule[i])

    # 获得满足GC 和均聚物的序列
    modify_count = []
    for i in range(len(dna_divide_all_satgchomo)):
        modify_count.append(modify_content(dna_divide_all_satgchomo[i]))
    min_index = modify_count.index(min(modify_count))    # 获得不期望基序最小的序列索引
    # return dna_divide_all[min(gc_count)], rule_satgchomo[min_index]
    return dna_divide_all_satgchomo[min_index][-8:], rule_satgchomo[min_index]

# 四进制转二进制是   输入：[[0, 1,3,2],[0, 1,3,2]], 输出 【"00011110"，"00011110"】
def four_to_two(rule_all):
    rule_list_all = []
    for rule in rule_all:
        matrix = ""
        for ru in rule:
            matrix += str(bin(ru))[2:].zfill(2)
        rule_list_all.append(matrix)
    return rule_list_all


def encode_divide_func(bit_segments_divide_all_left, bit_segments_divide1_all_right):
    dna_segments = []
    rule_all = []

    for i in range(len(bit_segments_divide_all_left)):
        # bit----to---dna
        rule = []
        dna_divide = ""
        for j in range(0, len(bit_segments_divide_all_left[i]), 8):
            # binary_to_dna(bit_segments_divide_all[i][j: j+8])
            dna_block, rule_temp = check(bit_segments_divide_all_left[i][j: j + 8],
                                         bit_segments_divide1_all_right[i][j: j + 8], dna_divide)
            dna_divide += dna_block
            rule.append(rule_temp)
        dna_segments.append(dna_divide)
        rule_all.append(rule)
    return dna_segments, rule_all

#s是字符串数据类型，取字符串的前32位，后面数据分割提取前四位
def process_string(s):
    chunk = s
    # 提取前32个字符
    prefix = chunk[:16]
    # 提取32位之后的数据
    suffix = chunk[16: -40]
    # 将32位之后的数据分为八位一组，取前四位
    processed_suffix = ''.join(suffix[i:i + 8][:4] for i in range(0, len(suffix), 8))
    result = prefix + processed_suffix + chunk[-40:]
    return result

#s是字符串数据类型(解码)，取字符串的前32位，后面数据分割提取前四位
def process_string_decode(s):
    chunk = s
    # 提取前32个字符
    prefix = chunk[:16]
    # 提取32位之后的数据
    suffix = chunk[16: -56]
    # 将32位之后的数据分为八位一组，取前四位
    processed_suffix = ''.join(suffix[i:i + 8][:4] for i in range(0, len(suffix), 8))
    result = prefix + processed_suffix + chunk[-56:]
    return result

#  data字符串数据类型
def divide_temp(data):
    temp_all = []
    for i in range(len(data)):
        temp = data[i]
        temp_group = process_string(temp)
        temp_all.append(temp_group)
    return temp_all

def cut(obj, sec):  # 按一定长度划分
    return [obj[i:i + sec] for i in range(0, len(obj), sec)]

def two_baseconversion(bitlist):  # 2进制转10进制
    tenbaseconversion = []
    for str_1 in bitlist:
        tenbaseconversion.append(int(str_1, 2))
    return tenbaseconversion

# 将二进制数据转为像素值形式   输出[[23, 25, 234],[23, 25, 234]]
def bit_to_demical_list(tempbinarymsg):
    demical_all = []
    for i in range(len(tempbinarymsg)):
        binarylen = 8
        temp8bitmsg1 = cut(tempbinarymsg[i], binarylen)
        orimsg = two_baseconversion(temp8bitmsg1)
        demical_all.append(orimsg)
    return demical_all

# 将二进制数据转为像素值形式   输出[[23, 25, 234],[23, 25, 234]]
def bit_to_demical_list1(tempbinarymsg):
    binarylen = 8
    temp8bitmsg1 = cut(tempbinarymsg, binarylen)
    orimsg = two_baseconversion(temp8bitmsg1)
    return orimsg

# 十进制转为二进行并进行补充
def to_bin(value, num):  # 十进制数据转8位2进制
    bin_chars = ""
    temp = value
    for i in range(num):
        bin_char = bin(temp % 2)[-1]
        temp = temp // 2
        bin_chars = bin_char + bin_chars
    return bin_chars.upper()  # 输出指定位宽的二进制字符串

def reed_solomon_model(symbol_size, message_length):
    return generalizedreedsolo.Generalized_Reed_Solomon(symbol_size,  # 伽罗华域的指数
                                                        message_length,  # 加入校验位后的长度
                                                        field_size=2,  # 伽罗华域的底数
                                                        payload_length=23,  # 原始消息长度
                                                        p_factor=1)  # 加速因子

def rs(data):
    # rs分割后的数据【'010010010101010'，'010010010101010'】
    binarylen = 8
    divide_binary_data = divide_temp(data)
    orimsg_divide = bit_to_demical_list(divide_binary_data)  # 将二进制数据转为像素值形式-----分割后的数据
    test_normal_rs = reed_solomon_model(binarylen, 25)  # 生成伽罗华域模型
    #   加入纠错吗的部分
    normal_msg = []
    k = 0
    for templist in orimsg_divide:
        if len(templist) == 23:
            normal_msg.append(test_normal_rs.encode(templist))
        else:
            supplementarybit = '0' * (23 - len(templist))
            templist.extend(map(int, list(supplementarybit)))
            normal_msg.append(test_normal_rs.encode(templist))

    # 获取纠错吗  并将其转为二进制数据类型
    rs_bit_list = []
    for i in range(len(normal_msg)):
        temp = normal_msg[i][-2 :]
        split = ""
        for j in range(len(temp)):
            value_1 = to_bin(temp[j], binarylen)
            split += str(value_1)
        rs_bit_list.append(split)
    return rs_bit_list


def encode(bit_segments):
    # 数据类型：['010101', '010101']
    bit_segments_divide_all_left, bit_segments_divide1_all_right = bits_to_divide(bit_segments)
    # 前半部分数据编码   # 后半部分数据编码
    # 将索引部分和有效数据部分进行编码
    dna_index_data_seg, rule_index_data = encode_divide_func(bit_segments_divide_all_left, bit_segments_divide1_all_right)

    a1 = dna_index_data_seg[-1]
    a2 = rule_index_data[-1]
    # step1： 将规则位 转为二进制数据类型      输出数据类型：【"00011110"，"00011110"】
    rule_bits_index_data = four_to_two(rule_index_data)
    a3 = rule_bits_index_data[-1]
    # 对规则位取其前32位二进制数据
    rule_bits_index_data_32 = []
    for i in rule_bits_index_data:
        rule_bits_index_data_32.append(i[:32])

    # 将规则位分为 前部分数据和后部分数据
    bit_rule_divide_all_left, bit_rule_divide1_all_right = bits_to_divide(rule_bits_index_data_32)
    # 将rule部分进行编码    rule_rule规则位编码产生的 四进制规则
    dna_rule_seg, rule_rule = encode_divide_func(bit_rule_divide_all_left, bit_rule_divide1_all_right)

    # 将规则位的规则位  转为二进制数据   rule_bits_rule两位左右
    rule_bits_rule = four_to_two(rule_rule)

    # 合成index  data  rule(34位二进制数据)    rule_rule（4位二进制数据）    补充的两位二进制数据
    # 同上，合成DNA序列
    # 将bit_segments的数据类型由[[0, 1, 1, 0, 1, 1, 1, 0], [0, 1, 1, 0, 1, 1, 1, 0]]转为['01101110', '01101110']
    bit_segments = result = [''.join(map(str, sublist)) for sublist in bit_segments]
    bit_segments_nors = []
    dna_segments_nors = []
    for i in range(len(bit_segments)):

        temp_bit = bit_segments[i] + rule_bits_index_data[i] + rule_bits_rule[i] + "00"
        # 将数据和规则为的前三十二位数据  合并
        # 将索引  数据   规则位  和 规则位产生的规则位以及补充的一位    共同构成temp_dna
        rule_dna_rule = binary_to_dna(temp_bit[-8:])
        temp_dna = dna_index_data_seg[i] + dna_rule_seg[i] + rule_dna_rule
        # 数据类型均为['010010100', '0101010010']
        bit_segments_nors.append(temp_bit)
        dna_segments_nors.append(temp_dna)

    # 添加纠错码    返回的是rs纠错码转为二进制数据   【'00000000011111111'，'00000000011111111'】 每个元素16位
    rs_bits = rs(bit_segments_nors)
    # rs数据分割
    bit_rs_divide_all_left, bit_rs_divide1_all_right = bits_to_divide(rs_bits)
    # rs编码
    dna_rs_seg, rule_rs = encode_divide_func(bit_rs_divide_all_left, bit_rs_divide1_all_right)
    # 将rs规则位 转为二进制数据类型      输出数据类型：【"00011110"，"00011110"】
    rule_bits_rs = four_to_two(rule_rs)
    # 将rs规则位 的规则位转为  碱基形式
    dna_segments_all = []
    for i in range(len(rule_bits_rs)):
        rs_dna_rule = binary_to_dna(rule_bits_rs[i])
        temp_dna_all = dna_segments_nors[i] + dna_rs_seg[i] + rs_dna_rule
        dna_segments_all.append(temp_dna_all)

    #print(dna_segments_all[-1], "encode_dna")
    a4 = dna_segments_all[-1]
    # 添加标识位
    last_dna_segments_all = []
    for i in range(len(dna_segments_all)):
        #temp_last = dna_segments_all[i][:32] + "AC" + dna_segments_all[i][32:64] + "AC" + dna_segments_all[i][64:96] + "AC" + dna_segments_all[i][96:128] + "AC" +dna_segments_all[i][128:]

        temp_last = dna_segments_all[i][:16] + "A" + dna_segments_all[i][16:32] + "AC" + dna_segments_all[i][32:48] + "A" + dna_segments_all[i][48:64] + "AC" + dna_segments_all[i][
                                                                                          64:80] + "A" + dna_segments_all[i][80:96] + "AC" + \
                    dna_segments_all[i][96:112] + "A" + dna_segments_all[i][112:128] + "AC" + dna_segments_all[i][128:]

        last_dna_segments_all.append(temp_last)
    return last_dna_segments_all

# 去除单个标识位
def remove_1_flag(dna_segments):

    #1将数据分割为33个为一块   2 判断其中的第17位时候分A，然后将数据块进行重组（完成工作）
    remove_dna = []
    miss_index = []
    for i in range(len(dna_segments)):
        temp_block = ""
        temp_lian = ""
        t = 0
        for j in range(0,132,33):
            if dna_segments[i][j + 16] == "A":
                t += 1
            temp_block = dna_segments[i][j : j+ 16] + dna_segments[i][j + 17 : j + 33]
            temp_lian += temp_block
        temp_lian += dna_segments[i][-37 :]
        if t < 2:
            miss_index.append(i)
        remove_dna.append(temp_lian)
    return remove_dna, miss_index

# 去除标识位
def remove_flag(dna_segments):
    # step1 :
    remove_dna_segment = []

    test1_no = 0
    test1_no_list = []
    test1_cor= 0
    test1_cor_list = []
    test1_block = []
    test1_data = []
    miss_index_in_del = []
    for i in range(len(dna_segments)):
         dna_segment = ""
         t = 33
         for j in range(4):
            if dna_segments[i][t:t + 2] == "AC":
                dna_segment += dna_segments[i][t - 33:t]
            elif dna_segments[i][t - 1:t + 1] == "AC":
                dna_segment += dna_segments[i][t - 32:t]
                dna_segment += "A"
                test1_cor += 1
                test1_cor_list.append(i)
                miss_index_in_del.append(i)


            elif dna_segments[i][t + 1:t + 3] == "AC":
                dna_segment += dna_segments[i][t - 33:t]
                test1_cor += 1
                test1_cor_list.append(i)
                miss_index_in_del.append(i)

            elif dna_segments[i][t:t + 2] != "AC" and dna_segments[i][t - 1:t + 1] != "AC" and dna_segments[i][t +1:t + 3] != "AC":
                dna_segment += dna_segments[i][t - 33:t]
                test1_block.append(dna_segments[i][t - 33:t])
                test1_data.append(dna_segments[i])
                test1_no += 1
                test1_no_list.append(i)
                # 将错误索引储存下来
                temp2 = bin(i)[2:].zfill(16)
                temp2 = int(temp2, 2)
                #temp2 = [int(bit) for bit in temp2]
                miss_index_in_del.append(temp2)
            t += 35
         dna_segment += dna_segments[i][-37 :]
         remove_dna_segment.append(dna_segment)

    remove_dna_segment,miss_1_index = remove_1_flag(remove_dna_segment)
    miss_index_in_del.extend(miss_1_index)

    return  remove_dna_segment, miss_index_in_del



# 将数据的字符串DNA序列转为  二进制数据
def mapping_dna_to_bit(dna_data, flag):
    mapp = {
        'A':'00',
        'T':'01',
        'C':'10',
        'G':'11'
    }
    mapping = {
        'A':'00',
        'T':'01',
        'C':'10',
        'G':'11'
    }
    mapping1 = {
        'T': '00',
        'A': '01',
        'G': '10',
        'C': '11'
    }
    mapping2 = {
        'G': '00',
        'C': '01',
        'A': '10',
        'T': '11'
    }
    mapping3 = {
        'C': '00',
        'G': '01',
        'T': '10',
        'A': '11'
    }
    if flag == "A":
        mapp = mapping
    elif flag == "T":
        mapp = mapping1
    elif flag == "C":
        mapp = mapping2
    else:
        mapp = mapping3

    bits_block = ""
    for i in range(len(dna_data)):
        bits_block += mapp[dna_data[i]]
    return bits_block

def dna_to_divide(bit_segments, flag):
        bit_segments_all = []
        for i in range(len(bit_segments)):
            data = bit_segments[i]
            # 分块
            blocks = [data[j:j + 8] for j in range(0, len(data), 8)]
            # 将每块的前四位和后四位分别合成一个列表
            result = [[block[0:4], block[4:]] for block in blocks]
            temp_block = ""
            for j in range(len(result)):
                temp_block_left = ''.join(map(str, result[j][0][:2])) + ''.join(map(str, result[j][1][:2]))
                temp_block_right = ''.join(map(str, result[j][0][2:])) + ''.join(map(str, result[j][1][2:]))
                # temp_block_left 按照规则一进行解码， temp_block_right 按照特地规则进行解码  数据类型位"0011"
                temp_block_left = mapping_dna_to_bit(temp_block_left, "A")

                temp_block_right = mapping_dna_to_bit(temp_block_right, flag[i][j])
                temp_block = temp_block + temp_block_left[:4] + temp_block_right[:4] + temp_block_left[4:] + temp_block_right[4:]
            bit_segments_all.append(temp_block)
        return bit_segments_all


def dna_to_bits(dna_data):
    data = []
    mapping = {
        'A':'00',
        'T':'01',
        'C':'10',
        'G':'11'
    }
    for i in dna_data:
        data.append(mapping[i])
    return data


def finalfile_re_zuhe(origin_data, divide_data):
    temp_all = []
    for i in range(len(origin_data)):
        temp_origin_data = origin_data[i][16: -56]
        temp_divide_data = divide_data[i][16 : -40]
        temp = ""
        for t in range(0, len(temp_divide_data), 4):
            temp = temp + temp_divide_data[t:t + 4] + temp_origin_data[2*t + 4: 2*t + 8]
        temp = divide_data[i][:16] + temp
        temp_all.append(temp)
    return temp_all


def decode_rs(bits_segments_all):
    #  对二进制数据进行重新分组，finalfile_divide是['0100101010','010100101']数据类型
    finalfile_divide_left = []

    # 获取每条链中 rs重点保护的数据
    for i in bits_segments_all:
        temp_data = process_string_decode(i)
        finalfile_divide_left.append(temp_data)

    # 将其转为十进制数据
    finalfile_ten_left = []
    for i in finalfile_divide_left:
        orimsg_divide = bit_to_demical_list1(i)  # 将二进制数据转为像素值形式-----分割后的数据
        finalfile_ten_left.append(orimsg_divide)

    correct = []
    for dis_m in finalfile_ten_left:
        model = reed_solomon_model(8, len(dis_m))
        # noinspection PyBroadException
        try:
            correct.append(model.decode(dis_m))  # 去掉纠错吗的  像素值
        except Exception:
            correct.append(dis_m[:23])

    # 自己添加的部分   当数据中不是20时，保证列表是20
    for i in range(len(correct)):
        if len(correct[i]) < 23:
            for r in range(23 - len(correct[i])):
                correct[i].append(0)
        if len(correct[i]) > 23:
            correct[i] = correct[i][:23]

    # step 2：将像素值形式其转为二进制数据，并重新进行排序
    bit_temp = []
    for msg in correct:
        split = ""
        for value_1 in msg:
            value_1 = to_bin(value_1, 8)
            split += str(value_1)
        bit_temp.append(split)
    bit_data_all = finalfile_re_zuhe(bits_segments_all, bit_temp)  # 对数据重新进行组合
    return bit_data_all

# 对序列进行筛选
def shaixuan(data):
    binary_list = data
    dict_list = {}

    for sub_list in binary_list:
        index = sub_list[:16]
        key = "".join(str(bit) for bit in index)
        if key not in dict_list:
            dict_list[key] = sub_list

    filtered_list = list(dict_list.values())
    return filtered_list

# 对解码失败的链进行补充
def data_subment(nested_list):
    data_sub = []
    miss_index = []
    y = 0
    # step1: 排序
    sorted_list = sorted(nested_list, key=lambda sublist: int(''.join(str(x) for x in sublist[:16])), reverse=False)
    sorted_list = shaixuan(sorted_list)
    # step 2 :获取文件大小，确定链的长度
    length = len(nested_list)
    z = 0
    for i in range(length):
        t = 0
        decimal_value = int(''.join(map(str, sorted_list[z][:16])), 2)
        if decimal_value == i:
            data_sub.append(sorted_list[z])
            z += 1
        else:
            miss_index.append(i)
            # 找到这个数组，并重新赋值
            y = y + 1
            temp0 = ""
            if i == 0:
                for j in range(256):
                    temp0 += "0"
                data_sub.append(temp0)
            elif 0 < i < 8:
                t = i - 1
            else:
                t = i - 8
            temp1 = data_sub[t][16:]
            temp2 = bin(i)[2:].zfill(16)
            temp2 = [int(bit) for bit in temp2]
            temp2.extend(temp1)
            data_sub.append(temp2)
            # data_sub.append(temp0)
    print(z)
    return data_sub, miss_index

# 整体性均值迭代
def data_subment1(nested_list):
    data_sub = []
    miss_index = []
    y = 0
    # step1: 排序
    sorted_list = sorted(nested_list, key=lambda sublist: int(''.join(str(x) for x in sublist[:16])), reverse=False)
    sorted_list = shaixuan(sorted_list)

    # 1将数据索引高出的数据筛选出来（完成）
    # 2对未超出部分进行判断和补充，超出部分直接进行随机处理
    # step 2 :获取文件大小，确定链的长度
    length = len(nested_list)
    z = 0
    for i in range(length):
        decimal_value1 = int(''.join(map(str, sorted_list[-1][:16])), 2)
        if i <= decimal_value1:
            decimal_value = int(''.join(map(str, sorted_list[z][:16])), 2)
            if decimal_value == i:
                data_sub.append(sorted_list[z])
                z += 1
            else:
                miss_index.append(i)
                list = []

                length1 = len(sorted_list[0])
                for i in range(length1):
                    temp = random.choice([0, 1])
                    list.append(temp)
                data_sub.append(list)
        else:
            list = []
            length1 = len(sorted_list[0])
            for i in range(length1):
                temp = random.choice([0, 1])
                list.append(temp)
            data_sub.append(list)
    print(z)
    return data_sub, miss_index

def decode(dna_segments):
    # 去除标识位
    dna_segments, miss_index_in_del = remove_flag(dna_segments)
    # 数据转化
    dna_rs = []
    dna_rs_flag = []
    dna_flag_16nt = []
    dna_flag_flag_2nt = []
    dna_flag_flag_1nt = []
    dna_flag_flag_4nt = []
    dna_index_data = []
    dna_index_data_flag = []
    for i in range(len(dna_segments)):
        # 分别获取rs的纠错位和rs的标识位
        dna_rs.append(dna_segments[i][-9: -1])
        dna_rs_flag.append(dna_segments[i][-1 :])
        # 分别获取标识位和标识位的标识位
        dna_flag_16nt.append(dna_segments[i][-29 : -13])
        dna_flag_flag_4nt.append(dna_segments[i][-13:-9])
        dna_flag_flag_2nt.append(dna_segments[i][-12:-10])
        dna_flag_flag_1nt.append(dna_segments[i][-13:-12])
        # 分别获取index，data 和他们的标识位
        dna_index_data.append(dna_segments[i][:136])
    # 解码      ---rs
    bits_rs = dna_to_divide(dna_rs, dna_rs_flag)    # rs
    # 标识位
    bits_flag = dna_to_divide(dna_flag_16nt, dna_flag_flag_2nt)
    bits_flag_2 = dna_to_bits(dna_flag_flag_1nt)               # 将第17位数据转为二进制形式

    bits_flag_34 = []
    dna_flag_17_rule0 = []
    bits_flag_flag_8 = []
    # 获得完整的17位
    for i in range(len(bits_flag)):
        temp_bits = bits_flag[i] + bits_flag_2[i]
        bits_flag_34.append(temp_bits)
        temp_dna = binary_to_dna(temp_bits)
        dna_flag_17_rule0.append(temp_dna)
        temp_bits_8 = mapping_dna_to_bit(dna_flag_flag_4nt[i],"A")
        bits_flag_flag_8.append(temp_bits_8)

    # index,data
    bits_index_data = dna_to_divide(dna_index_data, dna_flag_17_rule0)

    # 将所有二进制数据进行拼接
    bits_segments_all = []
    for i in range(len(bits_index_data)):
        last_temp = bits_index_data[i] + bits_flag[i] + bits_flag_flag_8[i] + bits_rs[i]
        bits_segments_all.append(last_temp)

    #a3 = bits_segments_all[-1]
    #print(a3, "decode——rs纠错前")
    # 最后进行rs纠错，就完成工作了
    bits_segments_all = decode_rs(bits_segments_all)
    #print(bits_segments_all[-1], "decode——rs纠错后")
    #a4 = bits_segments_all[-1]
    #  进行修改，解码失败的链进行补充, 并获取丢失链的索引
    bits_segments_all, miss_index = data_subment1(bits_segments_all)

    # 将两部分索引进行合并，并完成排序
    miss_index.extend(miss_index_in_del)
    print("miss_index", miss_index)
    unique_list = list(set(miss_index))    #去除列表中重复的部分
    miss_index = sorted(unique_list)
    #a5 = bits_segments_all[-1]
    #print(bits_segments_all[-1], "decode——数据链补充")

    # 去掉索引位
    for i in range(len(bits_segments_all)):
        bits_segments_all[i] = bits_segments_all[i][16:]
    a6 = bits_segments_all[-1]
    return bits_segments_all, miss_index

# data_r 通道的数据，    data_window 窗口的数据
def Enhance_image(data_r, data_window, row, column):
    size = len(data_window[0]) + len(data_window[1]) + len(data_window[2])
    if size == 9:
        data_med = np.median(data_window)
        data_sort = np.sort(data_window.reshape(-1))
        A2 = data_sort[1]
        A8 = data_sort[7]
        if (A2 < data_r[row, column] and data_r[row, column] <= data_med) or (data_med <= data_r[row, column] and data_r[row, column] < A8):
            data_r[row, column] = data_r[row, column]
            #data_r[row, column] = 0
        else:
            # 测试
            #data_med = 0
            data_r[row, column] = data_med
    return data_r


# 将错误索引传参,获得遍历后的索引数据   data_r       error_index
def enhance_image_main(data_r, error_index):
    for t in range(len(error_index)):
        index = error_index[t]
        # set the window size
        row = int(index / 8)
        colunm = index % 8
        # slide without stack     data_r.shape[]  0和1分别代表元祖的行和列
        # 代码的缺陷部分是边框部分处理不了
        for i in range(row, row + 1, 1):
            for j in range(colunm*32, colunm*32 + 32, 1):
                if 0 < i < data_r.shape[0] - 1 and 0 < j < data_r.shape[1] - 1:
                    data_r = Enhance_image(data_r, data_r[i - 1:i + 2, j - 1:j + 2],
                                  i, j)
                else:
                    continue
    return data_r


def bit_to_image(img, compressed_data , compressed_data1 ,compressed_data2, error_index_all):
    img = Image.open(img)
    red_channel, green_channel, blue_channel = img.split()
    img_np = np.array(red_channel)
    error_index = error_index_all[0]
    error_index1 = error_index_all[1]
    error_index2 = error_index_all[2]
    # 对数据进行补充
    wigth, height = img.size
    sum = wigth*height*8
    if len(compressed_data) < sum:
        compressed_data = str(compressed_data).ljust(sum, '0')
    if len(compressed_data1) < sum:
        compressed_data1 = str(compressed_data1).ljust(sum, '0')
    if len(compressed_data2) < sum:
        compressed_data2 = str(compressed_data2).ljust(sum, '0')
    quantized_coef_1d = np.array([int(compressed_data[i:i + 8], 2) for i in range(0, len(compressed_data), 8)])
    quantized_coef_1d1 = np.array([int(compressed_data1[i:i + 8], 2) for i in range(0, len(compressed_data1), 8)])
    quantized_coef_1d2 = np.array([int(compressed_data2[i:i + 8], 2) for i in range(0, len(compressed_data2), 8)])
    # 将量化系数重构为二维数组
    quantized_coef = np.reshape(quantized_coef_1d, img_np.shape)
    quantized_coef = enhance_image_main(quantized_coef, error_index)
    quantized_coef = quantized_coef.astype(np.uint8)
    print(quantized_coef, "quantized_coef")
    quantized_coef1 = np.reshape(quantized_coef_1d1, img_np.shape)
    quantized_coef1 = enhance_image_main(quantized_coef1, error_index1)
    quantized_coef1 = quantized_coef1.astype(np.uint8)
    quantized_coef2 = np.reshape(quantized_coef_1d2, img_np.shape)
    quantized_coef2 = enhance_image_main(quantized_coef2, error_index2)
    quantized_coef2 = quantized_coef2.astype(np.uint8)
    # 将numpy 转为通道数据
    reconstructed_red_channel = Image.fromarray(quantized_coef)
    reconstructed_red_channel1 = Image.fromarray(quantized_coef1)
    reconstructed_red_channel2 = Image.fromarray(quantized_coef2)
    # 数据合并
    reconstructed_image = Image.merge("RGB", (reconstructed_red_channel, reconstructed_red_channel1, reconstructed_red_channel2))
    # reconstructed_image = Image.merge("RGB", (red_channel, green_channel, blue_channel))
    reconstructed_image.show()
    #path = "no_err.jpg"
    #path1 = "this-work-Lisa-0.01-del.jpg"
    path = "no_err.jpg"
    path1 = "1this-work-Lisa-0.01-sus-del.jpg"
    reconstructed_image.save(path1)
    # 测出所有图像的数据
    image_quantity_test(path, path1)
    #inhance_image(path1)
    #path2 =  "2.jpg"
    #image_quantity_test(path, path2)



# 错误加入
def error_add(data,error_rate):
    # 数据转为list数据类型
    dna_sequences = [list(sequence) for sequence in data]
    # 测量链长度*错误率
    #chosen_count = int(len(dna_sequences) * error_rate)
    chosen_count = int(len(dna_sequences)*len(dna_sequences[0])*error_rate)
    # dna_sequences = random.sample(copy.deepcopy(encoded_data["dna"]), chosen_count)
    total_indices = [sequence_index for sequence_index in range(len(dna_sequences))]
    for insertion_iteration in range(int(chosen_count*0.1)):
        chosen_index = random.choice(total_indices)

        dna_sequences[chosen_index].insert(random.randint(0, len(dna_sequences[chosen_index]) - 1),
                                           random.choice(['A', 'C', 'G', 'T']))

    # mutation errors
    for mutation_iteration in range(int(chosen_count*0.8)):
        chosen_index = random.choice(total_indices)
        chosen_index_in_sequence = random.randint(0, len(dna_sequences[chosen_index]) - 1)
        chosen_nucleotide = dna_sequences[chosen_index][chosen_index_in_sequence]
        dna_sequences[chosen_index][chosen_index_in_sequence] = \
            random.choice(list(filter(lambda nucleotide: nucleotide != chosen_nucleotide,
                                      ['A', 'C', 'G', 'T'])))

    # deletion errors
    for deletion_iteration in range(int(chosen_count*0.1)):
        chosen_index = random.choice(total_indices)
        del dna_sequences[chosen_index][random.randint(0, len(dna_sequences[chosen_index]) - 1)]

    # 将数据转为【'ATATATT'，'ATATATT'】
    last_dna_sequences = [''.join(sequence) for sequence in dna_sequences]
    return last_dna_sequences

def is_string(data):
    return isinstance(data, str)

# path---图像的文件名
def pipeline(path, error_rate):
    # 获得三通道二进制数据
    iamge_channel_bits = []
    compressed_data, compressed_data1, compressed_data2 = image_to_bit(path)
    iamge_channel_bits.append(compressed_data)
    iamge_channel_bits.append(compressed_data1)
    iamge_channel_bits.append(compressed_data2)

    dna_segments_all = []
    bits_segment_all = []
    miss_index_all = []
    for i in range(3):
        # 将数据进行分割，并转为list数据类型
        bit_segments = bits_to_list(iamge_channel_bits[i], 256)
        a = bit_segments[-1][:256]
        # 添加索引
        bit_segments = connect_all(bit_segments, 16)   #16 索引的长度
        a1 = bit_segments[-1][:16]
        # 编码
        dna_segments = encode(bit_segments)
        error_prone(dna_segments)
        # 错误加入
        dna_segments = error_add(dna_segments, error_rate)
        dna_segments_all.append(dna_segments)
    for i in range(3):
        bits_segments, miss_index = decode(dna_segments_all[i])
        bits_segments1 = ""
        for j in range(len(bits_segments)):
            if is_string(bits_segments[j]) and len(bits_segments[j]) == 256:
                temp = bits_segments[j][:256]
                # 判断
                if len(bits_segments[j]) < 256:
                    lack_length = 256 - len(bits_segments[j])
                    for s in range(lack_length):
                        temp += "0"
                bits_segments1 += temp
            else:
                string_list = [str(i) for i in bits_segments[j]]
                # 将字符串列表连接成一个字符串
                temp = ''.join(string_list)
                temp = temp[:256]

                if len(bits_segments[j]) < 256:
                    lack_length = 256 - len(bits_segments[j])
                    for s in range(lack_length):
                        temp += "0"
                bits_segments1 += temp
        bits_segment_all.append(bits_segments1)
        miss_index_all.append(miss_index)
    # 将二进制数据转为  图像部分

    bit_to_image("Lisa.jpg", bits_segment_all[0], bits_segment_all[1], bits_segment_all[2], miss_index_all)
pipeline("Lisa.jpg", 0.02)










