from PIL import Image
import numpy as np
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

def index_map_dna(binary_data):
    mapping = {
        '00': 'A',
        '01': 'T',
        '10': 'C',
        '11': 'G'
    }
    dna_data = ""
    for i in range(0,len(binary_data), 2):
        dna_data += mapping[binary_data[i: i +2]]
    return dna_data

def dna_map_index(dna_data):
    mapping = {
        'A': '00',
        'T': '01',
        'C': '10',
        'G': '11'
    }
    bit_data = ""
    for i in range(0,len(dna_data)):
        bit_data += mapping[dna_data[i: i +1]]
    return bit_data

def encode_map_rule(split_sub , pre_single_dna , rule):
    if rule == 0 or rule == 1:
        if pre_single_dna == "A":
            if split_sub == "00":
                pre_single_dna1 = "C"
            elif split_sub == "01":
                pre_single_dna1 = "G"
            else:
                pre_single_dna1 = "T"
            return pre_single_dna1

        if pre_single_dna == "T":
            if split_sub == "00":
                pre_single_dna1 = "A"
            elif split_sub == "01":
                pre_single_dna1 = "C"
            else:
                pre_single_dna1 = "G"
            return pre_single_dna1

        if pre_single_dna == "G":
            if split_sub == "00":
                pre_single_dna1 = "T"
            elif split_sub == "01":
                pre_single_dna1 = "A"
            else:
                pre_single_dna1 = "C"
            return pre_single_dna1

        if pre_single_dna == "C":
            if split_sub == "00":
                pre_single_dna1 = "G"
            elif split_sub == "01":
                pre_single_dna1 = "T"
            else:
                pre_single_dna1 = "A"
            return pre_single_dna1

    if rule == 1:
        if pre_single_dna == "A":
            if split_sub == "00":
                pre_single_dna1 = "C"
            elif split_sub == "01":
                pre_single_dna1 = "G"
            else:
                pre_single_dna1 = "T"
            return pre_single_dna1

        if pre_single_dna == "T":
            if split_sub == "00":
                pre_single_dna1 = "A"
            elif split_sub == "01":
                pre_single_dna1 = "C"
            else:
                pre_single_dna1 = "G"
            return pre_single_dna1

        if pre_single_dna == "G":
            if split_sub == "00":
                pre_single_dna1 = "T"
            elif split_sub == "01":
                pre_single_dna1 = "A"
            else:
                pre_single_dna1 = "C"
            return pre_single_dna1

        if pre_single_dna == "C":
            if split_sub == "00":
                pre_single_dna1 = "G"
            elif split_sub == "01":
                pre_single_dna1 = "T"
            else:
                pre_single_dna1 = "A"
            return pre_single_dna1

    if rule == 2 or rule == 3:
        if pre_single_dna == "A":
            if split_sub == "10":
                pre_single_dna1 = "C"
            elif split_sub == "11":
                pre_single_dna1 = "G"
            else:
                pre_single_dna1 = "T"
            return pre_single_dna1

        if pre_single_dna == "T":
            if split_sub == "10":
                pre_single_dna1 = "A"
            elif split_sub == "11":
                pre_single_dna1 = "C"
            else:
                pre_single_dna1 = "G"
            return pre_single_dna1

        if pre_single_dna == "G":
            if split_sub == "10":
                pre_single_dna1 = "T"
            elif split_sub == "11":
                pre_single_dna1 = "A"
            else:
                pre_single_dna1 = "C"
            return pre_single_dna1

        if pre_single_dna == "C":
            if split_sub == "10":
                pre_single_dna1 = "G"
            elif split_sub == "11":
                pre_single_dna1 = "T"
            else:
                pre_single_dna1 = "A"
            return pre_single_dna1

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


def bit_to_image(img, compressed_data , compressed_data1 ,compressed_data2):
    img = Image.open(img)
    red_channel, green_channel, blue_channel= img.split()
    img_np = np.array(red_channel)
    quantized_coef_1d = np.array([int(compressed_data[i:i + 8], 2) for i in range(0, len(compressed_data), 8)])
    quantized_coef_1d1 = np.array([int(compressed_data1[i:i + 8], 2) for i in range(0, len(compressed_data1), 8)])
    quantized_coef_1d2 = np.array([int(compressed_data2[i:i + 8], 2) for i in range(0, len(compressed_data2), 8)])
    # 将量化系数重构为二维数组
    quantized_coef = np.reshape(quantized_coef_1d, img_np.shape)
    quantized_coef = quantized_coef.astype(np.uint8)
    quantized_coef1 = np.reshape(quantized_coef_1d1, img_np.shape)
    quantized_coef1 = quantized_coef1.astype(np.uint8)
    quantized_coef2 = np.reshape(quantized_coef_1d2, img_np.shape)
    quantized_coef2 = quantized_coef2.astype(np.uint8)
    # 将numpy 转为通道数据
    reconstructed_red_channel = Image.fromarray(quantized_coef)
    reconstructed_red_channel1 = Image.fromarray(quantized_coef1)
    reconstructed_red_channel2 = Image.fromarray(quantized_coef2)
    # 数据合并
    reconstructed_image = Image.merge("RGB", (reconstructed_red_channel, reconstructed_red_channel1, reconstructed_red_channel2))
    reconstructed_image.show()
    path1 = "1last-hl-Lisa-sus-ind-0.01.jpg"
    path = "Lisa.jpg"

    reconstructed_image.save(path1)
    # 将均值迭代的图像与原始图像进行测量
    image_quantity_test(path1, path)
    inhance_image(path1)
    path3 = "2.jpg"
    image_quantity_test(path3, path)


def encode_core(binary , rule):
    # 判断是否为0
    dna_data = []
    for i in range(len(binary)):
        if len(binary[i]) % 2 != 0:
            binary[i] = binary[i] + "0"
        binary_index = binary[i][:16]
        binary_data = binary[i][16:]
        # 将索引部分直接编码
        index_dna_temp = index_map_dna(binary_index)
        dna_data_block = ""
        dna_data_block_all = ""
        pre_single_dna = "A"
        for j in range(0, len(binary_data), 2):
            # 判断索引是否是否为0
            if j == 0:
            # 最开始默认值是A
                split_sub = binary_data[j: j + 2]
                single_dna = encode_map_rule(split_sub, pre_single_dna, rule)     # 对应规则
                dna_data_block = dna_data_block + single_dna
            else:
                pre_single_dna = single_dna
                split_sub = binary_data[j: j + 2]
                single_dna = encode_map_rule(split_sub, pre_single_dna, rule)
                dna_data_block = dna_data_block + single_dna
        dna_data_block_all = index_dna_temp + dna_data_block
        dna_data.append(dna_data_block_all)
    return dna_data


# 将十进制数据转为16位的二进制数据
def decimal_to_16_bit_binary_string(decimal_num):
    binary_str = bin(decimal_num)[2:].zfill(
        16)  # bin()函数将十进制数转换为二进制字符串，并添加了'0b'前缀。我们通过[2:]移除这个前缀，然后使用zfill()函数将字符串填充为16位。
    return binary_str


def check_rule(content):
    odd_data = ""
    mate_data = ""
    rule = 0
    """
    for i in range(0, len(content), 2):
        odd_data = odd_data + content[i]
    for i in range(0, len(content), 2):
        mate_data = mate_data + content[i]
    if odd_data.count("0") > len(odd_data) / 2:
        if odd_data.count("0") > len(odd_data) / 2:
            rule = 0
        else:
            rule = 1
    else:
        if odd_data.count("0") > len(odd_data) / 2:
            rule = 2
        else:
            rule = 3
    """
    for i in range(0, len(content), 2):
        odd_data = odd_data + content[i]
        mate_data = mate_data + content[i + 1]
    if odd_data.count("0") > len(odd_data) / 2:
        if odd_data.count("0") > len(odd_data) / 2:
            rule = 0
        else:
            rule = 1
    else:
        if mate_data.count("0") > len(odd_data) / 2:
            rule = 2
        else:
            rule = 3

    return rule

def add_index(content, segment_length):
    if len(content) % segment_length != 0:
        temp = len(content) % segment_length
        content = pad_string_right(content, temp)
    data = []
    t = 0
    for i in range(0, len(content), segment_length):
        index_emp = decimal_to_16_bit_binary_string(t)
        temp_data = index_emp + content[i: i + segment_length]
        data.append(temp_data)
        t += 1
    return data

# 添加标识位
def add_flag(dna_data):
    divide_data = []
    for i in range(len(dna_data)):
        temp = dna_data[i]
        temp_block = ""
        for j in range(3):
            temp_block += temp[j*34: j*34 + 34] + "AA"
        temp_block += temp[-34:]
        divide_data.append(temp_block)
    return divide_data


# 解码规则
def decode_map_rule(split_sub , pre_single_dna , rule):
    if rule == 0 or rule == 1:
        if pre_single_dna == "A":
            if split_sub == "C":
                pre_single_dna1 = "00"
            elif split_sub == "G":
                pre_single_dna1 = "01"
            else:
                if rule == 0:
                    pre_single_dna1 = "10"
                if rule == 1:
                    pre_single_dna1 = "11"
            return pre_single_dna1

        if pre_single_dna == "T":
            if split_sub == "A":
                pre_single_dna1 = "00"
            elif split_sub == "C":
                pre_single_dna1 = "01"
            else:
                if rule == 0:
                    pre_single_dna1 = "10"
                if rule == 1:
                    pre_single_dna1 = "11"
            return pre_single_dna1

        if pre_single_dna == "G":
            if split_sub == "T":
                pre_single_dna1 = "00"
            elif split_sub == "A":
                pre_single_dna1 = "01"
            else:
                if rule == 0:
                    pre_single_dna1 = "10"
                if rule == 1:
                    pre_single_dna1 = "11"
            return pre_single_dna1

        if pre_single_dna == "C":
            if split_sub == "G":
                pre_single_dna1 = "00"
            elif split_sub == "T":
                pre_single_dna1 = "01"
            else:
                if rule == 0:
                    pre_single_dna1 = "10"
                if rule == 1:
                    pre_single_dna1 = "11"
            return pre_single_dna1

    if rule == 2 or rule == 3:
        if pre_single_dna == "A":
            if split_sub == "C":
                pre_single_dna1 = "10"
            elif split_sub == "G":
                pre_single_dna1 = "11"
            else:
                if rule == 2:
                    pre_single_dna1 = "00"
                if rule == 3:
                    pre_single_dna1 = "01"
            return pre_single_dna1

        if pre_single_dna == "T":
            if split_sub == "A":
                pre_single_dna1 = "10"
            elif split_sub == "C":
                pre_single_dna1 = "11"
            else:
                if rule == 2:
                    pre_single_dna1 = "00"
                if rule == 3:
                    pre_single_dna1 = "01"
            return pre_single_dna1

        if pre_single_dna == "G":
            if split_sub == "T":
                pre_single_dna1 = "10"
            elif split_sub == "A":
                pre_single_dna1 = "11"
            else:
                if rule == 2:
                    pre_single_dna1 = "00"
                if rule == 3:
                    pre_single_dna1 = "01"
            return pre_single_dna1

        if pre_single_dna == "C":
            if split_sub == "G":
                pre_single_dna1 = "10"
            elif split_sub == "T":
                pre_single_dna1 = "11"
            else:
                if rule == 2:
                    pre_single_dna1 = "00"
                if rule == 3:
                    pre_single_dna1 = "01"
            return pre_single_dna1
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
    for i in range(len(dna_segments)):
         dna_segment = ""
         t = 34
         for j in range(3):
            if dna_segments[i][t:t + 2] == "AA":
                dna_segment += dna_segments[i][t - 34:t]
            elif dna_segments[i][t - 1:t + 1] == "AC":
                dna_segment += dna_segments[i][t - 33:t]
                dna_segment += "A"
                test1_cor += 1
                test1_cor_list.append(i)
                print(dna_segments[i], i , len(dna_segments[i]), "1")
            elif dna_segments[i][t + 1:t + 3] == "AC":
                dna_segment += dna_segments[i][t - 34:t]
                test1_cor += 1
                test1_cor_list.append(i)
                print(dna_segments[i], i, len(dna_segments[i]), "2")
            else:
                dna_segment += dna_segments[i][t - 34:t]
                test1_block.append(dna_segments[i][t - 34:t])
                test1_data.append(dna_segments[i])
                test1_no += 1
                test1_no_list.append(i)
            t += 36
         dna_segment += dna_segments[i][-34 :]
         length = len(dna_segment)
         if length != 136:
             print(dna_segments[i], "dna_segments")
             print(dna_segment, "dna_segment")
         remove_dna_segment.append(dna_segment)
    return  remove_dna_segment

# 解码
def decode_core(data, rule):
    """
    解码逻辑需要修改：针对DNA链进行处理
    编码方面同样需要修改
    数据全部按照：['ATATCG','ATATCG'],这种数据类型
    """
    last_data = []
    for i in range(len(data)):
        binary_bit_data = ""
        # 将数据部分分割成索引部分和index部分
        index_temp = data[i][:8]
        data[i] = data[i][8:]
        # 对索引部分进行解码
        index_bit_temp = dna_map_index(index_temp)

        # 对数据部分解码
        for j in range(0, len(data[i])):
            if j == len(data[i]) - 1:
                pre_single_dna = "A"
                split_two = decode_map_rule(data[i][len(data[i]) - 1 - j], pre_single_dna, rule)
            else:
                pre_single_dna = data[i][len(data[i]) - 2 - j]
                split_two = decode_map_rule(data[i][len(data[i]) - 1 - j], pre_single_dna, rule)
            binary_bit_data = split_two + binary_bit_data
        binary_bit_data = index_bit_temp + binary_bit_data
        last_data.append(binary_bit_data)
    # 去除索引
    #binary_last_data = ""
    #for i in range(0, len(binary_bit_data), 272):
    #    temp = binary_bit_data[i: i + 272]
    #    binary_last_data += temp[16:]
    return last_data


def pad_string_right(string, length):
    # 如果字符串长度已经大于等于指定的长度，则不需要补充0
    if len(string) >= length:
        return string
    else:
        # 计算需要补充的0的数量
        num_zeros = length - len(string)
        # 在字符串右侧补充0
        return string + '0' * num_zeros


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
    length = len(data)
    filtered_list1 = []
    for i in range(len(filtered_list)):
        index_demical = int(''.join(map(str, filtered_list[i][:16])), 2)
        if index_demical <= length:
            filtered_list1.append(filtered_list[i])
    return filtered_list1
# 对解码失败的链进行补充
def data_subment(nested_list):
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
                #data_sub.append(list)
        else:
            list = []
            length1 = len(sorted_list[0])
            for i in range(length1):
                temp = random.choice([0, 1])
                list.append(temp)
            data_sub.append(list)
    print(z)
    return data_sub, miss_index

# 错误加入
def error_add(data,error_rate):
    # 数据转为list数据类型
    dna_sequences = [list(sequence) for sequence in data]
    # 测量链长度*错误率
    #chosen_count = int(len(dna_sequences) * error_rate)
    chosen_count = int(len(dna_sequences)*len(dna_sequences[0])*error_rate)
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


# path---图像的文件名
def pipeline(path, error_rate):
    #compressed_data, compressed_data1, compressed_data2, compressed_data3 = image_to_bit(path)
    compressed_data, compressed_data1, compressed_data2= image_to_bit(path)
    t0 = len(compressed_data)
    # step 2 :在图像中随机加入错误
    colord_bit = []
    channel_dna =""
    colord_bit.append(compressed_data)
    colord_bit.append(compressed_data1)
    colord_bit.append(compressed_data2)
    #colord_bit.append(compressed_data3)
    colord_channel = []
    for three_passage in range(len(colord_bit)):
        length = len(colord_bit[three_passage])
        #  判断使用那种规则
        rule = check_rule(colord_bit[three_passage])
        content = add_index(colord_bit[three_passage], 256) # 添加索引
        # step2 编码规则
        dna_data_last = encode_core(content, rule)
        # step3 添加标识位----返回的是【'TAAT'，'TAAT'】
        dna_data = add_flag(dna_data_last)
        # step4:加入错误
        dna_data = error_add(dna_data , error_rate)

        error_prone(dna_data)
        # 解码
        dna_data = remove_flag(dna_data)
        ceshi1 = []
        test1 = 0
        for i in range(len(dna_data)):
            if len(dna_data[i]) != 136:
                test1 += 1
                temp = len(dna_data[i])
                ceshi1.append(temp)
        bit_data = decode_core(dna_data, rule)
        # 对索引进行排序
        #  进行修改，解码失败的链进行补充, 并获取丢失链的索引
        bits_segments_all, miss_index = data_subment(bit_data)
        # 去除索引


        bit_data_last = ""
        test1 = 0
        for i in range(len(bits_segments_all)):
            if len(bits_segments_all[i]) == 256:
                test1 += 1
            bits_segments_all[i] = ''.join(map(str, bits_segments_all[i]))
            bits_segments_all[i] = bits_segments_all[i][-256:]
            bit_data_last = bit_data_last + bits_segments_all[i]
        # 将list数据转为str数据
        bit_data_last = bit_data_last[:length]
        colord_channel.append(bit_data_last)

    bit_to_image(path, colord_channel[0], colord_channel[1], colord_channel[2])

pipeline("Lisa.jpg", 0.0000)






