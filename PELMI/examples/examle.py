import generalizedreedsolo
import numpy as np
import random
def main():
    show__example()


def cut(obj, sec):  # 按一定长度划分
    return [obj[i:i + sec] for i in range(0, len(obj), sec)]


def two_baseconversion(bitlist):  # 2进制转10进制
    tenbaseconversion = []
    for str_1 in bitlist:
        tenbaseconversion.append(int(str_1, 2))
    return tenbaseconversion


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
                                                        payload_length=32,  # 原始消息长度
                                                        p_factor=1)  # 加速因子


def disruption_model_1(length, error_rate, listvalue):  # 加入干扰
    disruptionnumber = length * error_rate  # 错误数量
    # print(listvalue)    非常棒的一点是   0。95代表替换错误的
    for i in range(int(disruptionnumber * 0.95)):
        strand = random.randint(0, len(listvalue) - 1)  # 序列随机
        location = random.randint(0, len(listvalue[strand]) - 1)  # 错误位置随机
        number = random.randint(0, 3)  # 碱基随机
        if number == 0:
            listvalue[strand][location] = 'A'  # 替换错误
        elif number == 1:
            listvalue[strand][location] = 'T'
        elif number == 2:
            listvalue[strand][location] = 'C'
        elif number == 3:
            listvalue[strand][location] = 'G'
    for j in range(int(disruptionnumber * 0.05)):
        strand = random.randint(0, len(listvalue) - 1)  # 序列随机
        location = random.randint(0, len(listvalue[strand]) - 1)  # 错误位置随机
        number_1 = random.randint(0, 3)  # 数字随机
        errortype = random.randint(1, 2)  # 错误类型随机
        if errortype == 1:  # 插入错误
            if number_1 == 0:
                listvalue[strand].insert(location, 'A')
            elif number_1 == 1:
                listvalue[strand].insert(location, 'T')
            elif number_1 == 2:
                listvalue[strand].insert(location, 'C')
            elif number_1 == 3:
                listvalue[strand].insert(location, 'G')
        elif errortype == 2:  # 删除错误
            del listvalue[strand][location]
    return listvalue


def model_1(length, msg, binarylen, ori_base_msg):
    for error_rate in range(1, 7, 1):
        dis_msg = disruption_model_1(length, error_rate * 0.01, msg)         # 错误的随机加入
        # print(dis_msg)
        # *****************纠错****************
        finalfile = []
        for dis_i in range(len(dis_msg)):
            sup_msg = ''
            for dis_j in dis_msg[dis_i]:
                if dis_j == 'A':
                    sup_msg = sup_msg + '00'
                elif dis_j == 'T':
                    sup_msg = sup_msg + '01'
                elif dis_j == 'C':
                    sup_msg = sup_msg + '10'
                elif dis_j == 'G':
                    sup_msg = sup_msg + '11'
            finalfile.append(sup_msg)
        # print(finalfile, "测试", len(finalfile))
        finalfile_ten = []
        for dis_l in finalfile:
            tempbitmsg = cut(dis_l, binarylen)
            finalfile_ten.append(two_baseconversion(tempbitmsg))
        # print(finalfile_ten)    将二进制数据转为十进制数据
        correct = []
        for dis_m in finalfile_ten:
            model = reed_solomon_model(binarylen, len(dis_m))
            # noinspection PyBroadException
            try:
                correct.append(model.decode(dis_m))    # 去掉纠错吗的  像素值
            except Exception:
                correct.append(([]))
                # print("纠错失败的序列号:", finalfile_ten.index(dis_m))

        # print(correct)
        error_correction_rate(error_rate, binarylen, correct, ori_base_msg)


def mapp_base(data):
    data_base = []
    for split_num in range(len(data)):
        data_base.append([])
        for i in range(0, len(data[split_num]), 2):
            if data[split_num][i:i + 2] == '00':
                data_base[split_num].append('A')
            elif data[split_num][i:i + 2] == '01':
                data_base[split_num].append('T')
            elif data[split_num][i:i + 2] == '10':
                data_base[split_num].append('C')
            elif data[split_num][i:i + 2] == '11':
                data_base[split_num].append('G')
    # print(len(data_base))
    return data_base


def error_correction_rate(error_rate, binarylen, listvalue, msg):
    base_orimsg = msg  # 原始碱基序列
    length = 0
    for base_split in base_orimsg:
        length = length + len(base_split)

    file_3 = open('temp.txt', 'a')
    file_3.seek(0)
    file_3.truncate()
    for temp_split in listvalue:
        if len(temp_split) != 0:
            for temp_value in temp_split:
                file_3.write(to_bin(temp_value, binarylen))
            file_3.write('\n')
        else:
            file_3.write('' + '\n')
    file_3.close()

    file_3 = open('temp.txt', 'r')
    data_decode = [line.strip('\n') for line in file_3.readlines()]  # 按行读取
    file_3.close()
    # print(data_decode)

    base_decode = []
    for temp_split_1 in data_decode:
        if len(temp_split_1) != 0:
            temp_list = []
            for i in range(0, len(temp_split_1), 2):
                if temp_split_1[i:i + 2] == '00':  # 映射成碱基
                    temp_list.append('A')
                elif temp_split_1[i:i + 2] == '01':
                    temp_list.append('T')
                elif temp_split_1[i:i + 2] == '10':
                    temp_list.append('C')
                elif temp_split_1[i:i + 2] == '11':
                    temp_list.append('G')
            base_decode.append(temp_list)
            # print(mapp_base(temp_split_1))
        else:
            base_decode.append([])
    # print(base_decode)

    right_num = 0
    for i in range(len(base_orimsg)):
        if len(base_orimsg[i]) >= len(base_decode[i]) and len(base_decode[i]) != 0:
            for j in range(len(base_decode[i])):
                if base_decode[i][j] == base_orimsg[i][j]:
                    right_num = right_num + 1
        elif len(base_orimsg[i]) < len(base_decode[i]):
            for k in range(len(base_orimsg[i])):
                if base_decode[i][k] == base_orimsg[i][k]:
                    right_num = right_num + 1
    c = np.abs(right_num / length)
    print('错误率:{:.1%}'.format(error_rate * 0.01), '纠错率:{:.4%}'.format(c))


def show__example():
    f = open('handle.txt', 'r')
    data = [line.strip('\n') for line in f.readlines()]  # 按行读取
    f.close()
    file_1 = open('binfile.txt', 'a')
    file_1.seek(0)
    file_1.truncate()
    file_2 = open('10file.txt', 'a')
    file_2.seek(0)
    file_2.truncate()
    tempbinarymsg = ''
    for binstr in data:
        tempbinarymsg = tempbinarymsg + binstr
    print("2进制原始消息长度:", len(tempbinarymsg))

    # 上面的处理方式    就是将文件转为字符的二进制数据    tempbinarymsg

    #
    binarylen = 0
    for i in range(2, 11):
        if len(tempbinarymsg) % i == 0 and i > binarylen:
            binarylen = i
    print("按", binarylen, "位读取,转换成10进制")
    temp8bitmsg = cut(tempbinarymsg, binarylen)
    # print(temp8bitmsg)#此处的tempbinarymsg刚好被8整除
    orimsg = two_baseconversion(temp8bitmsg)
    # print(orimsg)
    print("10进制原始消息长度:", len(orimsg))
    # 前面的逻辑是   将二进制数据转为像素值形式

    # 下面的逻辑是 将数据分为36个一组，形成二维数组
    temporimsg = cut(orimsg, 32)
    # temporimsg = cut(orimsg, 20)
    # print(temporimsg)
    # R S 进行编码的部分     binarylen指的是八个一组
    test_normal_rs = reed_solomon_model(binarylen, 34)  # 生成伽罗华域模型
    # test_normal_rs = reed_solomon_model(binarylen, 24)
    normal_msg = []
    for templist in temporimsg:
        file_2.write(str(temporimsg.index(templist)) + ' ' + str(templist) + '\n')
        if len(templist) == 32:
        # if len(templist) == 20:
            normal_msg.append(test_normal_rs.encode(templist))
        else:
            supplementarybit = '0' * (32 - len(templist))
            # supplementarybit = '0' * (20 - len(templist))
            templist.extend(map(int, list(supplementarybit)))
            normal_msg.append(test_normal_rs.encode(templist))
    # print(normal_msg)     上面完成编码部分
    file_2.close()
    length = 0
    for msg in normal_msg:
        length = length + len(msg)
        for value_1 in msg:
            file_1.write(to_bin(value_1, binarylen))
        file_1.write('\n')
    file_1.close()
    print("加密后二进制文件已在“binfile.txt”中生成")
    print("加密后10进制消息长度:", length)
    c = np.abs((length - len(orimsg)) / length)
    print('校验位占比:{:.4%}'.format(c))

    ori_base_msg = mapp_base(cut(tempbinarymsg, 256))# 原始碱基序列
    # ori_base_msg = mapp_base(cut(tempbinarymsg, 160))
    # print(ori_base_msg)

    file_5 = open('binfile.txt', 'r')
    data = [line.strip('\n') for line in file_5.readlines()]  # 按行读取
    file_5.close()
    base = mapp_base(data)  # 映射成碱基
    # print(base)
    num = 0
    for i in base:
        num = num + len(i)  # 碱基总数
    # 加入干扰 计算纠错率
    model_1(num, base, binarylen, ori_base_msg)


if __name__ == "__main__":
    main()

# 今天的任务
# 加入RS的逻辑    大致捋顺了，进行基本的修改  就可以啦
# 随机加入错误的逻辑
