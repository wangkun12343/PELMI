from PIL import Image
import numpy as np


# 图像转为二进制数据
def image_to_numpy(img1):
    img = Image.open(img1)
    red_channel, green_channel, blue_channel = img.split()
    # img.show()
    # 转换为Numpy数组
    img_np = np.array(red_channel)
    img_np1 = np.array(green_channel)
    img_np2 = np.array(blue_channel)
    return img_np, img_np1, img_np2


def image_to_bit(img_np, img_np1, img_np2):
    # 转换为一维数组
    quantized_coef_1d = img_np.flatten()
    quantized_coef_1d1 = img_np1.flatten()
    quantized_coef_1d2 = img_np2.flatten()
    # 一维数组转换为二进制数据
    compressed_data = ''.join([format(int(x), '08b') for x in quantized_coef_1d])
    compressed_data1 = ''.join([format(int(x), '08b') for x in quantized_coef_1d1])
    compressed_data2 = ''.join([format(int(x), '08b') for x in quantized_coef_1d2])
    return compressed_data, compressed_data1, compressed_data2


def bit_to_image(shape, compressed_data, compressed_data1, compressed_data2):
    # 对数据进行补充
    wigth, height = shape
    sum = wigth * height * 8
    if len(compressed_data) < sum:
        compressed_data = str(compressed_data).ljust(sum, '0')
    if len(compressed_data1) < sum:
        compressed_data1 = str(compressed_data1).ljust(sum, '0')
    if len(compressed_data2) < sum:
        compressed_data2 = str(compressed_data2).ljust(sum, '0')

    quantized_coef_1d = np.array([int(compressed_data[i:i + 8], 2) for i in range(0, len(compressed_data), 8)])
    quantized_coef_1d1 = np.array([int(compressed_data1[i:i + 8], 2) for i in range(0, len(compressed_data1), 8)])
    quantized_coef_1d2 = np.array([int(compressed_data2[i:i + 8], 2) for i in range(0, len(compressed_data2), 8)])
    print("quantized_coef_1d", len(quantized_coef_1d))
    print("quantized_coef_1d1", len(quantized_coef_1d1))
    print("quantized_coef_1d2", len(quantized_coef_1d2))

    # 将量化系数重构为二维数组
    quantized_coef = np.reshape(quantized_coef_1d, shape)
    quantized_coef = quantized_coef.astype(np.uint8)
    print(quantized_coef, "quantized_coef")

    quantized_coef1 = np.reshape(quantized_coef_1d1, shape)
    quantized_coef1 = quantized_coef1.astype(np.uint8)

    quantized_coef2 = np.reshape(quantized_coef_1d2, shape)
    quantized_coef2 = quantized_coef2.astype(np.uint8)

    # 将numpy 转为通道数据
    reconstructed_red_channel = Image.fromarray(quantized_coef)
    reconstructed_red_channel1 = Image.fromarray(quantized_coef1)
    reconstructed_red_channel2 = Image.fromarray(quantized_coef2)

    # 数据合并
    reconstructed_image = Image.merge("RGB", (
    reconstructed_red_channel, reconstructed_red_channel1, reconstructed_red_channel2))
    # reconstructed_image = Image.merge("RGB", (red_channel, green_channel, blue_channel))
    reconstructed_image.show()
    reconstructed_image.save("church_colord_0.01.jpg")


def Enhance_image(data, window_row, window_column):
    data_min = np.min(data)
    data_max = np.max(data)
    data_med = np.median(data)
    data_sort = np.sort(data.reshape(-1))
    A2 = data_sort[1]
    A8 = data_sort[7]

    for i in range(window_row):
        for j in range(window_column):
            if (A2 < data[i, j] and data[i, j] <= data_med) or (data_med <= data[i, j] and data[i, j] < A8):
                continue
            else:
                data[i, j] = data_med
    return data


if __name__ == '__main__':
    # read data from image
    data_r, data_g, data_b = image_to_numpy('image.jpg')

    """
    enhance the image
    """
    # set the window size
    window_row = 3
    window_column = 3

    # slide without stack
    for i in range(0, data_r.shape[0] - window_row + 1, window_row):
        for j in range(0, data_r.shape[1] - window_column + 1, window_column):
            Enhance_image(data_r[i:i + window_row, j:j + window_column],
                                                                          window_row, window_column)
            Enhance_image(data_g[i:i + window_row, j:j + window_column],
                                                                          window_row, window_column)
            Enhance_image(data_b[i:i + window_row, j:j + window_column],
                                                                          window_row, window_column)
            # data_r[i:i + window_row, j:j + window_column] = Enhance_image(data_r[i:i + window_row, j:j + window_column],
            #                                                               window_row, window_column)
            # data_g[i:i + window_row, j:j + window_column] = Enhance_image(data_g[i:i + window_row, j:j + window_column],
            #                                                               window_row, window_column)
            # data_b[i:i + window_row, j:j + window_column] = Enhance_image(data_b[i:i + window_row, j:j + window_column],
            #                                                               window_row, window_column)

    data_r_bit, data_g_bit, data_b_bit = image_to_bit(data_r, data_g, data_b)
    bit_to_image(data_r.shape, data_r_bit, data_g_bit, data_b_bit)

