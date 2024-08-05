from PIL import Image
import numpy as np
from skimage.metrics import structural_similarity as ssim
import cv2
from skimage.metrics import structural_similarity as ssim
from skimage.transform import resize
def psnr(img1, img2):
    mse = np.mean((img1-img2)**2)
    if mse == 0:
        return float('inf')
    else:
        return 20*np.log10(255/np.sqrt(mse))

def MSE(img1,img2):
    mse = np.mean( (img1 - img2) ** 2)
    return mse

def ssim1 (img1 ,img2):
    red_channel, green_channel, blue_channel = img1.split()
    img_np = np.array(red_channel)
    img_np1 = np.array(green_channel)
    img_np2 = np.array(blue_channel)
    red_channel2, green_channel2, blue_channel2 = img2.split()
    img_np_2 = np.array(red_channel2)
    img_np1_2 = np.array(green_channel2)
    img_np2_2 = np.array(blue_channel2)
    ssim_score1 = ssim(img_np, img_np_2, multichannel=True)
    ssim_score2 = ssim(img_np1, img_np1_2, multichannel=True)
    ssim_score3 = ssim(img_np2, img_np2_2, multichannel=True)
    ssim_score = (ssim_score1 + ssim_score2 + ssim_score3) / 3
    return ssim_score




def compute_msssim(path, path1):
    """
    计算两个图像之间的MS-SSIM。
    :param image1: 第一个图像。
    :param image2: 第二个图像。
    :return: MS-SSIM 分数。
    """
    # 读取图像
    image1 = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
    image2 = cv2.imread(path1, cv2.IMREAD_GRAYSCALE)
    msssim_score = 0
    weights = [0.1, 0.3, 0.5, 0.7, 0.9]  # 这些权重可以根据需要调整

    for weight in weights:
        # 为每个尺度调整图像大小
        resized_image1 = resize(image1, (int(image1.shape[0] * weight), int(image1.shape[1] * weight)),
                                anti_aliasing=True)
        resized_image2 = resize(image2, (int(image2.shape[0] * weight), int(image2.shape[1] * weight)),
                                anti_aliasing=True)
        # 计算并累加当前尺度的SSIM
        score, _ = ssim(resized_image1, resized_image2, data_range=resized_image1.max() - resized_image1.min(),
                        full=True)
        msssim_score += score / len(weights)
    return msssim_score

def image_quantity_test(path, path1):

    img1 = np.array(Image.open(path))
    img2 = np.array(Image.open(path1))
    c1 = Image.open(path1)  # 在测量ssim时，只需要打开图像即可
    c2 = Image.open(path)

    # 计算两张彩色图像之间的 SSIM  MES  PSNR
    mse = MSE(img1, img2)
    print(mse, "mse")
    ssim_score = ssim1(c1, c2)
    print("SSIM Score:", ssim_score)
    psnr_score = psnr(img1, img2)
    print(psnr_score, "psnr")
    msssim_score = compute_msssim(path,path1)
    print(msssim_score, "msssim_score")


if __name__ == "__main__":
    path = "no_err.jpg"
    #path = "Lisa.jpg"
    path1 = "0.02.jpg"    #图像的名称   如"Lisa.jpg"
    image_quantity_test(path, path1)