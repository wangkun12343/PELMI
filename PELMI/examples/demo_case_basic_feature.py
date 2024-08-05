import os
import Chamaeleo
from Chamaeleo.methods.default import BaseCodingAlgorithm
from Chamaeleo.methods.fixed import Church, Goldman, Grass, Blawat
from Chamaeleo.methods.flowed import DNAFountain, YinYangCode
from Chamaeleo.utils.pipelines import BasicFeaturePipeline

"""
下午的任务：运行通这个程序，理解其底层逻辑
step1 ：整理清楚代码的功能部分
（1）学会计算不同编码方案编码的基本特征
（2）choice 选择在不同编码方案中找到最佳结果
（3）conmbined 学会比较不同的编码方案
（4）individual 了解如何完成编码和解码的过程

step2 ：运行整体程序，加入错误，获得结果
问题一：如何加入错误   解决
问题二：如何输出 错误的图像

"""


if __name__ == "__main__":
    root_path = os.path.dirname(Chamaeleo.__file__)

    file_paths = {
        "Mona Lisa.jpg": os.path.join(root_path, "data", "pictures", "Mona Lisa.jpg")
    }

    coding_schemes = {
        "Base": BaseCodingAlgorithm(),
        "Church et al.": Church(), "Goldman et al.": Goldman(), "Grass et al.": Grass(), "Blawat et al.": Blawat(),
        "DNA Fountain": DNAFountain(redundancy=0.5), "Yin-Yang Code": YinYangCode()
    }
    needed_indices = [
        True,
        True, True, True, True,
        False, True
    ]

    pipeline = BasicFeaturePipeline(
        coding_schemes=coding_schemes,
        needed_indices=needed_indices,
        file_paths=file_paths,
        segment_length=128,
        index_length=16,
        need_logs=True)

    pipeline.calculate()
    pipeline.output_records(type="string")
