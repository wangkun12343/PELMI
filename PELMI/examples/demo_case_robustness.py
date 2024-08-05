import os
import Chamaeleo
from Chamaeleo.methods.default import BaseCodingAlgorithm
from Chamaeleo.methods.ecc import Hamming, ReedSolomon
from Chamaeleo.methods.fixed import Church
from Chamaeleo.utils.pipelines import RobustnessPipeline
import os
import Chamaeleo
from Chamaeleo.methods.default import BaseCodingAlgorithm
from Chamaeleo.methods.fixed import Church, Goldman, Grass, Blawat
from Chamaeleo.methods.flowed import DNAFountain, YinYangCode
from Chamaeleo.utils.pipelines import BasicFeaturePipeline

if __name__ == "__main__":
    root_path = os.path.dirname(Chamaeleo.__file__)
    #error_corrections = ReedSolomon()
    file_paths =  {
        "image1.jpg": os.path.join(root_path, "data", "pictures", "Lisa.jpg")
    }
    coding_schemes = {
        #"DNAFountain et al.": YinYangCode()
        "Grass et al.": Grass()
    }
    error_corrections = {
        "ReedSolomon": ReedSolomon()
    }

    needed_indices = [
        True

    ]
    out_path = "Lisa7.jpg"

    pipeline = RobustnessPipeline(
        coding_schemes=coding_schemes,
        error_corrections=error_corrections,
        needed_indices=needed_indices,
        file_paths=file_paths,
        nucleotide_insertion=0.000,
        nucleotide_mutation=0.004,
        nucleotide_deletion=0.000,
        output_path=out_path,
        sequence_loss=0.00,
        iterations=3,
        segment_length= 264,
        index_length=16,
        need_logs=True
    )

    pipeline.evaluate()
    pipeline.output_records(type="string")
