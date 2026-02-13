from skimage import io as skio
import numpy as np
from cellpose import models, core, io as cpio
from pathlib import Path
import os
import gc


def _readtiffs(rawdir: Path) -> dict[str, np.ndarray]:
	rawimgs = {}

	for file in os.listdir(rawdir):
		try: 
			label = Path(file).stem.replace('.ome', '')
			img = skio.imread(rawdir / file)[1] # get second channel
			rawimgs[label] = img
		except:
			print(f"\tCould not read \"{file}\"")
			continue
	
	return rawimgs


# for segcellpose(): implement model.eval() settings into params?
def segcellpose(rawdir: Path, segdir: Path, modeldir: Path, flow: float = 0.7, cellprob: float = 0.2, norm: dict = {"tile_norm_blocksize": 0}) -> None:
	os.environ["CELLPOSE_LOCAL_MODELS_PATH"] = str(modeldir)
	os.makedirs(modeldir, exist_ok=True)
	os.makedirs(segdir, exist_ok=True)

	# cpio.logger_setup()

	if core.use_gpu() == False:
		raise ImportError("No GPU access, change your runtime")

	print("Reading TIFFs...")
	rawimgs = _readtiffs(rawdir)
	
	print("Loading model...")
	model = models.CellposeModel(gpu=True)

	print("Segmenting...")
	for label, rawimg in rawimgs.items():
		masks, _, _ = model.eval(rawimg, batch_size=32, 
									flow_threshold=flow, 
									cellprob_threshold=cellprob,
									normalize=norm)
		gc.collect()
		print(f"\tSuccessfully segmented {label}")
		skio.imsave(segdir / f"{label}_seg.tiff", masks, check_contrast=False)


if __name__ == "__main__":
	FLUORDIR = Path("/Users/kvying/Files/ucsb-local/rothman-sam/fluorescence-pipeline")
	DATADIR = Path("/Users/kvying/Files/ucsb-local/rothman-data/embryo")
	RAWDIR = DATADIR / "11-13-2025"
	OUTDIR = RAWDIR / "out"
	SEGDIR = OUTDIR / "seg"
	MODELDIR = FLUORDIR / "model"

	flow = 0.7
	cellprob = 0.2

	segcellpose(RAWDIR, SEGDIR, MODELDIR, flow=flow, cellprob=cellprob)