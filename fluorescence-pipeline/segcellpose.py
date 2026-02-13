from skimage import io as skio
import numpy as np
from cellpose import models, core, io as cpio
from pathlib import Path
import os
import gc


def readtiffs(rawdir: Path, channel: int = 0) -> dict[str, np.ndarray]:
	imgs = {}

	for file in os.listdir(rawdir):
		try: 
			label = '.'.join(Path(file).stem.split('.')[:-1]) # ex. "file.seg.tiff" -> "file"
			img = skio.imread(rawdir / file)
			
			if len(img.shape) > 3:
				img = img[channel]
			
			imgs[label] = img
		except:
			print(f"\tCould not read \"{file}\"")
			continue
	
	return imgs


# for segcellpose(): implement model.eval() settings into params?
def segcellpose(rawdir: Path, segdir: Path, modeldir: Path, channel: int = 1, flow: float = 0.7, cellprob: float = 0.2, norm: dict = {"tile_norm_blocksize": 0}) -> None:
	os.environ["CELLPOSE_LOCAL_MODELS_PATH"] = str(modeldir)
	os.makedirs(modeldir, exist_ok=True)
	os.makedirs(segdir, exist_ok=True)

	# cpio.logger_setup()

	if core.use_gpu() == False:
		raise ImportError("No GPU access, change your runtime")

	print("Reading TIFFs...")
	rawimgs = readtiffs(rawdir, channel=channel)
	
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
		skio.imsave(segdir / f"{label}.seg.tiff", masks, check_contrast=False)
	
	return


if __name__ == "__main__":
	FLUORDIR = Path("/Users/kvying/Files/ucsb-local/rothman-sam/fluorescence-pipeline")
	DATADIR = Path("/Users/kvying/Files/ucsb-local/rothman-data/embryo")
	RAWDIR = DATADIR / "11-13-2025"
	OUTDIR = RAWDIR / "out"
	SEGDIR = OUTDIR / "seg"
	MODELDIR = FLUORDIR / "model"

	flow = 0.7
	cellprob = 0.2

	segcellpose(RAWDIR, SEGDIR, MODELDIR, channel=1, flow=flow, cellprob=cellprob)