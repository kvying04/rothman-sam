from skimage import io as skio, measure, morphology, segmentation
from matplotlib import pyplot as plt
import numpy as np
from cellpose import models, core, io as cpio, plot
from pathlib import Path
from tqdm import trange
from segcellpose import readtiffs
import json
import os
import concurrent.futures
from itertools import repeat


def writetiffs(imgs: dict[str, np.ndarray], writedir: Path, ext: str) -> None:
	if not ext.startswith('.'):
		ext = '.' + ext

	for name, img in imgs.items():
		skio.imsave(writedir / f"{name}{ext}.tiff", img)

	return


def _filtereccentricity_img(img: np.ndarray, maxecc: float) -> np.ndarray:
	eccs = measure.regionprops_table(img, properties=['label', 'eccentricity'])
	eccmap = dict(zip(eccs['label'], eccs['eccentricity']))
	noisemask = np.zeros(img.shape, dtype='uint16')

	for objlabel, ecc in eccmap.items():
		if ecc > maxecc:
			noisemask += np.where(img == objlabel, objlabel, 0).astype('uint16')
	
	eccimg = img - noisemask

	return eccimg


def _filtereccentricity(indir: Path | dict[str, np.ndarray], eccdir: Path | None = None, maxecc: float = 0.95) -> None | dict[str, np.ndarray]:
	if isinstance(indir, Path):
		imgs = readtiffs(indir)
	else:
		imgs = indir
	
	eccimgs = {}
	for name, img in imgs.items():
		eccimgs[name] = _filtereccentricity_img(img, maxecc)

	if eccdir is not None:
		writetiffs(eccimgs, eccdir, ext=".ecc")
		return
	else:
		return eccimgs


def _getobjmap(img: np.ndarray) -> dict[int, np.ndarray]:
	labels = np.unique(img).astype(int)
	objmap = {}

	for label in labels:
		if label != 0:
			objmap[label] = np.where(img == label, 1, 0).astype('uint16')
	
	return objmap


# BROKEN: feel like the even numbers aren't finding themselves on the overlap?
def _findneighbors_img(img: np.ndarray, mintouchprop: float) -> dict[str, set]:
	boundarymap = {}
	objs = _getobjmap(img)

	for objlabel, mask in objs.items():
		boundarymap[objlabel] = segmentation.find_boundaries(mask, mode='inner').astype('uint16')
		boundarymap[objlabel] = morphology.binary_dilation(boundarymap[objlabel], morphology.disk(1))

	del objs

	intersections = {}
	# for label1, mask1 in boundarymap.items(): # how much label1 touches label2
	# 	intersections[label1] = {}
	# 	for label2, mask2 in boundarymap.items():
	# 		touchprop = (mask1 & mask2).sum() / mask1.sum()
	# 		if touchprop > mintouchprop:
	# 			intersections[label1][label2] = touchprop

	for label, mask in boundarymap.items():
		intersections[label] = {}
		overlapmask = (mask & img) * img
		print(f"{label}: {(mask&img).sum()}") # size of intersecting part
		overlaplabels = np.unique(overlapmask).astype(int)
		print(f"{label}: {np.unique(img)}: {overlaplabels}")

		for overlaplabel in overlaplabels:
			if overlaplabel != 0:
				touchprop = (overlapmask[overlapmask == overlaplabel].sum()) / mask.sum()
				intersections[label][overlaplabel] = touchprop


	return intersections


def _findneighbors(imgs: dict[str, np.ndarray], mintouchprop: float = 0, export: Path | bool = False, maxworkers=None):
	neighbormaps = {}

	for name, img in imgs.items():
		print(name)
		neighbormaps[name] = _findneighbors_img(img, mintouchprop)
	# with concurrent.futures.ProcessPoolExecutor(max_workers=maxworkers) as exec:
	# 	for name, neighbormap in zip(imgs.keys(), exec.map(_findneighbors_img, imgs.values(), repeat(mintouchprop))):
	# 		print(name)
	# 		neighbormaps[name] = neighbormap

	if export:
		with open(export, 'w') as f:
			json.dump(neighbormaps, f)
	else:
		return neighbormaps


def _mergeclusters():
	return


def _filtersize():
	return


def noisereduction(segdir: Path, filterdir: Path, maxecc: float=0.95, mintouchprop: float=0, sizebounds: tuple[int, int]=(7000,11000), fluor: int=300, export: dict[str, bool]=None, maxworkers: int=os.cpu_count()) -> None:
	eccdir = filterdir / "ecc"
	size1dir = filterdir / "size1"
	clusterdir = filterdir / "cluster"
	size2dir = filterdir / "size2"

	os.makedirs(eccdir, exist_ok=True)
	os.makedirs(size1dir, exist_ok=True)
	os.makedirs(clusterdir, exist_ok=True)
	os.makedirs(size2dir, exist_ok=True)

	sizebounds = tuple(sorted(sizebounds))

	imgs = readtiffs(segdir)
	
	eccimgs = _filtereccentricity(imgs, None, maxecc)
	neighbormap = _findneighbors(imgs, mintouchprop=mintouchprop, export=export['neighbors'], maxworkers=maxworkers)

	return


if __name__ == "__main__":
	FLUORDIR = Path("/Users/kvying/Files/ucsb-local/rothman-sam/fluorescence-pipeline")
	DATADIR = Path("/Users/kvying/Files/ucsb-local/rothman-data/embryo")
	RAWDIR = DATADIR / "11-13-2025"
	OUTDIR = RAWDIR / "out"
	SEGDIR = OUTDIR / "seg"
	FILTERDIR = OUTDIR / "filtered"
	
	ecc = 0.95
	mintouchprop = 0
	sizebounds = (7000, 11000)
	fluorthreshold = 300

	export = {'ecc': False, 'size1': False, 'cluster': False, 'size2': True, 'neighbors': FILTERDIR / "neighbors.json"}

	noisereduction(SEGDIR, FILTERDIR, maxecc=ecc, mintouchprop=mintouchprop, sizebounds=sizebounds, fluor=fluorthreshold, export=export)