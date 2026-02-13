from skimage import io as skio, measure, morphology
from matplotlib import pyplot as plt
import numpy as np
from cellpose import models, core, io as cpio, plot
from pathlib import Path
from tqdm import trange
from segcellpose import readtiffs
import json
import os


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


def _findneighbors(): # outputs json of neighbors
	return


def _mergeclusters():
	return


def _filtersize():
	return


def noisereduction(segdir: Path, filterdir: Path, maxecc: float=0.95, touchprop: float=0, sizebounds: tuple[int, int]=(7000,11000), fluor: int=300) -> None:
	eccdir = filterdir / "ecc"
	size1dir = filterdir / "size1"
	clusterdir = filterdir / "cluster"
	size2dir = filterdir / "size2"

	os.makedirs(eccdir, exist_ok=True)
	os.makedirs(size1dir, exist_ok=True)
	os.makedirs(clusterdir, exist_ok=True)
	os.makedirs(size2dir, exist_ok=True)

	sizebounds = tuple(sorted(sizebounds))

	_filtereccentricity(segdir, eccdir, maxecc)

	return


if __name__ == "__main__":
	FLUORDIR = Path("/Users/kvying/Files/ucsb-local/rothman-sam/fluorescence-pipeline")
	DATADIR = Path("/Users/kvying/Files/ucsb-local/rothman-data/embryo")
	RAWDIR = DATADIR / "11-13-2025"
	OUTDIR = RAWDIR / "out"
	SEGDIR = OUTDIR / "seg"
	FILTERDIR = OUTDIR / "filtered"
	
	ecc = 0.95
	touchprop = 0
	sizebounds = (7000, 11000)
	fluorthreshold = 300

	noisereduction(SEGDIR, FILTERDIR, maxecc=ecc, touchprop=touchprop, sizebounds=sizebounds, fluor=fluorthreshold)