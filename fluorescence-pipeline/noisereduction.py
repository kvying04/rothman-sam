from skimage import io as skio, measure, morphology
from matplotlib import pyplot as plt
import numpy as np
from cellpose import models, core, io as cpio, plot
from pathlib import Path
from tqdm import trange
import os
import gc


def _denoise():
	return


def _filtereccentricity():
	return


def _findneighbors():
	return


def _mergeclusters():
	return


def _filtersize():
	return


def noisereduction():
	return


if __name__ == "__main__":
	FLUORDIR = Path("/Users/kvying/Files/ucsb-local/rothman-sam/fluorescence-pipeline")
	DATADIR = Path("/Users/kvying/Files/ucsb-local/rothman-data/embryo")
	RAWDIR = DATADIR / "11-13-2025"
	OUTDIR = RAWDIR / "out"
	MASKDIR = OUTDIR / "filtered"
	
	ecc = 0.95