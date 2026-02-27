import pandas as pd
import os
import requests
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor


def _gettranscript(iso: str) -> str:
	transcript = iso.split('.')[:-1] # all except last '.'
	if len(transcript) > 1:
		return '.'.join(transcript)
	else:
		return iso


def _getisomap(isos: list) -> dict[str, list]:
	isomap = {}

	for iso in isos:
		transcript = _gettranscript(iso)
		if transcript in isomap:
			isomap[transcript].append(iso)
		else:
			isomap[transcript] = [iso]
	
	return {transcript: isos for transcript, isos in isomap.items() if len(isos) > 1}


# def _queryiso(iso: str, session: requests.Session) -> bool:
# 	url = f"http://rest.wormbase.org/rest/widget/transcript/{iso}/overview"
# 	with session.get(url) as response:
# 		query = response.json()
# 	return 'reason' in query

def _queryiso(iso: str, session: requests.Session) -> bool:
    url = f"http://rest.wormbase.org/rest/widget/transcript/{iso}/overview"
    try:
        with session.get(url, timeout=10) as response:
            if not response.ok:
                return False
            
            content_type = response.headers.get('Content-Type', '')
            if 'application/json' not in content_type:
                return False

            query = response.json()
            return 'reason' in query
            
    except (requests.exceptions.RequestException, ValueError):
        return False


def _mpqueryiso(isos: list, session: requests.Session) -> dict[str, bool]:
	# n=1000:  14.451182s
	# n=2000:  19.883863s
	# n=5000:  83.820977s
	# n=10000: 153.033952s
	# n=20000: 201.934475s
	# n=50000: 567.078500s
	# n=len(transcripts): 634.591072s
	with ThreadPoolExecutor(max_workers=16) as exec:
		res = list(exec.map(lambda iso: _queryiso(iso, session), isos))
	
	return {iso:val for iso, val in zip(isos, res)}


def _makemergemap(isomap: dict[str, list], isoquery: dict[str, bool]) -> dict[str, list]:
	mergemap = {}

	for transcript, isos in isomap.items():
		bools = list(map(isoquery.get, isos))
		if any(bools):
			mergemap[transcript] = isos

	return mergemap


def _mergeisos(rawcounts: pd.DataFrame, mergemap: dict[str, list]) -> pd.DataFrame:
    invmap = {iso: transcript for transcript, isos in mergemap.items() for iso in isos}
    rawcounts['root'] = rawcounts['transcript'].map(invmap).fillna(rawcounts['transcript'])

    numeric_cols = rawcounts.select_dtypes(include=['number']).columns
    
    rawcounts = rawcounts.groupby('root')[numeric_cols].sum().reset_index()
    rawcounts[numeric_cols] = rawcounts[numeric_cols].astype(int)

    rawcounts.rename(columns={'root': 'transcript'}, inplace=True)

    return rawcounts


def resolveiso(tsvpath: Path) -> pd.DataFrame:
	print("Reading raw counts...")
	rawcounts = pd.read_csv(tsvpath, sep='\t')
	isos = list(rawcounts['transcript'])

	print("Organizing isoforms by transcript...")
	isomap = _getisomap(isos)
	
	print("Validating isoforms...")
	session = requests.Session()
	start = time.perf_counter()
	isoquery = _mpqueryiso(isos, session)
	end = time.perf_counter()
	print(f"{(end - start):.6f} seconds elapsed. ")
	
	print("Filtering for merge candidates...")
	mergemap = _makemergemap(isomap, isoquery)
	
	print("Merging candidates...")
	rawcounts = _mergeisos(rawcounts, mergemap)

	print("Done!")

	return rawcounts


if __name__ == "__main__":
	DIR = Path("/home/hanwenying/rothman-sam/deg-pipeline")
	OUTDIR = DIR / "out"
	os.makedirs(OUTDIR, exist_ok=True)

	resolved_rawcounts = resolveiso(DIR/'raw_rawcounts.tsv')
	print(resolved_rawcounts.head())
	resolved_rawcounts.to_csv(OUTDIR/'resolved_rawcounts.tsv', sep='\t')
