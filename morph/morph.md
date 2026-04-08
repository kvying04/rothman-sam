---
title: "Proposed Morphology Pipeline"
author: "Kevin Ying"
date: "Last Updated 2026-04-01"
---
# 0. Imaging Parameters
- 100µm resolution; <!-- worthwhile to go down to 50µm? binary masking is struggling to retain donut characteristics -->
- 3-second interval between frames; <!-- pref. ~1.5s intervals; do a run with this, see how diff. it is? -->
- 40-frame "stacks".

# 1. Quantifying "Donut-ness"
Initial approach: 

1. segment image into potential donuts;
2. Find fluorescence profile of each individual candidate and discriminate "blob" vs. "donut".

## 1.1. Denoised Segmentation of "Mitochondrial Objects"
### 1.1.1. Important Parameters
Basic filtering/denoising parameters need to be determined (similarly to Cellpose case). <!-- actually worthwile trying cellpose? -->

- Maximum eccentricity (`maxecc`; how "round" an object needs to be);
- Minimum size (`minsize`; how large an object needs to be).

### 1.1.2. Watershedding
Likely will need to watershed to deal with candidate units engaging in fusion/fission events. 

During initial segmentation, is worthwhile to watershed for candidates (assumed to be merged) that are $\text{size} \le 1.5 \cdot \text{minsize}$ and $\text{ecc} > 1.7 \cdot \text{maxecc}$.[^1]

## 1.2. Fluorescence Profile
For each object, generate cross-section ROI line across both major axes (perpendicular always). Then, find fluorescence profile for each axis. If a fluorescent "divot" is observed for both plots, the object can be considered a donut. 

Need to determine the length and magnitude of the divot to make sure that drops in fluorescence are not just noise. 

A rigorous definition for differentiation: 

- Uniformity of "ring" width:[^2]
  - Find the width of both sides of the fluorescing ring for both ROIs (4 total measurements);
  - For each ROI, $w_1 = w_2 \pm n_a$ for some leniency $n_a$.
  - Compare ROIs by looking at $\max(w_1,w_2)_{\text{axis=0}} = \max(w_1, w_2)_{\text{axis=1}} \pm n_b$ for some other leniency $n_b$.
- Substantiality of divot:
  - Ensure that, for each ROI, $\frac{w_1 + w_2}{(w_1 + w_2) + w_{\text{divot}}} < 0.5$.[^3]
- Fluorescence difference of donut:
  - Want to see that $\min(\max({\text{fluor}_{l_1}}), \max({\text{fluor}_{l_2}})) > \max(\text{fluor}_{\text{divot}})$ for both axes by some margin. 

## 1.3. Adjusting for Time-Series Changes


<!-- 
# 2. Object-wise Participation in Network Events
## 2.1. 

# 3. Representing Donut Mitochondria Graphically -->

[^1]: Coefficients subject to change. 
[^2]: Extremely skewed donuts have not been observed (at least by me). 
[^3]: Based on heuristic that the donut ring thickness is never more than half of the diameter of the whole object. 