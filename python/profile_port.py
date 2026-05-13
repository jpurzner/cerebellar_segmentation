"""Run pipeline with per-stage timing to find the bottleneck."""
import time
import numpy as np
import scipy.ndimage as ndi
import tifffile
from skimage import exposure, filters, morphology
from skimage.morphology import disk

INPUT = "/Users/jpurzner/Dropbox/images/edu_repeat/p27/2018_05_22_s1_2_p27-0021/2018_05_22_s1_2_p27-0021_fused_crop.tif"
PX_UM = 0.5119049

stack = tifffile.imread(INPUT)
print(f"shape={stack.shape}")
p27_raw  = stack[0]
neun_raw = stack[1]
dapi_raw = stack[2]
H, W = p27_raw.shape

def t(name):
    if not hasattr(t, "last"): t.last = time.time(); print(f"--- start ---"); return
    now = time.time()
    print(f"  {now - t.last:6.1f}s  {name}")
    t.last = now

def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)

t("init")

a_nf_pre = exposure.equalize_adapthist(to_unit(p27_raw),  clip_limit=0.01)
t("CLAHE p27")
b_nf_pre = exposure.equalize_adapthist(to_unit(dapi_raw), clip_limit=0.01)
t("CLAHE dapi")
c_nf_pre = exposure.equalize_adapthist(to_unit(neun_raw), clip_limit=0.01)
t("CLAHE neun")

# gaussian blurs
sig = 5.0
a_nfg = filters.gaussian(a_nf_pre, sigma=sig, preserve_range=True)
t(f"gaussian sigma={sig}")
b_nfg = filters.gaussian(b_nf_pre, sigma=sig, preserve_range=True)
t(f"gaussian sigma={sig}")
c_nfg = filters.gaussian(c_nf_pre, sigma=sig, preserve_range=True)
t(f"gaussian sigma={sig}")

# scipy.ndimage equivalents
ag2 = ndi.gaussian_filter(a_nf_pre.astype(np.float32), sigma=sig)
t(f"ndi.gaussian sigma={sig}")
ag3 = ndi.gaussian_filter(a_nf_pre.astype(np.float32), sigma=30.0)
t(f"ndi.gaussian sigma=30")
ag4 = ndi.gaussian_filter(a_nf_pre.astype(np.float32), sigma=63.0)
t(f"ndi.gaussian sigma=63 (was used for adapthresh)")

# morphological ops with various disk sizes
mask = a_nf_pre > 0.5
for r in [5, 10, 20, 40]:
    closed = morphology.closing(mask, disk(r))
    t(f"morph.closing disk({r})")

# imfill
fill = ndi.binary_fill_holes(closed)
t("ndi.binary_fill_holes")

# skeletonize
sk = morphology.skeletonize(closed)
t("skeletonize")

# rank filter (used for ordfilt2)
m = ndi.maximum_filter(a_nf_pre, size=5)
t("ndi.maximum_filter size=5")

# fft convolve
from scipy.signal import fftconvolve
ks = np.ones((50, 50), dtype=np.float32)
out = fftconvolve(a_nf_pre.astype(np.float32), ks, mode='same')
t("fftconvolve 50x50 kernel")
