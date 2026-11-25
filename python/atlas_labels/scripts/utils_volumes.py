"""Loaders + axis conventions for LAMBADA + Allen volumes.

Empirically determined for the P14 volumes:
  LAMBADA P14 ref (NIfTI, RAS, identity affine, 25 µm iso)
      numpy shape (382, 480, 328) = (R-L, A-P, S-I)
  ADMBA   P14 ref (MetaImage, .mhd labels itself RAI but stored Allen-style)
      numpy shape (390, 469, 794) = (D-V, L-R, A-P)
      spacing (25, 16.752, 16.752) µm along those numpy axes
"""
import json
from pathlib import Path
import numpy as np
import nibabel as nib
import SimpleITK as sitk

ROOT = Path('/Users/jpurzner/Dropbox/Q_research/py_projects/cerebellar_atlas_labels')

LAMBADA_AXES = {'sagittal': 0, 'coronal': 1, 'axial': 2}    # numpy axis to slice for each view
LAMBADA_SPACING_UM = (25.0, 25.0, 25.0)                     # along numpy axes (R-L, A-P, S-I)

ADMBA_AXES = {'axial': 0, 'sagittal': 1, 'coronal': 2}      # numpy axis to slice
ADMBA_SPACING_UM = (25.0, 16.752, 16.752)                   # along numpy axes (D-V, L-R, A-P)


def load_lambada_ref(age='P14'):
    p = ROOT/f'LAMBADA/LAMBADA-{age}/LAMBADA_25um_reference_{age}_v1.0.nii.gz'
    img = nib.load(p)
    return img.get_fdata(), img


def load_admba(age='P14', what='atlas'):
    """what in {'atlas','annotation'}"""
    sub = 'atlasVolume/atlasVolume/atlasVolume.mhd' if what == 'atlas' else 'annotation/annotation.mhd'
    p = ROOT/f'allen_reference/ADMBA_{age}'/sub
    img = sitk.ReadImage(str(p))
    return sitk.GetArrayFromImage(img), img


def admba_cerebellar_ids(structure_graph_path=None):
    """Return set of label ids whose structure (or any ancestor) name contains 'cerebell' in ADMBA graph 17."""
    if structure_graph_path is None:
        structure_graph_path = ROOT/'allen_reference/ADMBA_structure_graph_17.json'
    data = json.loads(Path(structure_graph_path).read_text())
    cb = set()
    def walk(node, inside_cb):
        name = (node.get('name') or '').lower()
        acr  = (node.get('acronym') or '')
        is_cb = inside_cb or ('cerebell' in name) or (acr.upper() in {'CB','CBX','CBN','VERM','HEM','CbV','CbH','CbVCx'})
        if is_cb:
            cb.add(node['id'])
        for c in node.get('children', []) or []:
            walk(c, is_cb)
    for top in data.get('msg', []):
        walk(top, False)
    return cb


def to_display_sagittal_lambada(vol2d):
    """LAMBADA sagittal slice lam[mid_x,:,:] has shape (A-P, S-I).
       Display: anterior right, superior up. Use .T + origin='lower':
         rows = S-I (axis 1), cols = A-P (axis 0)."""
    return vol2d.T   # used with origin='lower'


def to_display_sagittal_admba(vol2d):
    """ADMBA sagittal slice adm[:,mid_y,:] has shape (D-V, A-P).
       D-V (axis 0) increases inferiorly → flip so that 'up' on screen = superior.
       For anterior-right we need to know the file's A-P direction; we'll
       choose orientation later visually."""
    return vol2d   # display caller applies origin/flip
