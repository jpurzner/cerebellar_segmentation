"""
Replay the user's paint / erase strokes from the New IGL Step page onto the 4 um grid, exactly as
the page draws them (cell (i, j) is inside a stroke when its centre (i + 0.5, j + 0.5) lies within
r_um / 4 cells of one of the stroke's segments; later strokes override earlier ones), and give
the section's corrected IGL: the proposed IGL (MATLAB's IGL plus the new step's additions)
with the painted cells added and the erased cells removed.

    from igl_detect_edits import replay, corrected_igl
"""
import numpy as np


def replay(strokes, H, W):
    """user layer: 0 untouched, 1 painted, 2 erased."""
    user = np.zeros((H, W), np.uint8)
    yy, xx = np.mgrid[0:H, 0:W]
    cy, cx = yy + 0.5, xx + 0.5
    for st in strokes or []:
        val = 1 if st["mode"] == "paint" else 2
        r = st["r_um"] / 4.0
        P = [(p[0] * W, p[1] * H) for p in st["pts"]]
        segs = [(P[0], P[0])] if len(P) == 1 else list(zip(P[:-1], P[1:]))
        for (ax, ay), (bx, by) in segs:
            x0, x1 = max(int(np.floor(min(ax, bx) - r - 1)), 0), min(int(np.ceil(max(ax, bx) + r + 1)), W - 1)
            y0, y1 = max(int(np.floor(min(ay, by) - r - 1)), 0), min(int(np.ceil(max(ay, by) + r + 1)), H - 1)
            if x1 < x0 or y1 < y0:
                continue
            sx, sy = cx[y0:y1 + 1, x0:x1 + 1], cy[y0:y1 + 1, x0:x1 + 1]
            dx, dy = bx - ax, by - ay
            L2 = dx * dx + dy * dy
            t = np.clip(((sx - ax) * dx + (sy - ay) * dy) / L2, 0, 1) if L2 > 0 else np.zeros_like(sx)
            ex, ey = sx - (ax + t * dx), sy - (ay + t * dy)
            m = ex * ex + ey * ey <= r * r
            user[y0:y1 + 1, x0:x1 + 1][m] = val
    return user


def corrected_igl(proposed, strokes):
    """proposed: grid IGL proposal (bool); returns the IGL after the user's edits."""
    user = replay(strokes, *proposed.shape)
    return (proposed & (user != 2)) | (user == 1)
