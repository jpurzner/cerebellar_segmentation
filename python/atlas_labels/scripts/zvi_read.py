"""
Minimal reader for Zeiss AxioVision .zvi files (OLE2 compound files), no extra packages.
Returns each image plane (Image/Item(n)/Contents) as a 2-D array, plus the UTF-16 strings of the file.

    from zvi_read import read_zvi
    planes, info = read_zvi(path)       # planes: list of (H, W) uint16 arrays
"""
import re
import struct

import numpy as np


class CFB:
    def __init__(self, path):
        self.b = b = open(path, "rb").read()
        assert b[:8] == bytes.fromhex("D0CF11E0A1B11AE1"), "not an OLE2 file"
        self.ss = 1 << struct.unpack_from("<H", b, 0x1E)[0]
        self.mss = 1 << struct.unpack_from("<H", b, 0x20)[0]
        nfat, self.dir0 = struct.unpack_from("<II", b, 0x2C)
        self.cut, mfat0, nmfat, difat0, ndifat = struct.unpack_from("<IIIII", b, 0x38)
        difat = list(struct.unpack_from("<109I", b, 0x4C))
        s = difat0
        for _ in range(ndifat):
            v = struct.unpack_from(f"<{self.ss // 4}I", b, self._off(s))
            difat += list(v[:-1])
            s = v[-1]
        self.fat = []
        for s in difat[:nfat]:
            self.fat += struct.unpack_from(f"<{self.ss // 4}I", b, self._off(s))
        self.dirs = self._dirs()
        root = self.dirs[0]
        self.ministream = self._chain(root["start"], root["size"])
        self.minifat = list(struct.unpack(f"<{len(m) // 4}I", m)) if (m := self._chain(mfat0, nmfat * self.ss)) else []

    def _off(self, s):
        return (s + 1) * self.ss

    def _chain(self, s, size):
        out = bytearray()
        while s < 0xFFFFFFFA and len(out) < size:
            out += self.b[self._off(s):self._off(s) + self.ss]
            s = self.fat[s]
        return bytes(out[:size])

    def _minichain(self, s, size):
        out = bytearray()
        while s < 0xFFFFFFFA and len(out) < size:
            out += self.ministream[s * self.mss:(s + 1) * self.mss]
            s = self.minifat[s]
        return bytes(out[:size])

    def _dirs(self):
        raw = self._chain(self.dir0, 1 << 30)
        out = []
        for i in range(len(raw) // 128):
            e = raw[i * 128:(i + 1) * 128]
            nl = struct.unpack_from("<H", e, 64)[0]
            name = e[:max(nl - 2, 0)].decode("utf-16le", "ignore")
            typ = e[66]
            left, right, child = struct.unpack_from("<iii", e, 68)
            start, size = struct.unpack_from("<IQ", e, 116)
            out.append(dict(name=name, type=typ, left=left, right=right, child=child, start=start, size=size & 0xFFFFFFFF))
        return out

    def stream(self, i):
        d = self.dirs[i]
        return self._minichain(d["start"], d["size"]) if d["size"] < self.cut else self._chain(d["start"], d["size"])

    def paths(self):
        """(path, dir index) for every stream."""
        res = []

        def walk(i, pre):
            if i < 0 or i >= len(self.dirs):
                return
            d = self.dirs[i]
            walk(d["left"], pre)
            walk(d["right"], pre)
            p = pre + "/" + d["name"]
            if d["type"] == 2:
                res.append((p, i))
            elif d["type"] in (1, 5):
                walk(d["child"], p)
        walk(self.dirs[0]["child"], "")
        return res


def read_zvi(path, width=None, height=None):
    c = CFB(path)
    strings = sorted(set(x.decode("utf-16le") for x in re.findall(rb"(?:[\x20-\x7e]\x00){4,}", c.b)))
    planes = []
    for p, i in sorted(c.paths()):
        if re.fullmatch(r"/Image/Item\(\d+\)/Contents", p):
            s = c.stream(i)
            # the Contents header holds width, height, depth, pixel format; pixels are the stream's tail
            w, h = width, height
            if w is None:
                for off in range(0, 200, 2):
                    try:
                        a, bb = struct.unpack_from("<ii", s, off)
                    except struct.error:
                        break
                    if 100 < a < 20000 and 100 < bb < 20000 and a * bb * 2 <= len(s) and len(s) - a * bb * 2 < 8192:
                        w, h = a, bb
                        break
            n = w * h * 2
            planes.append((p, np.frombuffer(s[-n:], "<u2").reshape(h, w)))
    return planes, dict(strings=strings, streams=[p for p, _ in c.paths()])


if __name__ == "__main__":
    import sys
    pl, info = read_zvi(sys.argv[1])
    print(len(info["streams"]), "streams;", [p for p in info["streams"] if "Contents" in p][:6])
    for p, a in pl:
        print(p, a.shape, a.dtype, a.min(), a.max(), round(float(a.mean()), 1))
