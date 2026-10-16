#!/usr/bin/env python3
"""make_logo.py -- write the SLEvolve logo, favicon and a 32 px preview strip.

Letters are Raleway (SIL OFL 1.1; default ~/Library/Fonts/PhotoDirector/Raleway.ttf, or set
SLEVOLVE_LOGO_FONT) outlines converted to paths with matplotlib's TextPath, so the SVGs need no font.
Each glyph also gets a light outline stroke of the same colour (STROKE font units) so the regular
Raleway outlines read as semibold.  If Raleway is not found the script falls back to Inter
(tools/fonts/Inter-{SemiBold,Bold}.ttf, kept for a possible later weight change; no extra stroke).  Glyph positions (advances + the font's kerning) come from matplotlib's
text_to_path layout of the whole word; each letter is then coloured separately.

  wordmark: 'SL' teal, 'E' amber, 'volve' grey, all the same weight; below it a single teal
  sea-level wave placed AFTER the word (right of 'volve'), centred on the x-height, decaying to the
  right (amplitude and wavelength shrink).

Outputs (in OUTDIR): logo.svg, logo.png (1600 px wide, transparent), logo_dark.png (400 px wide on
#111A22), logo_small.png (320 px wide, transparent), favicon.svg, favicon_{64,128,512}.png,
favicon.png (=512), logo_preview.png (wordmark 32 px tall on white and on dark).

    python tools/make_logo.py docs/src/img [raleway|semibold|bold]      (needs numpy and matplotlib)
"""
import sys, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath, text_to_path
from matplotlib.font_manager import FontProperties, get_font
from matplotlib.path import Path
from matplotlib.patches import PathPatch, FancyBboxPatch
from matplotlib.transforms import Affine2D

OUT = sys.argv[1] if len(sys.argv) > 1 else "docs/src/img"
WEIGHT = sys.argv[2] if len(sys.argv) > 2 else "raleway"
FONTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts")
RALEWAY = os.environ.get("SLEVOLVE_LOGO_FONT", os.path.expanduser("~/Library/Fonts/PhotoDirector/Raleway.ttf"))
if WEIGHT == "raleway" and not os.path.exists(RALEWAY):
    print("Raleway not found; falling back to Inter SemiBold"); WEIGHT = "semibold"
FONTFILE = RALEWAY if WEIGHT == "raleway" else os.path.join(
    FONTDIR, {"semibold": "Inter-SemiBold.ttf", "bold": "Inter-Bold.ttf"}[WEIGHT])
STROKE = 5.0 if WEIGHT == "raleway" else 0.0   # same-colour outline (font units), same for every letter
FP = FontProperties(fname=FONTFILE)
TEAL, ACCENT, GREY = "#1E88A8", "#F29E38", "#8794A1"
DARK = "#111A22"
TRACK = 1.0          # extra letter gap (font units at size 100) on top of advance + kerning
SIZE = 100

def word(text, colours, track=TRACK):
    """glyph paths placed by the font's advances and kerning; returns [(path, colour)]"""
    # pen positions from FreeType advances + kerning (do not rely on get_glyphs_with_font, whose
    # x positions are all 0 in matplotlib 3.11, which stacked every glyph at x=0)
    from matplotlib.ft2font import LoadFlags, Kerning
    font = get_font(FONTFILE)
    font.set_size(text_to_path.FONT_SCALE, text_to_path.DPI)
    pen, prev, xs = 0.0, None, []
    for ch in text:
        gi = font.get_char_index(ord(ch))
        if prev is not None:
            pen += font.get_kerning(prev, gi, Kerning.DEFAULT) / 64
        xs.append(pen)
        pen += font.load_char(ord(ch), flags=LoadFlags.NO_HINTING).linearHoriAdvance / 65536
        prev = gi
    items = []
    for i, (ch, c, x0) in enumerate(zip(text, colours, xs)):
        x = x0 * SIZE / text_to_path.FONT_SCALE + i * track
        tp = TextPath((x, 0), ch, size=SIZE, prop=FP)
        # plain Path with the offset baked into the vertices (no lazy TextPath offset)
        items.append((Path(np.array(tp.vertices), np.array(tp.codes)), c))
    return items

def ink(items):
    ex = [p.get_extents() for p, _ in items]
    return min(e.x0 for e in ex), max(e.x1 for e in ex), min(e.y0 for e in ex), max(e.y1 for e in ex)

def wave(x0, x1, ymid, amp0, n=600):
    """decaying wave: amplitude ~ exp(-1.8 t), phase accelerates so the wavelength shrinks; ~3.5 crests"""
    t = np.linspace(0, 1, n)
    phase = 2 * np.pi * (1.8 * t + 1.7 * t ** 2)        # 3.5 cycles total, local wavelength shrinking
    y = ymid + amp0 * np.exp(-1.8 * t) * np.sin(phase)
    return Path(np.c_[x0 + (x1 - x0) * t, y])

def wordmark():
    items = word("SLEvolve", [TEAL, TEAL, ACCENT] + [GREY] * 5)
    x0, x1, y0, y1 = ink(items)
    amp, gap, length, lw = 20.0, 9.0, 108.0, 5.0          # font units at size 100
    xs = x1 + STROKE / 2 + gap + lw / 2
    w = wave(xs, xs + length, 26.0, amp)                   # centred ~ x-height middle
    s2 = STROKE / 2
    return items, w, lw, (x0 - s2, xs + length + lw / 2, min(y0 - s2, 26 - amp - lw), y1 + s2)

def render(fname, width_px=None, height_px=None, bg=None, pad=0.03):
    items, w, lw, (x0, x1, y0, y1) = wordmark()
    W, H = x1 - x0, y1 - y0
    p = pad * W
    Wt, Ht = W + 2 * p, H + 2 * p
    if width_px is None:
        width_px = height_px * Wt / Ht
    dpi = 100
    fig = plt.figure(figsize=(width_px / dpi, width_px / dpi * Ht / Wt), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()
    ax.set_xlim(x0 - p, x1 + p); ax.set_ylim(y0 - p, y1 + p); ax.set_aspect("equal")
    if bg:
        fig.patch.set_facecolor(bg)
    pts_per_unit = width_px / Wt * 72 / dpi
    for path, c in items:
        ax.add_patch(PathPatch(path, fc=c, ec=c if STROKE else "none", lw=STROKE * pts_per_unit,
                               joinstyle="round"))
    ax.add_patch(PathPatch(w, fc="none", ec=TEAL, lw=lw * pts_per_unit, capstyle="round", joinstyle="round"))
    fig.savefig(fname, transparent=bg is None, facecolor=bg if bg else "none", dpi=dpi,
                metadata={"Date": None} if fname.endswith(".svg") else None)
    plt.close(fig)

def favicon(fname, px):
    fig = plt.figure(figsize=(px / 100, px / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.add_patch(FancyBboxPatch((0.02, 0.02), 0.96, 0.96, boxstyle="round,pad=0,rounding_size=0.18",
                                fc=TEAL, ec="none"))
    items = word("SLE", ["white"] * 3, track=0.0)
    a, b, c, d = ink(items)
    s = 0.74 / (b - a + STROKE)
    for path, col in items:
        v = (path.vertices - [(a + b) / 2, (c + d) / 2]) * s + 0.5
        ax.add_patch(PathPatch(Path(v, path.codes), fc=col, ec=col if STROKE else "none",
                               lw=STROKE * s * px * 0.72, joinstyle="round"))
    fig.savefig(fname, transparent=True, dpi=100,
                metadata={"Date": None} if fname.endswith(".svg") else None)
    plt.close(fig)

def preview(fname):
    import matplotlib.image as mpimg, tempfile
    tmp = tempfile.mkdtemp()
    a, b = os.path.join(tmp, "w.png"), os.path.join(tmp, "d.png")
    render(a, height_px=32, bg="white", pad=0.02); render(b, height_px=32, bg=DARK, pad=0.02)
    A, B = mpimg.imread(a)[..., :3], mpimg.imread(b)[..., :3]
    h, wd = max(A.shape[0], B.shape[0]), A.shape[1] + B.shape[1] + 16
    S = np.ones((h + 16, wd + 16, 3)); S[:, (wd + 16) // 2:] = np.array([0x11, 0x1A, 0x22]) / 255
    S[8:8 + A.shape[0], 8:8 + A.shape[1]] = A
    S[8:8 + B.shape[0], 16 + A.shape[1] + 8 - 8:16 + A.shape[1] + B.shape[1]] = B
    mpimg.imsave(fname, S)

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    plt.rcParams["svg.fonttype"] = "path"
    render(os.path.join(OUT, "logo.svg"), width_px=640)
    render(os.path.join(OUT, "logo.png"), width_px=1600)
    render(os.path.join(OUT, "logo_small.png"), width_px=320)
    render(os.path.join(OUT, "logo_dark.png"), width_px=400, bg=DARK, pad=0.06)
    favicon(os.path.join(OUT, "favicon.svg"), 64)
    for n in (64, 128, 512):
        favicon(os.path.join(OUT, f"favicon_{n}.png"), n)
    favicon(os.path.join(OUT, "favicon.png"), 512)
    preview(os.path.join(OUT, "logo_preview.png"))
    print("wrote logo files to", OUT, "with", os.path.basename(FONTFILE))
