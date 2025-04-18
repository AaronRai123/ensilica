"""visualization.py – clean, professional plasmid + linear map renderers
==================================================================
* **render_plasmid_map()** – returns base‑64 SVG of a circular map.
* **visualize_linear_map()** – returns base‑64 PNG of a linear map.

Improvements vs previous version
--------------------------------
* Unified, Google‑material colour palette via `FEATURE_COLOURS`.
* Smart lane algorithm for linear view (unbounded lanes, spills up/down).
* Text labels auto‑flip for readability (no upside‑down text).
* Outer tick marks every 1 kb on plasmid (configurable).
* High‑DPI output (`SVG 120 dpi`, `PNG 300 dpi`).
* Utilities: `_arc_path`, `_arrow_polygon`, `_text_rotate` for clarity.

Dependencies
------------
* `biopython>=1.80`
* `numpy`, `matplotlib` (linear only)

Usage
-----
```python
svg_b64 = render_plasmid_map(record, highlight_region={"start":500,"end":800})
png_b64 = visualize_linear_map(record)
```
"""
from __future__ import annotations

import base64
import io
import math
from typing import Dict, Any, Iterable, List, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
DPI_SVG = 120
DPI_PNG = 300
PLASMID_RADIUS = 250
PLASMID_WIDTH = 30
TICK_INTERVAL = 1000  # bp

FEATURE_COLOURS = {
    "CDS": "#4285F4", "gene": "#4285F4",
    "promoter": "#9C27B0", "terminator": "#E91E63",
    "rep_origin": "#2196F3", "origin": "#2196F3",
    "misc_feature": "#607D8B", "primer_bind": "#FF5722",
    "RBS": "#00BCD4", "regulatory": "#9C27B0",
    "resistance": "#F44336", "marker": "#FF9800",
}
DEFAULT_COLOUR = "#BDBDBD"
FONT_FAMILY = "Arial, Helvetica, sans‑serif"

# ---------------------------------------------------------------------------
#   PUBLIC API
# ---------------------------------------------------------------------------

def visualize_construct(record: SeqRecord, highlight_region: Dict[str, Any] | None = None) -> str:
    """
    Wrapper function for render_plasmid_map for backwards compatibility
    
    Args:
        record: SeqRecord object containing the construct
        highlight_region: Optional region to highlight
        
    Returns:
        Base64 encoded SVG data
    """
    return render_plasmid_map(record, highlight_region)

def render_plasmid_map(record: SeqRecord, highlight_region: Dict[str, Any] | None = None,
                       width: int = 800, height: int = 800) -> str:
    """Return **base‑64 SVG** circular plasmid map."""
    seq_len = len(record.seq)
    cx, cy = width / 2, height / 2

    def _arc_path(rad_out: float, rad_in: float, ang0: float, ang1: float) -> str:
        # ensure ang1 > ang0
        if ang1 < ang0:
            ang1 += 2 * math.pi
        large = int((ang1 - ang0) > math.pi)
        # outer arc endpoints
        x1, y1 = cx + rad_out * math.cos(ang0), cy + rad_out * math.sin(ang0)
        x2, y2 = cx + rad_out * math.cos(ang1), cy + rad_out * math.sin(ang1)
        # inner arc endpoints
        x3, y3 = cx + rad_in * math.cos(ang1), cy + rad_in * math.sin(ang1)
        x4, y4 = cx + rad_in * math.cos(ang0), cy + rad_in * math.sin(ang0)
        return (
            f"M {x1:.2f} {y1:.2f} "
            f"A {rad_out:.2f} {rad_out:.2f} 0 {large} 1 {x2:.2f} {y2:.2f} "
            f"L {x3:.2f} {y3:.2f} "
            f"A {rad_in:.2f} {rad_in:.2f} 0 {large} 0 {x4:.2f} {y4:.2f} Z"
        )

    def _arrow_polygon(angle: float, size: float = 14) -> str:
        # draw small triangular arrow centred on path at given angle
        ax, ay = cx + (PLASMID_RADIUS) * math.cos(angle), cy + (PLASMID_RADIUS) * math.sin(angle)
        pa = angle + math.pi / 2
        p1 = (ax + size * math.cos(pa), ay + size * math.sin(pa))
        p2 = (ax + size * math.cos(angle), ay + size * math.sin(angle))
        p3 = (ax - size * math.cos(pa), ay - size * math.sin(pa))
        return " ".join(f"{x:.2f},{y:.2f}" for x, y in (p1, p2, p3))

    def _text_rotate(angle: float) -> float:
        deg = math.degrees(angle) % 360
        return deg + 180 if 90 < deg < 270 else deg

    svg_parts: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"'
        f' style="font-family:{FONT_FAMILY}">',
        f'<rect width="100%" height="100%" fill="white"/>',
        # backbone
        f'<circle cx="{cx}" cy="{cy}" r="{PLASMID_RADIUS}" fill="none" stroke="#666" stroke-width="12"/>'
    ]

    # ticks every 1 kb
    for bp in range(0, seq_len, TICK_INTERVAL):
        ang = 2 * math.pi * bp / seq_len
        x1, y1 = cx + (PLASMID_RADIUS + 6) * math.cos(ang), cy + (PLASMID_RADIUS + 6) * math.sin(ang)
        x2, y2 = cx + (PLASMID_RADIUS + 18) * math.cos(ang), cy + (PLASMID_RADIUS + 18) * math.sin(ang)
        svg_parts.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#888"/>' )
        tx, ty = cx + (PLASMID_RADIUS + 28) * math.cos(ang), cy + (PLASMID_RADIUS + 28) * math.sin(ang)
        svg_parts.append(
            f'<text x="{tx:.1f}" y="{ty:.1f}" font-size="10" text-anchor="middle"'
            f' dominant-baseline="middle" transform="rotate({_text_rotate(ang):.1f} {tx:.1f} {ty:.1f})">{bp}</text>'
        )

    # FEATURE ARCS -----------------------------------------------------------
    for idx, feat in enumerate(record.features):
        if feat.type == "source":
            continue
        start, end = int(feat.location.start), int(feat.location.end)
        if end == start:
            continue
        colour = FEATURE_COLOURS.get(feat.type, DEFAULT_COLOUR)
        ang0, ang1 = (2 * math.pi * start / seq_len), (2 * math.pi * end / seq_len)
        path = _arc_path(PLASMID_RADIUS + PLASMID_WIDTH / 2, PLASMID_RADIUS - PLASMID_WIDTH / 2, ang0, ang1)
        svg_parts.append(f'<path d="{path}" fill="{colour}" opacity="0.85" stroke="#333" stroke-width="0.6"/>')

        # arrow if > 5 % length and strand info
        length_bp = (end - start) % seq_len
        if feat.strand and length_bp > seq_len * 0.05:
            arrow_ang = ang0 + (ang1 - ang0) * (0.7 if feat.strand == 1 else 0.3)
            svg_parts.append(f'<polygon points="{_arrow_polygon(arrow_ang)}" fill="white" stroke="#333" stroke-width="0.6"/>')

        # label
        label = next((feat.qualifiers.get(k, [""])[0] for k in ("label", "gene", "product", "note") if k in feat.qualifiers), feat.type)
        mid = (ang0 + ang1) / 2 % (2 * math.pi)
        lx, ly = cx + (PLASMID_RADIUS + PLASMID_WIDTH / 2 + 24) * math.cos(mid), cy + (PLASMID_RADIUS + PLASMID_WIDTH / 2 + 24) * math.sin(mid)
        svg_parts.append(
            f'<text x="{lx:.1f}" y="{ly:.1f}" font-size="12" text-anchor="middle" dominant-baseline="middle"'
            f' transform="rotate({_text_rotate(mid):.1f} {lx:.1f} {ly:.1f})">{label}</text>'
        )

    # highlight region -------------------------------------------------------
    if highlight_region and {"start", "end"}.issubset(highlight_region):
        hs, he = highlight_region["start"], highlight_region["end"]
        ang0, ang1 = 2 * math.pi * hs / seq_len, 2 * math.pi * he / seq_len
        path = _arc_path(PLASMID_RADIUS + PLASMID_WIDTH + 10, PLASMID_RADIUS + PLASMID_WIDTH + 6, ang0, ang1)
        svg_parts.append(f'<path d="{path}" fill="none" stroke="#FF5722" stroke-dasharray="8,4" stroke-width="6" opacity="0.8"/>')
        mid = (ang0 + ang1) / 2
        tx, ty = cx + (PLASMID_RADIUS + PLASMID_WIDTH + 28) * math.cos(mid), cy + (PLASMID_RADIUS + PLASMID_WIDTH + 28) * math.sin(mid)
        svg_parts.append(
            f'<text x="{tx:.1f}" y="{ty:.1f}" font-size="14" font-weight="bold" fill="#FF5722" text-anchor="middle"'
            f' transform="rotate({_text_rotate(mid):.1f} {tx:.1f} {ty:.1f})">INSERT</text>'
        )

    # title
    svg_parts.append(f'<text x="{cx}" y="{cy - PLASMID_RADIUS - 40}" text-anchor="middle" font-size="26" font-weight="bold">{record.id}</text>')
    svg_parts.append(f'<text x="{cx}" y="{cy - PLASMID_RADIUS - 18}" text-anchor="middle" font-size="16" fill="#555">{seq_len} bp</text>')

    svg_parts.append('</svg>')
    svg = "\n".join(svg_parts)
    return base64.b64encode(svg.encode()).decode()

# ---------------------------------------------------------------------------
#   LINEAR MAP (matplotlib)
# ---------------------------------------------------------------------------

def visualize_linear_map(record: SeqRecord, highlight_region: Dict[str, Any] | None = None) -> str:
    seq_len = len(record.seq)
    fig, ax = plt.subplots(figsize=(12, 3), dpi=DPI_PNG)
    ax.set_ylim(-2, 2)
    ax.set_xlim(0, seq_len)
    ax.axis("off")
    ax.hlines(0, 0, seq_len, colors="#444", linewidth=2)

    # lane allocation
    lanes: List[Tuple[int, int, float]] = []  # (end, lane_idx, y)
    lane_height = 0.6

    def _next_lane(s: int) -> float:
        for idx, (end, lane_idx, y) in enumerate(lanes):
            if s > end:
                lanes[idx] = (s, lane_idx, y)
                return y
        y = (len(lanes) + 1) * lane_height * (-1 if len(lanes) % 2 else 1)
        lanes.append((s, len(lanes), y))
        return y

    for feat in sorted(record.features, key=lambda f: int(f.location.start)):
        if feat.type == "source":
            continue
        s, e = int(feat.location.start), int(feat.location.end)
        if e == s:
            continue
        y = _next_lane(s)
        colour = FEATURE_COLOURS.get(feat.type, DEFAULT_COLOUR)
        ax.add_patch(matplotlib.patches.Rectangle((s, y - 0.15), e - s, 0.3, color=colour, alpha=0.85))
        if feat.strand:
            head_x = e if feat.strand == 1 else s
            tri = np.array([[head_x, y], [head_x - 8 * feat.strand, y + 0.15], [head_x - 8 * feat.strand, y - 0.15]])
            ax.add_patch(matplotlib.patches.Polygon(tri, color=colour, alpha=0.85))
        label = next((feat.qualifiers.get(k, [""])[0] for k in ("label", "gene", "product", "note") if k in feat.qualifiers), feat.type)
        ax.text((s + e) / 2, y + 0.22 * (1 if y > 0 else -1), label, ha="center", va="center", fontsize=8, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))

    # highlight region
    if highlight_region and {"start", "end"}.issubset(highlight_region):
        hs, he = highlight_region["start"], highlight_region["end"]
        ax.axvspan(hs, he, color="#FFEB3B", alpha=0.3)

    # ticks every 1 kb
    for bp in range(0, seq_len + 1, TICK_INTERVAL):
        ax.vlines(bp, -0.2, 0.2, colors="#666", linewidth=1)
        ax.text(bp, -0.4, str(bp), ha="center", va="top", fontsize=7)

    ax.set_title(f"Linear Map • {record.id} ({seq_len} bp)", fontsize=12, pad=12)
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=DPI_PNG, bbox_inches="tight", transparent=True)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()
