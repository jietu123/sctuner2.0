from __future__ import annotations

import math
import random
from html import escape
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class SVG:
    def __init__(self, width: int, height: int, title: str):
        self.width = width
        self.height = height
        self.parts = [
            '<?xml version="1.0" encoding="UTF-8"?>',
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
            f'width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            f'<title>{escape(title)}</title>',
            '<desc>Fully editable vector reconstruction. Text, points, heatmaps, arrows, and panels are separate SVG elements.</desc>',
        ]

    @staticmethod
    def attrs(**attrs) -> str:
        converted = []
        for key, value in attrs.items():
            if value is None:
                continue
            key = key.rstrip('_').replace('_', '-')
            converted.append(f'{key}="{escape(str(value), quote=True)}"')
        return ' '.join(converted)

    def raw(self, value: str):
        self.parts.append(value)

    def start(self, tag: str, **attrs):
        self.parts.append(f'<{tag} {self.attrs(**attrs)}>')

    def end(self, tag: str):
        self.parts.append(f'</{tag}>')

    def el(self, tag: str, **attrs):
        self.parts.append(f'<{tag} {self.attrs(**attrs)}/>')

    def text(self, x, y, value, size=16, weight=400, anchor='start', fill='#111827', **attrs):
        defaults = {
            'x': x,
            'y': y,
            'font_family': 'Arial, Helvetica, sans-serif',
            'font_size': size,
            'font_weight': weight,
            'text_anchor': anchor,
            'fill': fill,
        }
        defaults.update(attrs)
        self.parts.append(f'<text {self.attrs(**defaults)}>{escape(str(value))}</text>')

    def multiline(self, x, y, lines, size=16, weight=400, anchor='middle', fill='#111827', line_height=1.25, **attrs):
        defaults = {
            'x': x,
            'y': y,
            'font_family': 'Arial, Helvetica, sans-serif',
            'font_size': size,
            'font_weight': weight,
            'text_anchor': anchor,
            'fill': fill,
        }
        defaults.update(attrs)
        tspans = []
        for i, line in enumerate(lines):
            dy = 0 if i == 0 else size * line_height
            tspans.append(f'<tspan x="{x}" dy="{dy}">{escape(str(line))}</tspan>')
        self.parts.append(f'<text {self.attrs(**defaults)}>{"".join(tspans)}</text>')

    def finish(self, path: Path):
        self.parts.append('</svg>')
        path.write_text('\n'.join(self.parts) + '\n', encoding='utf-8')


NAVY = '#102a72'
BLUE = '#42658b'
PANEL_BLUE = '#eef5ff'
GRID = '#d9dee8'
TEXT = '#111827'
PURPLE = '#6f42b5'
PURPLE_LIGHT = '#b99add'
RED = '#ef5350'
ORANGE = '#f59e0b'
GREEN = '#62b53f'
TEAL = '#25a9b0'
PINK = '#ee5b91'
BROWN = '#b36a3c'
GREY = '#bfc3c7'


def layer(svg: SVG, layer_id: str, label: str):
    svg.start('g', id=layer_id, inkscape_groupmode='layer', inkscape_label=label)


def panel(svg: SVG, x, y, w, h, title, title_size=21, radius=2, border='#252b33', header_h=54):
    svg.el('rect', x=x, y=y, width=w, height=h, rx=radius, fill='#ffffff', stroke=border, stroke_width=1.5)
    svg.el('rect', x=x, y=y, width=w, height=header_h, rx=radius, fill=PANEL_BLUE)
    svg.el('line', x1=x, y1=y + header_h, x2=x + w, y2=y + header_h, stroke=border, stroke_width=1.5)
    svg.text(x + w / 2, y + 35, title, size=title_size, weight=700, anchor='middle', fill='#111111')


def arrow(svg: SVG, x1, y, x2, color='#111111', width=3, head=18):
    shaft_end = x2 - head
    svg.el('line', x1=x1, y1=y, x2=shaft_end, y2=y, stroke=color, stroke_width=width)
    svg.el('path', d=f'M {shaft_end} {y-head*0.72} L {x2} {y} L {shaft_end} {y+head*0.72} Z', fill='white', stroke=color, stroke_width=width)


def blue_arrow(svg: SVG, x1, y, x2):
    svg.el('line', x1=x1, y1=y, x2=x2 - 15, y2=y, stroke=BLUE, stroke_width=6)
    svg.el('path', d=f'M {x2-15} {y-12} L {x2} {y} L {x2-15} {y+12} Z', fill=BLUE)


def heat_color(value: float) -> str:
    stops = [
        (0.0, (33, 19, 78)),
        (0.25, (93, 29, 129)),
        (0.5, (161, 42, 143)),
        (0.75, (224, 95, 93)),
        (1.0, (255, 210, 80)),
    ]
    value = max(0.0, min(1.0, value))
    for (a, ca), (b, cb) in zip(stops, stops[1:]):
        if a <= value <= b:
            t = (value - a) / (b - a)
            rgb = tuple(round(ca[i] + t * (cb[i] - ca[i])) for i in range(3))
            return '#%02x%02x%02x' % rgb
    return '#ffd250'


def diverging_color(value: float) -> str:
    value = max(-1.0, min(1.0, value))
    if value < 0:
        t = value + 1
        a, b = (69, 62, 175), (245, 245, 245)
    else:
        t = value
        a, b = (245, 245, 245), (222, 45, 38)
    rgb = tuple(round(a[i] + t * (b[i] - a[i])) for i in range(3))
    return '#%02x%02x%02x' % rgb


def draw_heatmap(svg: SVG, x, y, rows, cols, cw, ch, values, palette='heat', stroke='#ffffff', stroke_width=.45):
    for r in range(rows):
        for c in range(cols):
            value = values[r][c]
            fill = heat_color(value) if palette == 'heat' else diverging_color(value)
            svg.el('rect', x=x + c * cw, y=y + r * ch, width=cw, height=ch,
                   fill=fill, stroke=stroke, stroke_width=stroke_width)
    svg.el('rect', x=x, y=y, width=cols*cw, height=rows*ch, fill='none', stroke='#4b5563', stroke_width=1)


def gaussian_clusters(svg: SVG, centers, colors, points_per=48, sx=18, sy=18, radius=4.5, seed=1, opacity=1.0):
    rng = random.Random(seed)
    for idx, ((cx, cy), color) in enumerate(zip(centers, colors)):
        svg.start('g', id=f'cluster-{seed}-{idx}', fill=color, opacity=opacity)
        for _ in range(points_per):
            angle = rng.random() * math.tau
            distance = min(2.0, abs(rng.gauss(.75, .38)))
            px = cx + math.cos(angle) * distance * sx
            py = cy + math.sin(angle) * distance * sy
            svg.el('circle', cx=f'{px:.2f}', cy=f'{py:.2f}', r=radius)
        svg.end('g')


def hex_circle(svg: SVG, cx, cy, radius, spacing=10, dot=4.2, mode='purple', seed=1):
    rng = random.Random(seed)
    dy = spacing * .8660254
    rows = int(2 * radius / dy) + 2
    points = []
    for r in range(-rows, rows + 1):
        py = cy + r * dy
        offset = spacing / 2 if r % 2 else 0
        for c in range(-rows, rows + 1):
            px = cx + c * spacing + offset
            if (px-cx)**2 + (py-cy)**2 <= radius**2:
                points.append((px, py))
    svg.start('g', id=f'hex-map-{seed}')
    palette = [RED, ORANGE, GREEN, TEAL, '#4f86d9', '#8a5bc3']
    for px, py in points:
        if mode == 'purple':
            d = math.hypot(px-cx, py-cy) / radius
            val = max(0, min(1, .15 + .58 * rng.random() + .25*(1-d)))
            fill = heat_color(.15 + .5 * val)
        elif mode == 'spatial':
            angle = (math.atan2(py-cy, px-cx) + math.tau) % math.tau
            sector = int((angle + .32) / math.tau * len(palette)) % len(palette)
            fill = palette[sector]
        else:
            fill = PURPLE
        svg.el('circle', cx=f'{px:.2f}', cy=f'{py:.2f}', r=dot, fill=fill, stroke='#ffffff', stroke_width=.45)
    svg.end('g')


def tissue_points(svg: SVG, x, y, w, h, seed=1, mode='purple', spacing=9.2, dot=3.5):
    rng = random.Random(seed)
    colors = [PINK, GREEN, ORANGE, '#3d78bd', PURPLE, BROWN, GREY]
    pts = []
    dy = spacing * .866
    row = 0
    py = y
    while py <= y+h:
        px = x + (spacing/2 if row % 2 else 0)
        while px <= x+w:
            nx = (px - (x+w*.51))/(w*.5)
            ny = (py - (y+h*.50))/(h*.46)
            boundary = nx*nx + ny*ny
            notch = (px > x+w*.65 and py > y+h*.62 and py < y+h*.84)
            left_cut = px < x+w*.08 and py < y+h*.26
            top_bulge = math.exp(-((nx+.35)/.36)**2) * .18
            if boundary < 1 + top_bulge and not notch and not left_cut:
                pts.append((px, py, nx, ny))
            px += spacing
        py += dy
        row += 1
    svg.start('g', id=f'tissue-map-{seed}')
    for px, py, nx, ny in pts:
        if mode == 'purple':
            fill = heat_color(.18 + .45*rng.random() + .15*(1-abs(ny)))
        elif mode == 'corr':
            field = math.sin(px*.065) + math.cos(py*.055) + rng.gauss(0, .72)
            fill = diverging_color(max(-1, min(1, field/2.25)))
        else:
            score = [
                nx + .9*ny,
                -nx + .4*ny,
                ny - .3*nx,
                -ny + .1*nx,
                math.sin(nx*3)+math.cos(ny*3),
                nx*nx+ny,
            ]
            fill = colors[max(range(len(score)), key=lambda i: score[i])]
        svg.el('circle', cx=f'{px:.2f}', cy=f'{py:.2f}', r=dot, fill=fill, stroke='#ffffff', stroke_width=.35)
    svg.end('g')


def number_badge(svg: SVG, x, y, number):
    svg.el('circle', cx=x, cy=y, r=16, fill=NAVY)
    svg.text(x, y+7, number, size=20, weight=700, anchor='middle', fill='#ffffff')


def vertical_gradient_legend(svg: SVG, x, y, h, top, bottom, high='High', low='Low'):
    gid = f'gradient-{int(x)}-{int(y)}'
    svg.raw(f'<defs><linearGradient id="{gid}" x1="0" y1="1" x2="0" y2="0"><stop offset="0" stop-color="{bottom}"/><stop offset="1" stop-color="{top}"/></linearGradient></defs>')
    svg.el('rect', x=x, y=y, width=16, height=h, fill=f'url(#{gid})', stroke='#4b5563', stroke_width=.8)
    svg.text(x+8, y-10, high, size=13, anchor='middle')
    svg.text(x+8, y+h+20, low, size=13, anchor='middle')


def build_figure_1():
    svg = SVG(2172, 724, 'Construction of in silico simulated ST data')
    svg.el('rect', x=0, y=0, width=2172, height=724, fill='#ffffff')
    layer(svg, 'title-layer', 'Title')
    svg.text(25, 76, 'b', size=68, weight=700, fill='#000000')
    svg.raw('<text x="1086" y="82" font-family="Arial, Helvetica, sans-serif" font-size="46" font-weight="700" text-anchor="middle" fill="#000000">Construction of <tspan font-style="italic">in silico</tspan> simulated ST data</text>')
    svg.end('g')

    y, h = 142, 432
    layer(svg, 'inputs', 'Inputs')
    panel(svg, 26, y, 416, h, 'Inputs', title_size=22)
    svg.el('line', x1=237, y1=y+54, x2=237, y2=y+h, stroke='#60656d', stroke_width=1.2)
    svg.text(132, 241, '(a) scRNA-seq reference', size=16, weight=700, anchor='middle')
    svg.text(341, 241, '(b) Real ST scaffold', size=16, weight=700, anchor='middle')
    gaussian_clusters(svg,
                      [(84, 332), (145, 302), (191, 349), (65, 408), (126, 411), (189, 415), (103, 477), (171, 489), (206, 390)],
                      [RED, RED, '#91b51d', ORANGE, TEAL, '#3c83c9', GREEN, '#d64ac0', '#f58a8f'],
                      points_per=25, sx=18, sy=18, radius=4.6, seed=12)
    svg.el('line', x1=52, y1=532, x2=109, y2=532, stroke='#000', stroke_width=2)
    svg.el('path', d='M 109 532 L 100 527 L 100 537 Z', fill='#000')
    svg.el('line', x1=52, y1=532, x2=52, y2=469, stroke='#000', stroke_width=2)
    svg.el('path', d='M 52 469 L 47 478 L 57 478 Z', fill='#000')
    svg.text(78, 553, 'UMAP1', size=13, anchor='middle')
    svg.text(42, 501, 'UMAP2', size=13, anchor='middle', transform='rotate(-90 42 501)')
    svg.el('circle', cx=341, cy=395, r=112, fill='#f8f8fc', stroke='#c6cad4', stroke_width=1.2)
    hex_circle(svg, 341, 395, 101, spacing=11.5, dot=4.8, mode='purple', seed=31)
    svg.end('g')

    layer(svg, 'build-profiles', 'Build type profiles and spatial truth')
    panel(svg, 497, y, 448, h, 'Build type profiles & spatial truth', title_size=21)
    svg.text(517, 233, '•  derive type profiles', size=17)
    svg.text(517, 261, '•  build clustered spatial truth', size=17)
    svg.text(605, 313, 'Type profiles', size=17, weight=600, anchor='middle')
    svg.text(819, 313, 'Spatial truth', size=17, weight=600, anchor='middle')
    type_colors = [RED, ORANGE, GREEN, TEAL, '#5da4d5', PURPLE]
    for i, color in enumerate(type_colors):
        svg.el('rect', x=530+i*27, y=330, width=25, height=10, fill=color)
    vals = [[max(0, min(1, .15 + .56*random.Random(80+r*17+c).random() + (.25 if c == r % 6 else 0))) for c in range(6)] for r in range(11)]
    draw_heatmap(svg, 530, 346, 11, 6, 27, 17, vals)
    svg.text(611, 556, 'Cell types', size=15, anchor='middle')
    svg.text(516, 447, 'Genes', size=15, anchor='middle', transform='rotate(-90 516 447)')
    hex_circle(svg, 819, 433, 108, spacing=11, dot=4.7, mode='spatial', seed=52)
    svg.end('g')

    layer(svg, 'synthesis', 'Synthesize simulated ST expression')
    panel(svg, 995, y, 633, h, 'Synthesize simulated ST expression', title_size=21)
    svg.text(1015, 233, '•  mix type profiles', size=17)
    svg.text(1015, 261, '•  sample spot fractions', size=17)
    svg.text(1083, 334, 'Type profiles', size=16, weight=600, anchor='middle')
    for i, color in enumerate(type_colors):
        svg.el('rect', x=1031+i*20, y=350, width=18, height=9, fill=color)
    vals2 = [[max(0, min(1, .12 + .63*random.Random(140+r*11+c).random() + (.2 if c == (r+2) % 6 else 0))) for c in range(6)] for r in range(8)]
    draw_heatmap(svg, 1031, 364, 8, 6, 20, 20, vals2)
    svg.text(1091, 546, 'Cell types', size=15, anchor='middle')
    svg.text(1018, 444, 'Genes', size=15, anchor='middle', transform='rotate(-90 1018 444)')
    svg.text(1188, 442, '×', size=40, weight=400, anchor='middle')
    svg.text(1260, 334, 'Spot fractions', size=16, weight=600, anchor='middle')
    frac = [[random.Random(220+r*9+c).random() for c in range(6)] for r in range(6)]
    for r in range(6):
        s = sum(frac[r]); frac[r] = [v/s*2.2 for v in frac[r]]
    for r in range(6):
        for c in range(6):
            v = min(1, frac[r][c])
            fill = f'#{int(236-160*v):02x}{int(248-90*v):02x}{int(239-165*v):02x}'
            svg.el('rect', x=1202+c*20, y=352+r*25, width=20, height=25, fill=fill, stroke='#cfd5d8', stroke_width=.7)
    svg.el('rect', x=1202, y=352, width=120, height=150, fill='none', stroke='#4b5563', stroke_width=1)
    svg.text(1262, 526, 'Spots', size=15, anchor='middle')
    svg.text(1190, 429, 'Cell types', size=15, anchor='middle', transform='rotate(-90 1190 429)')
    svg.el('rect', x=1227, y=540, width=70, height=12, fill='#e8f4e9', stroke='#d0d6d0', stroke_width=.5)
    for i in range(7):
        svg.el('rect', x=1227+i*10, y=540, width=10, height=12, fill=f'#{int(235-160*i/6):02x}{int(246-80*i/6):02x}{int(236-165*i/6):02x}')
    svg.text(1217, 552, '0', size=13, anchor='end'); svg.text(1307, 552, '1', size=13)
    svg.text(1353, 442, '×', size=40, anchor='middle')
    svg.text(1407, 334, 'Library size', size=16, weight=600, anchor='middle')
    for i in range(6):
        fill = f'#{220-i*18:02x}{239-i*15:02x}{248-i*6:02x}'
        svg.el('rect', x=1397, y=350+i*26, width=20, height=26, fill=fill, stroke='#5d6875', stroke_width=.7)
    svg.text(1407, 534, 'Spots', size=15, anchor='middle')
    svg.text(1461, 442, '=', size=37, weight=700, anchor='middle')
    svg.multiline(1547, 310, ['Simulated ST', 'matrix'], size=16, weight=600)
    simvals = [[random.Random(370+r*8+c).random() for c in range(6)] for r in range(6)]
    for r in range(6):
        for c in range(6):
            v = simvals[r][c]
            fill = f'#{int(235-145*v):02x}{int(246-95*v):02x}{int(252-42*v):02x}'
            svg.el('rect', x=1491+c*19, y=349+r*26, width=19, height=26, fill=fill, stroke='#7e91a2', stroke_width=.6)
    svg.el('rect', x=1491, y=349, width=114, height=156, fill='none', stroke='#4b5563', stroke_width=1)
    svg.text(1548, 531, 'Spots', size=15, anchor='middle')
    svg.end('g')

    layer(svg, 'outputs', 'Outputs')
    panel(svg, 1679, y, 340, h, 'Outputs', title_size=22)
    svg.text(1758, 309, 'Simulated ST map', size=17, weight=600, anchor='middle')
    svg.el('circle', cx=1758, cy=425, r=84, fill='#f8f8fc', stroke='#c6cad4', stroke_width=1.2)
    hex_circle(svg, 1758, 425, 77, spacing=9.5, dot=3.8, mode='purple', seed=88)
    svg.text(1850, 375, '•  simulated ST matrix', size=15)
    svg.text(1850, 411, '•  cell-to-spot truth', size=15)
    svg.text(1850, 447, '•  spot-type fractions', size=15)
    svg.end('g')

    layer(svg, 'flow-arrows', 'Flow arrows')
    arrow(svg, 448, 360, 489)
    arrow(svg, 951, 360, 987)
    arrow(svg, 1635, 360, 1671)
    svg.end('g')
    svg.finish(ROOT / '1_editable.svg')


def step_panel(svg: SVG, x, y, w, h, number, title_lines):
    svg.el('rect', x=x, y=y, width=w, height=h, rx=18, fill='#ffffff', stroke='#4169b1', stroke_width=1.2)
    number_badge(svg, x+38, y+38, number)
    svg.multiline(x+w/2+16, y+43, title_lines, size=19, weight=700, anchor='middle', fill=NAVY, line_height=1.15)


def build_figure_1b():
    svg = SVG(1706, 922, 'SVTuner workflow')
    svg.el('rect', x=0, y=0, width=1706, height=922, fill='#ffffff')
    layer(svg, 'title-layer', 'Title')
    svg.text(853, 53, 'SVTuner workflow', size=42, weight=700, anchor='middle', fill='#000000')
    svg.end('g')
    y, h = 83, 800
    xs = [15, 339, 635, 1067, 1397]
    ws = [288, 269, 404, 295, 294]

    layer(svg, 'step-1-input-data', '1 Input data')
    step_panel(svg, xs[0], y, ws[0], h, 1, ['Input data'])
    svg.multiline(159, 189, ['ST expression', '(spots × genes)'], size=16, weight=600)
    vertical_gradient_legend(svg, 32, 258, 104, '#56349a', '#f5f2f9')
    tissue_points(svg, 68, 233, 222, 170, seed=101, mode='purple', spacing=9.4, dot=3.8)
    svg.text(159, 460, '+', size=52, weight=700, anchor='middle', fill='#000000')
    svg.multiline(159, 510, ['scRNA-seq reference', '(cells × genes)'], size=16, weight=600)
    gaussian_clusters(svg,
                      [(70, 613), (157, 584), (237, 632), (58, 715), (112, 672), (190, 736), (244, 744), (188, 685)],
                      [PINK, GREEN, ORANGE, GREY, '#3d78bd', ORANGE, PURPLE, GREY],
                      points_per=28, sx=18, sy=20, radius=4.5, seed=202)
    svg.end('g')

    layer(svg, 'step-2-marker-construction', '2 Marker construction')
    step_panel(svg, xs[1], y, ws[1], h, 2, ['Marker construction', 'for each cell type'])
    legend_y = 211
    marker_colors = [PINK, GREEN, ORANGE, '#ffffff']
    marker_text = ['Type 1 markers', 'Type 2 markers', 'Type 3 markers', 'Type K markers']
    for i, (color, label) in enumerate(zip(marker_colors, marker_text)):
        yy = legend_y + i*46
        svg.el('circle', cx=392, cy=yy, r=8, fill=color, stroke='#343a40', stroke_width=1)
        svg.text(414, yy+6, label, size=16)
    svg.text(422, 339, '...', size=18, anchor='middle')
    for i, color in enumerate([PINK, GREEN, BROWN, '#4f7ebd', PURPLE]):
        svg.el('rect', x=385+i*29, y=454, width=27, height=18, fill=color, opacity=.75, stroke='#ffffff', stroke_width=.5)
    values = []
    rng = random.Random(454)
    for r in range(7):
        row = []
        for c in range(7):
            diag = math.exp(-((r-c)/1.2)**2)
            row.append(max(-1, min(1, -0.55 + 1.5*diag + rng.uniform(-.32,.32))))
        values.append(row)
    draw_heatmap(svg, 369, 484, 7, 7, 25, 27, values, palette='diverge', stroke='#707782', stroke_width=.55)
    svg.text(356, 581, 'Genes', size=14, weight=600, anchor='middle', transform='rotate(-90 356 581)')
    vertical_gradient_legend(svg, 564, 526, 113, '#d92522', '#4b46ad')
    svg.multiline(476, 719, ['Marker expression', 'profile (avg. in scRNA-seq)'], size=16, weight=600)
    svg.end('g')

    layer(svg, 'step-3-support-score', '3 Compute ST support score')
    step_panel(svg, xs[2], y, ws[2], h, 3, ['Compute ST support score', 'for each cell type'])
    svg.multiline(735, 190, ['Marker profile', '(type j)'], size=16, weight=600)
    for i in range(6):
        svg.el('rect', x=661+i*25, y=229, width=25, height=25,
               fill=diverging_color(.95-i*.3), stroke='#7c4a42', stroke_width=.8)
    svg.text(736, 291, '×', size=34, anchor='middle')
    svg.text(687, 315, 'ST spots', size=15, weight=600)
    # Three miniature editable ST spot maps.
    for i, cx in enumerate([679, 746, 813]):
        tissue_points(svg, cx-27, 335, 55, 88, seed=510+i, mode='corr' if i < 2 else 'purple', spacing=8, dot=2.2)
        if i < 2:
            svg.text(cx+34, 382, '·', size=25, anchor='middle')
    svg.text(842, 382, '...', size=18)
    svg.multiline(939, 190, ['Spot-wise correlation', '(Fisher r→z)'], size=16, weight=600)
    tissue_points(svg, 862, 229, 110, 220, seed=530, mode='corr', spacing=8.3, dot=3.1)
    vertical_gradient_legend(svg, 1002, 273, 111, '#e02b27', '#4e4ab0')
    svg.multiline(837, 481, ['Aggregate to support score', '(max or top-k mean)'], size=16, weight=600)
    svg.text(836, 540, '▼', size=20, anchor='middle')
    gid = 'support-score-gradient'
    svg.raw(f'<defs><linearGradient id="{gid}" x1="0" y1="0" x2="1" y2="0"><stop offset="0" stop-color="#d7e9fb"/><stop offset="0.5" stop-color="#ffffff"/><stop offset="1" stop-color="#e4e4e4"/></linearGradient></defs>')
    svg.el('rect', x=699, y=543, width=253, height=15, fill=f'url(#{gid})', stroke='#545b66', stroke_width=1)
    svg.text(699, 534, '0', size=15, anchor='middle'); svg.text(952, 534, '1', size=15, anchor='middle')
    svg.el('line', x1=836, y1=538, x2=836, y2=575, stroke='#111', stroke_width=3)
    svg.text(836, 592, '0.28', size=16, anchor='middle')
    svg.text(835, 628, 'Support score for all cell types', size=16, weight=600, anchor='middle')
    types = [('Type 1', .16, PINK), ('Type 2', .32, GREEN), ('Type 3', .55, ORANGE), ('...', .57, PURPLE), ('Type K', .82, GREY)]
    chart_x, chart_y, chart_w = 716, 657, 268
    svg.el('line', x1=chart_x, y1=chart_y-16, x2=chart_x, y2=810, stroke='#20252b', stroke_width=1.2)
    svg.el('line', x1=chart_x, y1=810, x2=chart_x+chart_w, y2=810, stroke='#20252b', stroke_width=1.2)
    for i, (name, score, color) in enumerate(types):
        yy = chart_y+i*34
        svg.text(chart_x-16, yy+5, name, size=15, anchor='end')
        svg.el('line', x1=chart_x, y1=yy, x2=chart_x+chart_w, y2=yy, stroke='#aeb4bd', stroke_width=1, stroke_dasharray='5 5')
        svg.el('circle', cx=chart_x+score*chart_w, cy=yy, r=7, fill=color, stroke='#ffffff', stroke_width=1)
    for score, label in [(0,'0'),(.5,'0.5'),(1,'1.0')]:
        xx = chart_x+score*chart_w
        svg.el('line', x1=xx, y1=810, x2=xx, y2=816, stroke='#20252b', stroke_width=1)
        svg.text(xx, 836, label, size=14, anchor='middle')
    svg.text(chart_x+chart_w/2, 862, 'Support score', size=16, weight=600, anchor='middle')
    svg.end('g')

    layer(svg, 'step-4-classification', '4 Classify cell type status')
    step_panel(svg, xs[3], y, ws[3], h, 4, ['Classify cell type status'])
    statuses = [
        (209, '#f4f9ed', '#75b84e', ['Supported', '(strong evidence)']),
        (326, '#fff9ee', ORANGE, ['Weak / Grey zone', '(uncertain)']),
        (445, '#fff4f3', '#ff625b', ['Missing-like', '(weak marker identity', '& ST evidence)']),
        (575, '#f3effa', '#8b6dbd', ['Unsupported', '(no evidence)']),
    ]
    for yy, fill, stroke, lines in statuses:
        box_h = 84 if len(lines) == 2 else 92
        svg.el('rect', x=1093, y=yy, width=240, height=box_h, rx=9, fill=fill, stroke=stroke, stroke_width=1.2)
        svg.multiline(1213, yy+36, lines, size=18, weight=600, line_height=1.28)
    svg.el('rect', x=1082, y=704, width=265, height=130, rx=10, fill='#f7f7f6', stroke='#d7d7d4', stroke_width=1)
    svg.text(1214, 735, 'Conservative rescue & protection', size=15, weight=600, anchor='middle')
    # Shield, magnifier, swap, and group icons are intentionally simple editable paths.
    svg.el('path', d='M1100 772 L1117 765 L1134 772 L1132 792 Q1128 807 1117 812 Q1106 807 1102 792 Z', fill='none', stroke='#111', stroke_width=2)
    svg.el('circle', cx=1171, cy=786, r=13, fill='none', stroke='#111', stroke_width=2)
    svg.el('circle', cx=1171, cy=786, r=9, fill='none', stroke='#111', stroke_width=1.5)
    svg.el('line', x1=1181, y1=796, x2=1192, y2=807, stroke='#111', stroke_width=2)
    for cx in [1223, 1266]:
        svg.el('circle', cx=cx, cy=778, r=6, fill='none', stroke='#111', stroke_width=1.7)
        svg.el('path', d=f'M {cx-9} 803 Q {cx-9} 786 {cx} 786 Q {cx+9} 786 {cx+9} 803 Z', fill='none', stroke='#111', stroke_width=1.7)
    svg.text(1245, 797, '↔', size=24, anchor='middle')
    for cx, cy in [(1301,786),(1322,786),(1311,770)]:
        svg.el('circle', cx=cx, cy=cy, r=5, fill='none', stroke='#111', stroke_width=1.5)
        svg.el('circle', cx=cx, cy=cy+11, r=7, fill='none', stroke='#111', stroke_width=1.5)
    svg.end('g')

    layer(svg, 'step-5-filtered-pool', '5 Generate filtered mapping pool')
    step_panel(svg, xs[4], y, ws[4], h, 5, ['Generate plugin_type', 'and filtered mapping pool'])
    svg.text(1544, 205, 'Action for each cell type', size=17, weight=600, anchor='middle')
    actions = [
        (246, GREEN, ['Keep', '(enter mapping pool)']),
        (312, ORANGE, ['Relabel', '(map to similar type)']),
        (378, PINK, ['Mark Unknown', '(set as Unknown_sc_only)']),
        (444, PURPLE_LIGHT, ['Drop / Filter', '(remove from pool)']),
        (510, GREY, ['...']),
    ]
    for yy, color, lines in actions:
        svg.el('circle', cx=1446, cy=yy, r=8.5, fill=color, stroke='#666', stroke_width=.7)
        svg.multiline(1476, yy+5, lines, size=16, weight=600 if len(lines)==1 else 500, anchor='start', line_height=1.25)
    svg.el('line', x1=1414, y1=551, x2=1674, y2=551, stroke='#777', stroke_width=1.3, stroke_dasharray='6 6')
    svg.multiline(1544, 592, ['Filtered scRNA-seq pool', '(plugin_type)'], size=17, weight=600)
    gaussian_clusters(svg,
                      [(1477, 682), (1547, 731), (1477, 787), (1590, 795)],
                      [PINK, '#3d78bd', GREEN, PURPLE], points_per=30, sx=17, sy=20, radius=4.4, seed=808)
    svg.el('circle', cx=1624, cy=674, r=22, fill='none', stroke='#222', stroke_width=1.5, stroke_dasharray='7 5')
    svg.el('circle', cx=1659, cy=721, r=22, fill='none', stroke='#222', stroke_width=1.5, stroke_dasharray='7 5')
    svg.end('g')

    layer(svg, 'workflow-arrows', 'Workflow arrows')
    blue_arrow(svg, 307, 442, 336)
    blue_arrow(svg, 611, 442, 632)
    blue_arrow(svg, 1042, 442, 1064)
    blue_arrow(svg, 1366, 442, 1394)
    svg.end('g')
    svg.finish(ROOT / 'FIG1B_editable.svg')


if __name__ == '__main__':
    build_figure_1()
    build_figure_1b()
    print(ROOT / '1_editable.svg')
    print(ROOT / 'FIG1B_editable.svg')
