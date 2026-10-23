from pathlib import Path
import json

import numpy as np
from PIL import Image, ImageDraw, ImageFont

root = Path(__file__).resolve().parent
output = root / 'figures'
output.mkdir(exist_ok=True)
(root / 'panels').mkdir(exist_ok=True)
records = json.loads((root / 'manifest.json').read_text())
assert len(records) == 20
gnat = json.loads((root / 'gnat_annotations.json').read_text())['models']
phenotypes = json.loads((root / 'phenotype_labels.json').read_text())
red_labels = set(phenotypes['detected'])
assert len(red_labels) == 10 and red_labels <= {record['label'] for record in records}

# Combine independently rendered selection outlines at identical camera settings.
# Single-region panels are direct ChimeraX renders. For two regions, only the
# pixels changed relative to the unselected baseline are layered in source order.
for record in records:
    label = record['label']
    baseline = np.array(Image.open(root / 'layers' / f'{label}_A_base.png').convert('RGB'))
    result = baseline.copy()
    for index in range(len(gnat[label]['regions'])):
        layer = np.array(Image.open(root / 'layers' / f'{label}_A_region{index + 1}.png').convert('RGB'))
        mask = np.any(layer != baseline, axis=2)
        result[mask] = layer[mask]
    Image.fromarray(result).save(root / 'panels' / f'{label}_A.png')

font_path = '/System/Library/Fonts/Supplemental/Arial.ttf'
bold_path = '/System/Library/Fonts/Supplemental/Arial Bold.ttf'


def font(size, bold=False):
    return ImageFont.truetype(bold_path if bold else font_path, size)


width, height = 6000, 5420
cell_width, cell_height = 1160, 1170
margin = 100
top = 440
bounds = {}
for record in records:
    image = np.array(Image.open(root / 'panels' / f'{record["label"]}_A.png').convert('RGB'))
    ys, xs = np.where(np.min(image, axis=2) < 230)
    assert len(xs) > 1000
    assert min(xs) > 5 and max(xs) < 1395 and min(ys) > 5 and max(ys) < 1395, record['label']
    box = (min(xs), min(ys), max(xs) + 1, max(ys) + 1)
    size = int(max(box[2] - box[0], box[3] - box[1]) * 1.13)
    bounds[record['label']] = (box, size)

canvas = Image.new('RGB', (width, height), 'white')
draw = ImageDraw.Draw(canvas)
for index, record in enumerate(records):
    column, row = index % 5, index // 5
    x = margin + column * cell_width
    y = top + row * cell_height
    draw.text((x + 18, y), chr(97 + index), font=font(58, True), fill='#142333')
    draw.text(
        (x + 100, y + 3), record['label'], font=font(53, True),
        fill='#D73027' if record['label'] in red_labels else '#142333',
    )
    regions = gnat[record['label']]['regions']
    region_text = '; '.join(f'{start}–{end}' for start, end in regions)
    if record['label'] == 'SatA':
        status = f'3PP9 chain A | GNAT (aa): {region_text}'
    elif index == 1:
        status = f'REF | GNAT (aa): {region_text}'
    else:
        status = f'GNAT (aa): {region_text}'
    draw.text((x + 100, y + 78), status, font=font(29), fill='#657080')

    box, size = bounds[record['label']]
    center_x = (box[0] + box[2]) / 2
    center_y = (box[1] + box[3]) / 2
    crop = (
        int(center_x - size / 2), int(center_y - size / 2),
        int(center_x - size / 2) + size, int(center_y - size / 2) + size,
    )
    raw = Image.open(root / 'panels' / f'{record["label"]}_A.png').convert('RGB')
    padded = Image.new('RGB', (size, size), 'white')
    left, upper = max(0, crop[0]), max(0, crop[1])
    right, lower = min(1400, crop[2]), min(1400, crop[3])
    padded.paste(raw.crop((left, upper, right, lower)), (left - crop[0], upper - crop[1]))
    panel = padded.resize((1000, 1000), Image.Resampling.LANCZOS)
    panel_path = root / 'panels' / f'{record["label"]}_A_slide.png'
    panel.save(panel_path)
    canvas.paste(panel, (x + 60, y + 125))

y = 5160
draw.line((100, y - 45, 5900, y - 45), fill='#ced6dd', width=3)
draw.text((100, y), 'pLDDT fill (AF3 models)', font=font(38, True), fill='#243442')
for index, (color, label) in enumerate([
    ('#0053d6', '90–100'), ('#65cbf3', '70–<90'), ('#ffdb13', '50–<70'), ('#ff7d45', '<50'),
]):
    x = 620 + index * 540
    draw.rectangle((x, y + 4, x + 55, y + 46), fill=color)
    draw.text((x + 72, y), label, font=font(35), fill='#243442')
draw.rectangle((2870, y + 4, 2925, y + 46), fill='#AEB5BD')
draw.text((2945, y), 'SatA: 3PP9 experimental', font=font(35), fill='#243442')
draw.text((4240, y), 'Red names: NTC resistance (n = 10)', font=font(35), fill='#D73027')

y = 5250
draw.rectangle((100, y + 6, 185, y + 48), outline='#9467bd', width=12)
draw.text((205, y), 'GNAT region 1 (outline)', font=font(38), fill='#243442')
draw.rectangle((1420, y + 6, 1505, y + 48), outline='#009e73', width=12)
draw.text((1525, y), 'GNAT region 2: HIT_000035 (outline)', font=font(38), fill='#243442')

canvas.save(output / 'Figure_aligned.png', dpi=(600, 600))
canvas.save(output / 'Figure_aligned.tiff', compression='tiff_lzw', dpi=(600, 600))
canvas.resize((1500, 1355)).save(root / 'preview_A.png')
