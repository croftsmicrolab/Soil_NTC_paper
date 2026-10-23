from pathlib import Path
import hashlib
import json

import numpy as np
from chimerax.core.commands import run

root = Path(__file__).resolve().parent
out = root
(out / 'layers').mkdir(parents=True, exist_ok=True)
records = json.loads((root / 'input_manifest.json').read_text())
annotations = json.loads((root / 'gnat_annotations.json').read_text())['models']
ref_record = next(record for record in records if record['label'] == 'SatA')


def open_record(record):
    source = Path(record['source'])
    local_file = root / 'models' / source.parent.name / source.name
    assert hashlib.sha256(local_file.read_bytes()).hexdigest() == record['sha256']
    return next(model for model in run(session, f'open "{local_file}"') if hasattr(model, 'atoms'))


run(session, 'windowsize 1000 1000')
run(session, 'set bgColor white')
run(session, 'camera ortho')
run(session, 'lighting soft')
run(session, 'graphics silhouettes true width 1')
run(session, 'graphics quality 2')

ref = open_record(ref_record)
run(session, 'hide atoms')
run(session, f'show {ref.atomspec}/A cartoons')
run(session, f'view {ref.atomspec}/A orient')
run(session, 'zoom 0.90')
camera_rotation = session.main_view.camera.position.matrix[:, :3].copy()
ref_position = ref.position.matrix.copy()

results = []
for record in records:
    label = record['label']
    ref.display = True
    alignment_command = None
    alignment = None
    if label == 'SatA':
        model = ref
    else:
        model = open_record(record)
        alignment_command = (
            f'matchmaker {model.atomspec}/A to {ref.atomspec}/A '
            'pairing ss cutoffDistance 2.0 showAlignment false'
        )
        match_result = run(session, alignment_command)
        assert match_result and len(match_result) == 1, f'MatchMaker failed for {label}'
        match_result = match_result[0]
        assert len(match_result['final match atoms']) >= 3
        assert np.allclose(ref.position.matrix, ref_position), 'Reference moved'
        alignment = {
            'pruned_rmsd_angstrom': float(match_result['final RMSD']),
            'pruned_pairs': int(len(match_result['final match atoms'])),
            'all_pairs_rmsd_angstrom': float(match_result['full RMSD']),
            'all_pairs': int(len(match_result['full match atoms'])),
        }
        ref.display = False

    model.display = True
    run(session, f'hide {model.atomspec} atoms')
    run(session, f'show {model.atomspec}/A cartoons')
    ca_residues = [residue for residue in model.residues if residue.chain_id == 'A' and residue.find_atom('CA') is not None]
    assert ca_residues

    if label == 'SatA':
        for residue in ca_residues:
            residue.ribbon_color = (174, 181, 189, 255)
        alignment = {
            'pruned_rmsd_angstrom': 0.0,
            'pruned_pairs': len(ca_residues),
            'all_pairs_rmsd_angstrom': 0.0,
            'all_pairs': len(ca_residues),
        }
    else:
        numbers = [residue.number for residue in ca_residues]
        assert numbers == list(range(1, len(numbers) + 1)), f'Unexpected numbering for {label}'
        for residue in ca_residues:
            plddt = float(residue.find_atom('CA').bfactor)
            residue.ribbon_color = (
                (0, 83, 214, 255) if plddt >= 90 else
                (101, 203, 243, 255) if plddt >= 70 else
                (255, 219, 19, 255) if plddt >= 50 else
                (255, 125, 69, 255)
            )

    run(session, 'select clear')
    run(session, f'view {model.atomspec}/A')
    run(session, 'zoom 0.82')
    assert np.allclose(session.main_view.camera.position.matrix[:, :3], camera_rotation), 'Camera rotated'
    assert model == ref or not ref.display
    base_path = out / 'layers' / f'{label}_A_base.png'
    run(session, f'save "{base_path}" width 1400 height 1400 supersample 3')

    residue_numbers = {residue.number for residue in ca_residues}
    for index, (start, end) in enumerate(annotations[label]['regions']):
        assert any(start <= number <= end for number in residue_numbers)
        run(session, f'select {model.atomspec}/A:{start}-{end}')
        color = '#9467bd' if index == 0 else '#009e73'
        run(session, f'graphics selection color {color} width 7')
        region_path = out / 'layers' / f'{label}_A_region{index + 1}.png'
        run(session, f'save "{region_path}" width 1400 height 1400 supersample 3')

    run(session, 'select clear')
    updated = dict(record)
    updated.update(
        display_transform=model.position.matrix.tolist(),
        camera_transform=session.main_view.camera.position.matrix.tolist(),
        view_width_angstrom=float(session.main_view.camera.field_width),
        alignment_command=alignment_command,
        alignment_completed=bool(alignment_command),
        reference_hidden_in_panel=label != 'SatA',
        model_ca=len(ca_residues),
        alignment=alignment,
    )
    results.append(updated)
    (out / 'manifest.json').write_text(json.dumps(results, indent=2) + '\n')
    print('EXPORTED', label, alignment, flush=True)
    if model != ref:
        run(session, f'close {model.atomspec}')

run(session, 'exit')
