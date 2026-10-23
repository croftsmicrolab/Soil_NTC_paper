from pathlib import Path
import json

from chimerax.core.commands import run

root = Path(__file__).resolve().parent
records = json.loads((root / 'input_manifest.json').read_text())

record_by_label = {record['label']: record for record in records}
record_by_label['2QIR'] = {
    'label': '2QIR',
    'source': str(root / 'models/2QIR/2QIR_AAC_only.pdb'),
    'chain': 'A',
}

labels = ['SatA', '2QIR'] + [record['label'] for record in records if record['label'] != 'SatA']
references = [
    ('SatA', 'PDB 3PP9 chain A'),
    ('2QIR', 'PDB 2QIR chain A (AAC-only structure)'),
]


def open_atomic(record):
    source = Path(record['source'])
    path = root / 'models' / source.parent.name / source.name
    model = next(model for model in run(session, f'open "{path}"') if hasattr(model, 'atoms'))
    return model


def ca_count(model, chain_id='A'):
    return sum(
        1 for residue in model.residues
        if residue.chain_id == chain_id and residue.find_atom('CA') is not None
    )


all_results = {}
for reference_label, reference_description in references:
    reference_record = record_by_label[reference_label]
    reference_chain = reference_record.get('chain', 'A')
    reference_model = open_atomic(reference_record)
    reference_ca = ca_count(reference_model, reference_chain)
    table = []

    for label in labels:
        record = record_by_label[label]
        chain = record.get('chain', 'A')
        if label == reference_label:
            row = {
                'Protein': label,
                'Pruned_RMSD_A': 0.0,
                'Pruned_pairs': reference_ca,
                'All_pairs_RMSD_A': 0.0,
                'All_pairs': reference_ca,
                'Model_CA': reference_ca,
            }
        else:
            model = open_atomic(record)
            model_ca = ca_count(model, chain)
            command = (
                f'matchmaker {model.atomspec}/{chain} to '
                f'{reference_model.atomspec}/{reference_chain} '
                'pairing ss cutoffDistance 2.0 showAlignment false'
            )
            result = run(session, command)
            assert result and len(result) == 1, f'MatchMaker failed: {label} to {reference_label}'
            result = result[0]
            row = {
                'Protein': label,
                'Pruned_RMSD_A': float(result['final RMSD']),
                'Pruned_pairs': int(len(result['final match atoms'])),
                'All_pairs_RMSD_A': float(result['full RMSD']),
                'All_pairs': int(len(result['full match atoms'])),
                'Model_CA': model_ca,
                'MatchMaker_command': command,
            }
            run(session, f'close {model.atomspec}')
        table.append(row)
        print(reference_label, row, flush=True)

    all_results[reference_label] = {
        'reference_description': reference_description,
        'reference_source': reference_record['source'],
        'reference_CA': reference_ca,
        'rows': table,
    }
    run(session, f'close {reference_model.atomspec}')


def format_table(rows):
    lines = [
        f'{"Protein":<16}{"Pruned_RMSD_A":>15}{"Pruned_pairs":>14}'
        f'{"All_pairs_RMSD_A":>20}{"All_pairs":>12}{"Model_CA":>11}',
        '-' * 88,
    ]
    for row in rows:
        lines.append(
            f'{row["Protein"]:<16}'
            f'{row["Pruned_RMSD_A"]:>15.6f}'
            f'{row["Pruned_pairs"]:>14d}'
            f'{row["All_pairs_RMSD_A"]:>20.6f}'
            f'{row["All_pairs"]:>12d}'
            f'{row["Model_CA"]:>11d}'
        )
    return lines


text_lines = [
    'Pairwise structural alignment report',
    'Proteins: SatA, 2QIR, SatB, and 18 AF3 candidate structures (n = 21)',
    'Alignment: UCSF ChimeraX MatchMaker; pairing ss; cutoffDistance 2.0',
    'RMSD units: angstroms (A)',
    '',
]
for index, (reference_label, reference_description) in enumerate(references, start=1):
    text_lines.extend([
        f'=== Table {index}. Reference: {reference_label} = {reference_description} ===',
        f'Reference_source: {all_results[reference_label]["reference_source"]}',
        '',
        *format_table(all_results[reference_label]['rows']),
        '',
    ])

output_text = root / 'RMSD_to_SatA_with_2QIR.txt'
output_json = root / 'RMSD_SatA_2QIR_matchmaker_results.json'
output_text.write_text('\n'.join(text_lines).rstrip() + '\n')
output_json.write_text(json.dumps(all_results, indent=2) + '\n')
print('WROTE', output_text, flush=True)
print('WROTE', output_json, flush=True)
run(session, 'exit')
