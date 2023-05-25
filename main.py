from Bio.PDB import PDBParser, PDBList
from flask import Flask, render_template, request
import py3Dmol
import requests
# from Bio.ExPASy import Prosite, get_prosite_raw
import re
import numpy as np

app = Flask(__name__)


def get_sequences(pdb_id):
    # return ['ELVMTQTPLSLPVSLGDQASISCRSSQSLVHSNGNTYLHWYLQKPGQSPKLLIYKVSNRFSGVPDRFSGSGSGTDFTLKISRVEAEDLGVYFCSQSTHVPPTFGGGTKLEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRG']
    link = f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id}'
    data = requests.get(link).json()[pdb_id.lower()]
    return [elem['sequence'] for elem in data if 'sequence' in elem]


def convert_to_regex(prosite_pattern):
    return prosite_pattern \
        .replace('{', '[^') \
        .replace('}', ']') \
        .replace('-', '') \
        .replace('X', '[A-Z]') \
        .replace('x', '[a-z]') \
        .replace('Z', '[EQ]') \
        .replace('B', '[DN]') \
        .replace('J', '[IL]') \
        .replace('O', '[KQRN]') \
        .replace('U', '[CYS]')


def get_highlighted_positions(pattern_id, protein_sequence):
    # https://prosite.expasy.org doesn't work from Russia, use vpn if you want to use it instead of reading from files
    # prosite_raw = get_prosite_raw(pattern_id)
    # record = Prosite.read(prosite_raw)
    # prosite_pattern = record.pattern

    # read from files
    prosite_raw = open('prosite/%s.txt' % pattern_id, 'r')
    prosite_pattern = prosite_raw.readlines()[4][5:-1]

    prosite_regex = convert_to_regex(prosite_pattern)
    matches = [match for match in re.finditer(prosite_regex, protein_sequence)]
    highlighted_positions = []
    [[highlighted_positions.append(i) for i in range(match.start(), match.end())] for match in matches]
    return highlighted_positions, matches, prosite_pattern


def show_3d_molecule(pdb_id, positions, chosen_chain):
    view = py3Dmol.view(query='pdb:' + pdb_id)
    view.setStyle({'cartoon': {'color': 'white'}})
    for position in positions:
        view.addStyle(
            {'chain': chosen_chain, 'resi': position, 'elem': 'C'},
            {'cartoon': {'color': '#FF0000', 'radius': 0.2}},
        )
        view.addStyle(
            {'chain': chosen_chain, 'resi': position},
            {'cartoon': {'radius': 0.2}},
        )
    return view._make_html()


def get_pseudosequence(highlighted_positions, sequence):
    return ''.join([el if (i in highlighted_positions) else '-' for i, el in enumerate(sequence)])


def get_starts_and_ends(matches):
    return [(match.start(), match.end()) for match in matches]


pdbParser = PDBParser()


def get_mhc_chains(pdb_id):
    pdb_file = PDBList().retrieve_pdb_file(pdb_id, pdir='./pdb', file_format='pdb')
    return pdbParser.get_structure(pdb_id, pdb_file).get_chains()


def get_distances_matrix_colored(mhc_chains, distance_threshold):
    chain_L_residues, chain_H_residues = None, None
    for chain in mhc_chains:
        if chain.id == 'L':
            chain_L_residues = [r for r in chain.get_residues()]
        elif chain.id == 'H':
            chain_H_residues = [r for r in chain.get_residues()]
    distances_matrix_colored = [[(0, False) for _ in range(len(chain_H_residues))] for __ in range(len(chain_L_residues))]
    for i in range(len(chain_L_residues)):
        for j in range(len(chain_H_residues)):
            try:
                point_1 = np.array(chain_L_residues[i]['CA'].get_coord())
                point_2 = np.array(chain_H_residues[j]['CA'].get_coord())
                dist = np.linalg.norm(point_1 - point_2)
                distances_matrix_colored[i][j] = (dist, 'lightgreen' if dist <= distance_threshold else 'lightgrey')
            except:
                distances_matrix_colored[i][j] = ('NA', 'lightgrey')
    return distances_matrix_colored


@app.route('/', methods=['GET', 'POST'])
def index():
    pdb_id = '1A4J'
    prosite_pattern_id = 'PS00290'
    chosen_chain = 'L'
    distance_threshold = 30

    if request.method == 'POST':
        pdb_id = request.form['pdbId'] if request.form['pdbId'] else pdb_id
        prosite_pattern_id = request.form['proSitePatternId'] if request.form['proSitePatternId'] \
            else prosite_pattern_id
        chosen_chain = request.form['selectedChain'] if request.form['selectedChain'] else chosen_chain
        distance_threshold = float(request.form['distanceThreshold']) if request.form['distanceThreshold'] \
            else distance_threshold

    sequence = get_sequences(pdb_id)[0]
    positions, matches, prosite_pattern = get_highlighted_positions(prosite_pattern_id, sequence)
    pseudosequence = get_pseudosequence(positions, sequence)
    intervals = get_starts_and_ends(matches)
    matrix = get_distances_matrix_colored(get_mhc_chains(pdb_id), distance_threshold)
    molecule = show_3d_molecule(pdb_id, positions, chosen_chain)

    return render_template('index.html',
                           pdbId=pdb_id,
                           proSitePatternId=prosite_pattern_id,
                           selectedChain=chosen_chain,
                           distanceThreshold=distance_threshold,
                           molecule=molecule,
                           pattern=prosite_pattern,
                           pseudosequence=pseudosequence,
                           intervals=intervals,
                           chains=[chain for chain in get_mhc_chains(pdb_id)],
                           matrix=matrix)


if __name__ == '__main__':
    app.run()
