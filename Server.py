import os
import traceback
import base64
import urllib.request
import json
import math
import re
from flask import Flask, request, jsonify
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

app = Flask(__name__)
CORS(app)

APP_VERSION = "3.6.0-KM-Birthday-Pro"

COMMON_ISOMERS = {
    "C2H6O": [
        {"name": "Ethanol (Этанол)", "smiles": "CCO"},
        {"name": "Dimethyl ether (Диметиловый эфир)", "smiles": "COC"}
    ],
    "C3H8O": [
        {"name": "1-Propanol (Пропанол-1)", "smiles": "CCCO"},
        {"name": "Isopropanol (Изопропанол)", "smiles": "CC(O)C"},
        {"name": "Methoxyethane (Мэтоксиэтан)", "smiles": "CCOC"}
    ],
    "C4H10": [
        {"name": "n-Butane (н-Бутан)", "smiles": "CCCC"},
        {"name": "Isobutane (Изобутан)", "smiles": "CC(C)C"}
    ],
    "C6H6": [
        {"name": "Benzene (Бензол)", "smiles": "c1ccccc1"},
        {"name": "Dewar Benzene (Бензол Дьюара)", "smiles": "C1=CC2C=CC12"},
        {"name": "Prismane (Призман)", "smiles": "C12C3C1C4C2C34"}
    ],
    "CH4O": [
        {"name": "Methanol (Метанол)", "smiles": "CO"}
    ],
    "C2H4O2": [
        {"name": "Acetic acid (Уксусная кислота)", "smiles": "CC(=O)O"},
        {"name": "Methyl formate (Метилформиат)", "smiles": "COC=O"}
    ],
    "C3H6O": [
        {"name": "Acetone (Ацетон)", "smiles": "CC(=O)C"},
        {"name": "Propanal (Пропаналь)", "smiles": "CCC=O"},
        {"name": "Oxetane (Оксетан)", "smiles": "C1CCO1"}
    ],
    "FrHe": [
        {"name": "Francium Helide (Комплекс «Соседи»)", "smiles": "[He]-[Fr]"}
    ],
    "HeFr": [
        {"name": "Francium Helide (Комплекс «Соседи»)", "smiles": "[He]-[Fr]"}
    ]
}

COMMON_SOLVENTS = {
    "O", "ClCCl", "C#N", "CC#N", "CCO", "CO", "C1CCOC1", "ClC(Cl)Cl", 
    "CS(=O)C", "CCCCCC", "CCCCC", "c1ccccc1", "Cc1ccccc1", "CCOCC", "CC(=O)C", "CN(C)C=O"
}

COVALENT_RADII = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
    'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
    'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.20, 'Kr': 1.16,
    'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
    'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87, 'Lu': 1.87,
    'Hf': 1.75, 'Ta': 1.70, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48,
    'Fr': 2.60, 'Ra': 2.21
}

MAX_VALENCE = {
    'H': 1, 'He': 1, 'Li': 1, 'Be': 2, 'B': 3, 'C': 4, 'N': 4, 'O': 2, 'F': 1,
    'Na': 1, 'Mg': 2, 'Al': 3, 'Si': 4, 'P': 5, 'S': 6, 'Cl': 1,
    'K': 1, 'Ca': 2, 'Br': 1, 'I': 1, 'Fr': 1
}

TRANSITION_METALS = {
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Fr', 'He'
}

def export_mol_block_safe(mol) -> str:
    try:
        return Chem.MolToMolBlock(mol, forceV3000=False, kekulize=True)
    except Exception:
        return Chem.MolToMolBlock(mol, forceV3000=False, kekulize=False)

def perceive_bonds_by_covalent_radii(mol, connect_distant=False):
    conf = mol.GetConformer()
    rw_mol = Chem.RWMol(mol)
    num_atoms = rw_mol.GetNumAtoms()
    positions = [conf.GetAtomPosition(i) for i in range(num_atoms)]
    symbols = [rw_mol.GetAtomWithIdx(i).GetSymbol() for i in range(num_atoms)]

    # Проход 1: Стандартное ковалентное расстояние
    for i in range(num_atoms):
        r_i = COVALENT_RADII.get(symbols[i], 1.4)
        pos_i = positions[i]
        for j in range(i + 1, num_atoms):
            r_j = COVALENT_RADII.get(symbols[j], 1.4)
            pos_j = positions[j]
            dist = math.sqrt((pos_i.x - pos_j.x)**2 + (pos_i.y - pos_j.y)**2 + (pos_i.z - pos_j.z)**2)

            if 0.35 <= dist <= (r_i + r_j + 0.55):
                if rw_mol.GetBondBetweenAtoms(i, j) is None:
                    max_i = MAX_VALENCE.get(symbols[i], 6)
                    max_j = MAX_VALENCE.get(symbols[j], 6)
                    if rw_mol.GetAtomWithIdx(i).GetDegree() < max_i and rw_mol.GetAtomWithIdx(j).GetDegree() < max_j:
                        rw_mol.AddBond(i, j, Chem.BondType.SINGLE)

    # Проход 2: Если пользователь расставил атомы свободно на 3D-сетке (connect_distant=True),
    # соединяем одиночные атомы с ближайшими соседями с учетом валентности
    if connect_distant and num_atoms > 1:
        pairs = []
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                pos_i, pos_j = positions[i], positions[j]
                dist = math.sqrt((pos_i.x - pos_j.x)**2 + (pos_i.y - pos_j.y)**2 + (pos_i.z - pos_j.z)**2)
                pairs.append((dist, i, j))
        pairs.sort(key=lambda x: x[0])

        # Сначала обеспечиваем связность для атомов со степенью 0
        for dist, i, j in pairs:
            if dist > 8.0:
                continue
            atom_i = rw_mol.GetAtomWithIdx(i)
            atom_j = rw_mol.GetAtomWithIdx(j)
            max_i = MAX_VALENCE.get(symbols[i], 6)
            max_j = MAX_VALENCE.get(symbols[j], 6)
            if rw_mol.GetBondBetweenAtoms(i, j) is None:
                if (atom_i.GetDegree() == 0 or atom_j.GetDegree() == 0) and atom_i.GetDegree() < max_i and atom_j.GetDegree() < max_j:
                    rw_mol.AddBond(i, j, Chem.BondType.SINGLE)

        # Если все еще есть несвязанные фрагменты, соединяем ближайшие доступные центры
        frags = Chem.GetMolFrags(rw_mol.GetMol())
        while len(frags) > 1:
            connected_any = False
            frag_sets = [set(f) for f in frags]
            for dist, i, j in pairs:
                frag_i = next(idx for idx, s in enumerate(frag_sets) if i in s)
                frag_j = next(idx for idx, s in enumerate(frag_sets) if j in s)
                if frag_i != frag_j:
                    atom_i = rw_mol.GetAtomWithIdx(i)
                    atom_j = rw_mol.GetAtomWithIdx(j)
                    max_i = MAX_VALENCE.get(symbols[i], 6)
                    max_j = MAX_VALENCE.get(symbols[j], 6)
                    if atom_i.GetDegree() < max_i and atom_j.GetDegree() < max_j:
                        rw_mol.AddBond(i, j, Chem.BondType.SINGLE)
                        connected_any = True
                        break
            if not connected_any:
                break
            frags = Chem.GetMolFrags(rw_mol.GetMol())

    return rw_mol.GetMol()

def generate_metal_complex_sdf(smiles_str: str) -> str:
    if smiles_str in ["[He]-[Fr]", "[Fr]-[He]", "[He][Fr]", "[Fr][He]"]:
        return (
            "Neighbors_He_Fr\n  ACE_3D    00000000003D\n\n"
            "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
            "   -1.1500    0.6000    0.0000 He  0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    1.1500    0.6000    0.0000 Fr  0  0  0  0  0  0  0  0  0  0  0  0\n"
            "  1  2  1  0  0  0  0\n"
            "M  END\n$$$$\n"
        )

    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is not None:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                pass
            return export_mol_block_safe(mol)
    except Exception:
        pass

    mol = Chem.MolFromSmiles(smiles_str, sanitize=False)
    if mol is None:
        raise ValueError("Invalid SMILES string provided.")
    
    metal_idxs = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym in ['Pt', 'Pd', 'Fe', 'Ru', 'Rh', 'Ir', 'Ni', 'Fr', 'He']:
            metal_idxs[atom.GetIdx()] = sym
            atom.SetAtomicNum(50)
            atom.SetFormalCharge(0)
            
        if sym == 'N' and atom.GetDegree() >= 4: atom.SetFormalCharge(1)
        if sym == 'O' and atom.GetDegree() >= 3: atom.SetFormalCharge(1)
        if sym == 'P' and atom.GetDegree() >= 4: atom.SetFormalCharge(1)

    mol.UpdatePropertyCache(strict=False)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        except Exception:
            pass 
    
    try:
        mol = Chem.AddHs(mol)
    except Exception:
        pass

    try:
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True 
        params.maxAttempts = 1000
        embed_success = AllChem.EmbedMolecule(mol, params)
        if embed_success == 0:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        else:
            AllChem.EmbedMolecule(mol, useRandomCoords=True, maxAttempts=2000)
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        pass

    for idx, sym in metal_idxs.items():
        atom = mol.GetAtomWithIdx(idx)
        if sym == 'Pt': atom.SetAtomicNum(78)
        elif sym == 'Pd': atom.SetAtomicNum(46)
        elif sym == 'Fe': atom.SetAtomicNum(26)
        elif sym == 'Fr': atom.SetAtomicNum(87)
        elif sym == 'He': atom.SetAtomicNum(2)

    return export_mol_block_safe(mol)

@app.route('/api/version', methods=['GET'])
def get_version():
    return jsonify({
        'version': APP_VERSION,
        'edition': 'Kirill Mikhailovich Lyapishev Special Edition',
        'status': 'online'
    })

@app.route('/api/meme_image', methods=['GET'])
def get_meme_image():
    for filename in ['image.png', 'image.jpg', 'image.jpeg', 'image.webp']:
        if os.path.exists(filename):
            with open(filename, 'rb') as f:
                b64 = base64.b64encode(f.read()).decode('utf-8')
            mime = 'image/png' if filename.endswith('.png') else 'image/jpeg'
            return jsonify({'image': f"data:{mime};base64,{b64}", 'status': 'success'})
    return jsonify({'error': 'image.png not found in directory', 'status': 'missing'}), 404

@app.route('/api/generate_sdf', methods=['POST'])
def generate_sdf():
    data = request.get_json(force=True)
    smiles = data.get('smiles', '').strip()
    if not smiles:
        return jsonify({'error': 'No SMILES provided'}), 400
    try:
        sdf_block = generate_metal_complex_sdf(smiles)
        return jsonify({'sdf': sdf_block, 'status': 'success'})
    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e), 'status': 'failed'}), 400

@app.route('/api/generate_2d', methods=['POST'])
def generate_2d():
    data = request.get_json(force=True)
    smiles = data.get('smiles', '').strip()
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string', 'status': 'failed'}), 400
            
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'N' and atom.GetDegree() >= 4: atom.SetFormalCharge(1)
            if atom.GetSymbol() == 'O' and atom.GetDegree() >= 3: atom.SetFormalCharge(1)
            if atom.GetSymbol() == 'P' and atom.GetDegree() >= 4: atom.SetFormalCharge(1)
                
        mol.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            pass
            
        AllChem.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(500, 400)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        b64_svg = base64.b64encode(svg.encode('utf-8')).decode('utf-8')
        return jsonify({'image': f"data:image/svg+xml;base64,{b64_svg}", 'status': 'success'})
    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e), 'status': 'failed'}), 500

@app.route('/api/sdf_to_smiles', methods=['POST'])
def sdf_to_smiles():
    data = request.get_json(force=True)
    sdf = data.get('sdf', '')
    if not sdf:
        return jsonify({'error': 'No SDF provided'}), 400
    try:
        mol = Chem.MolFromMolBlock(sdf, sanitize=False)
        if mol is None:
            return jsonify({'error': 'Invalid SDF structure'}), 400
        mol.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            pass
        smiles = Chem.MolToSmiles(mol)
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        weight = round(Chem.rdMolDescriptors.CalcExactMolWt(mol), 3)
        return jsonify({'smiles': smiles, 'formula': formula, 'weight': weight, 'status': 'success'})
    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e), 'status': 'failed'}), 500

@app.route('/api/clean_3d', methods=['POST'])
def clean_3d():
    data = request.get_json(force=True)
    sdf = data.get('sdf', '')
    if not sdf:
        return jsonify({'error': 'No SDF provided'}), 400
    try:
        mol = Chem.MolFromMolBlock(sdf, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            pass
        
        try:
            ff = AllChem.UFFGetMoleculeForceField(mol)
            if ff is not None:
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() in TRANSITION_METALS:
                        ff.AddFixedPoint(atom.GetIdx())
                ff.Minimize(maxIts=500)
        except Exception:
            pass

        new_sdf = export_mol_block_safe(mol)
        return jsonify({'sdf': new_sdf, 'status': 'success'})
    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e), 'status': 'failed'}), 500

@app.route('/api/auto_bond', methods=['POST'])
def auto_bond():
    data = request.get_json(force=True)
    elements = data.get('elements', [])
    existing_sdf = data.get('sdf', '').strip()

    if set(elements) == {'He', 'Fr'} and len(elements) == 2:
        sdf_he_fr = generate_metal_complex_sdf("[He]-[Fr]")
        return jsonify({
            'mode': 'isomers',
            'formula': 'HeFr',
            'isomers': [{'name': 'Francium Helide («Соседи»)', 'smiles': '[He]-[Fr]', 'sdf': sdf_he_fr}],
            'status': 'success'
        })

    result_isomers = []

    # 1. ВСЕГДА строим 3D-связи напрямую из текущих координат сцены (исправляет баг неформирования связей)
    if existing_sdf:
        try:
            mol_3d = Chem.MolFromMolBlock(existing_sdf, sanitize=False, removeHs=False)
            if mol_3d is not None and mol_3d.GetNumAtoms() > 1:
                mol_3d = perceive_bonds_by_covalent_radii(mol_3d, connect_distant=True)
                mol_3d.UpdatePropertyCache(strict=False)
                try:
                    Chem.SanitizeMol(mol_3d, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
                except Exception:
                    pass

                try:
                    ff = AllChem.UFFGetMoleculeForceField(mol_3d)
                    if ff is not None:
                        for atom in mol_3d.GetAtoms():
                            if atom.GetSymbol() in TRANSITION_METALS:
                                ff.AddFixedPoint(atom.GetIdx())
                        ff.Minimize(maxIts=400)
                except Exception:
                    pass

                connected_sdf = export_mol_block_safe(mol_3d)
                derived_sm = Chem.MolToSmiles(mol_3d)
                result_isomers.append({
                    'name': 'Собранная 3D-структура (Auto-Bonded)',
                    'smiles': derived_sm,
                    'sdf': connected_sdf
                })
        except Exception:
            traceback.print_exc()

    # 2. Дополнительно ищем известные химические изомеры по брутто-формуле
    counts = {}
    for el in elements:
        counts[el] = counts.get(el, 0) + 1

    formula_parts = []
    if 'C' in counts:
        c_cnt = counts['C']
        formula_parts.append(f"C{c_cnt if c_cnt > 1 else ''}")
    if 'H' in counts:
        h_cnt = counts['H']
        formula_parts.append(f"H{h_cnt if h_cnt > 1 else ''}")
    for el in sorted(counts.keys()):
        if el not in ('C', 'H'):
            formula_parts.append(f"{el}{counts[el] if counts[el] > 1 else ''}")
    formula = "".join(formula_parts)

    candidates = list(COMMON_ISOMERS.get(formula, []))

    # Если пользователь набросал только тяжелые атомы без водорода (например C, C, O),
    # проверим также насыщенную водородом органическую формулу
    if not candidates and 'C' in counts and 'H' not in counts:
        c_num = counts['C']
        n_num = counts.get('N', 0)
        hal_num = sum(counts.get(x, 0) for x in ['F', 'Cl', 'Br', 'I'])
        for unsat in [0, 1, 2]:
            h_est = (2 * c_num + 2 + n_num - hal_num) - (2 * unsat)
            if h_est > 0:
                sat_parts = [f"C{c_num if c_num > 1 else ''}", f"H{h_est if h_est > 1 else ''}"]
                for el in sorted(counts.keys()):
                    if el not in ('C', 'H'):
                        sat_parts.append(f"{el}{counts[el] if counts[el] > 1 else ''}")
                sat_formula = "".join(sat_parts)
                if sat_formula in COMMON_ISOMERS:
                    candidates.extend(COMMON_ISOMERS[sat_formula])

    if not candidates and len(elements) <= 25:
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/{formula}/property/CanonicalSMILES,IUPACName/JSON"
            req = urllib.request.Request(url, headers={'User-Agent': 'ChemistryEngine/1.0'})
            with urllib.request.urlopen(req, timeout=3) as resp:
                pdata = json.loads(resp.read().decode())
                props = pdata.get('PropertyTable', {}).get('Properties', [])
                seen = set()
                for p in props[:6]:
                    sm = p.get('CanonicalSMILES')
                    if sm and sm not in seen:
                        seen.add(sm)
                        candidates.append({
                            'name': p.get('IUPACName', sm),
                            'smiles': sm
                        })
        except Exception:
            pass

    for cand in candidates:
        try:
            sdf = generate_metal_complex_sdf(cand['smiles'])
            result_isomers.append({
                'name': cand['name'],
                'smiles': cand['smiles'],
                'sdf': sdf
            })
        except Exception:
            pass

    if not result_isomers:
        return jsonify({'error': f'Не удалось сформировать связи для: {formula}'}), 400

    return jsonify({
        'mode': 'isomers',
        'formula': formula,
        'isomers': result_isomers,
        'status': 'success'
    })

@app.route('/api/manage_hydrogens', methods=['POST'])
def manage_hydrogens():
    data = request.get_json(force=True)
    sdf = data.get('sdf', '')
    action = data.get('action', 'add')

    if not sdf:
        return jsonify({'error': 'No SDF provided'}), 400
    try:
        mol = Chem.MolFromMolBlock(sdf, sanitize=False, removeHs=False)
        mol.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            pass

        if action == 'add':
            mol = Chem.AddHs(mol, addCoords=True)
            try:
                ff = AllChem.UFFGetMoleculeForceField(mol)
                if ff is not None:
                    for atom in mol.GetAtoms():
                        if atom.GetSymbol() != 'H':
                            ff.AddFixedPoint(atom.GetIdx())
                    ff.Minimize(maxIts=500)
            except Exception:
                pass
        elif action == 'strip_all':
            mol = Chem.RemoveHs(mol, sanitize=False)
        elif action == 'strip_nonpolar':
            rw_mol = Chem.RWMol(mol)
            to_remove = []
            for atom in rw_mol.GetAtoms():
                if atom.GetSymbol() == 'H':
                    neighbors = atom.GetNeighbors()
                    if len(neighbors) == 1 and neighbors[0].GetSymbol() == 'C':
                        to_remove.append(atom.GetIdx())
            for idx in sorted(to_remove, reverse=True):
                rw_mol.RemoveAtom(idx)
            mol = rw_mol.GetMol()
        elif action == 'strip_polar':
            rw_mol = Chem.RWMol(mol)
            to_remove = []
            for atom in rw_mol.GetAtoms():
                if atom.GetSymbol() == 'H':
                    neighbors = atom.GetNeighbors()
                    if len(neighbors) == 1 and neighbors[0].GetSymbol() != 'C':
                        to_remove.append(atom.GetIdx())
            for idx in sorted(to_remove, reverse=True):
                rw_mol.RemoveAtom(idx)
            mol = rw_mol.GetMol()

        out_sdf = export_mol_block_safe(mol)
        return jsonify({'sdf': out_sdf, 'status': 'success'})
    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e), 'status': 'failed'}), 500

@app.route('/api/calculate_charges', methods=['POST'])
def calculate_charges():
    data = request.get_json(force=True)
    sdf = data.get('sdf', '')
    if not sdf:
        return jsonify({'error': 'No SDF provided'}), 400
    try:
        mol = Chem.MolFromMolBlock(sdf, sanitize=False, removeHs=False)
        mol.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            pass
        
        try:
            AllChem.ComputeGasteigerCharges(mol)
        except Exception:
            pass

        charges = []
        for atom in mol.GetAtoms():
            val = 0.0
            if atom.HasProp('_GasteigerCharge'):
                try:
                    c = float(atom.GetProp('_GasteigerCharge'))
                    if not math.isnan(c) and not math.isinf(c):
                        val = round(c, 3)
                except Exception:
                    val = 0.0
            charges.append(val)

        return jsonify({'charges': charges, 'status': 'success'})
    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e), 'status': 'failed'}), 500

def parse_cif_to_sdf(cif_text: str, strip_solvents: bool = False) -> str:
    a = b = c = 1.0
    alpha = beta = gamma = 90.0

    for line in cif_text.splitlines():
        line = line.strip()
        if line.startswith('_cell_length_a'): a = float(re.findall(r"[\d.]+", line)[0])
        elif line.startswith('_cell_length_b'): b = float(re.findall(r"[\d.]+", line)[0])
        elif line.startswith('_cell_length_c'): c = float(re.findall(r"[\d.]+", line)[0])
        elif line.startswith('_cell_angle_alpha'): alpha = float(re.findall(r"[\d.]+", line)[0])
        elif line.startswith('_cell_angle_beta'): beta = float(re.findall(r"[\d.]+", line)[0])
        elif line.startswith('_cell_angle_gamma'): gamma = float(re.findall(r"[\d.]+", line)[0])

    ar = math.radians(alpha)
    br = math.radians(beta)
    gr = math.radians(gamma)
    v = math.sqrt(max(0.0001, 1 - math.cos(ar)**2 - math.cos(br)**2 - math.cos(gr)**2 + 2*math.cos(ar)*math.cos(br)*math.cos(gr)))

    atoms = []
    in_atom_loop = False
    headers = []

    for line in cif_text.splitlines():
        line = line.strip()
        if not line or line.startswith('#'): continue
        if line.startswith('loop_'):
            in_atom_loop = False
            headers = []
            continue
        if line.startswith('_atom_site_'):
            in_atom_loop = True
            headers.append(line)
            continue
        if in_atom_loop and not line.startswith('_'):
            parts = line.split()
            if len(parts) >= len(headers) and len(headers) >= 4:
                row = dict(zip(headers, parts))
                elem = row.get('_atom_site_type_symbol', '')
                if not elem:
                    elem = re.sub(r'[^a-zA-Z]', '', row.get('_atom_site_label', 'C'))
                else:
                    elem = re.sub(r'[^a-zA-Z]', '', elem)
                elem = elem.capitalize()

                if '_atom_site_Cartn_x' in row:
                    x = float(re.findall(r"[-+]?[\d.]+", row['_atom_site_Cartn_x'])[0])
                    y = float(re.findall(r"[-+]?[\d.]+", row['_atom_site_Cartn_y'])[0])
                    z = float(re.findall(r"[-+]?[\d.]+", row['_atom_site_Cartn_z'])[0])
                elif '_atom_site_fract_x' in row:
                    fx = float(re.findall(r"[-+]?[\d.]+", row['_atom_site_fract_x'])[0])
                    fy = float(re.findall(r"[-+]?[\d.]+", row['_atom_site_fract_y'])[0])
                    fz = float(re.findall(r"[-+]?[\d.]+", row['_atom_site_fract_z'])[0])
                    x = a * fx + b * math.cos(gr) * fy + c * math.cos(br) * fz
                    y = b * math.sin(gr) * fy + c * ((math.cos(ar) - math.cos(br)*math.cos(gr)) / math.sin(gr)) * fz
                    z = c * (v / math.sin(gr)) * fz
                else:
                    continue

                atoms.append((elem, x, y, z))

    if not atoms:
        raise ValueError("Could not extract atom coordinates from CIF file.")

    rw_mol = Chem.RWMol()
    conf = Chem.Conformer()
    for i, (el, x, y, z) in enumerate(atoms):
        idx = rw_mol.AddAtom(Chem.Atom(el if len(el) <= 2 else 'C'))
        conf.SetAtomPosition(idx, (x, y, z))
    rw_mol.AddConformer(conf)
    mol = rw_mol.GetMol()

    mol = perceive_bonds_by_covalent_radii(mol, connect_distant=False)

    if strip_solvents:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        retained = []
        for f in frags:
            sm = Chem.MolToSmiles(f)
            if sm not in COMMON_SOLVENTS and f.GetNumAtoms() > 3:
                retained.append(f)
        if retained:
            mol = retained[0]
            for extra in retained[1:]:
                mol = Chem.CombineMols(mol, extra)

    return export_mol_block_safe(mol)

@app.route('/api/import_crystallography', methods=['POST'])
def import_crystallography():
    data = request.get_json(force=True)
    raw_text = data.get('content', '')
    ext = data.get('format', 'cif').lower()
    strip_solvents = data.get('strip_solvents', False)

    if not raw_text:
        return jsonify({'error': 'No file content received'}), 400

    try:
        if ext in ['pdb', 'ent']:
            mol = Chem.MolFromPDBBlock(raw_text, sanitize=False, removeHs=False)
            if mol.GetNumBonds() == 0:
                mol = perceive_bonds_by_covalent_radii(mol, connect_distant=False)
            if strip_solvents:
                frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
                retained = [f for f in frags if Chem.MolToSmiles(f) not in COMMON_SOLVENTS and f.GetNumAtoms() > 3]
                if retained:
                    mol = retained[0]
                    for extra in retained[1:]:
                        mol = Chem.CombineMols(mol, extra)
            sdf_out = export_mol_block_safe(mol)
        else:
            sdf_out = parse_cif_to_sdf(raw_text, strip_solvents)

        return jsonify({'sdf': sdf_out, 'status': 'success'})
    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e), 'status': 'failed'}), 500

if __name__ == '__main__':
    print("🚀 Chemistry Supercomputer Server [Kirill Mikhailovich Edition] running on http://127.0.0.1:5000")
    app.run(host='127.0.0.1', port=5000, debug=True)
