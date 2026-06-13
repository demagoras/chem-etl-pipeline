import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Draw, Lipinski, rdDepictor

# Mute RDKit parsing warnings
RDLogger.DisableLog('rdApp.*')

# Raw data
smiles_list = [
    # Standard molecules that should pass
    "CC(=O)Oc1ccccc1C(=O)O",      # Aspirin
    "c1ccccc1O",                  # Phenol
    "CC(=O)Nc1ccccc1",            # Acetanilide
    "Nc1ccccc1",                  # Aniline
    "CCO",                        # Ethanol
    "CCN(CC)CC",                  # Triethylamine
    "C1CCOC1",                    # Tetrahydrofuran (THF)
    
    # Weirder ones
    "CCCCCC",                     # Hexane
    "C1CCCCC1",                   # Cyclohexane
    "C=C=C=C=C=C=C=C",            # Cumulene chain

    # Real drugs
    "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",  # Diazepam / Valium
    "CC(C)cc1ccc(cc1)C(C)C(=O)O",          # Ibuprofen
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",        # Caffeine

    # Duplicates
    "CC(=O)Oc1ccccc1C(=O)O",      # Aspirin
    "OCC",                        # Ethanol
    "c1(O)ccccc1",                # Phenol
    "Nc1ccccc1",                  # Aniline

    # Immediate fail
    "NOT_A_SMILES_STRING",
    "C1CCCC",                     # Unclosed ring
    "CC(=O)O_INVALID_!",
    "H2O",

    # These all fail Lipinski's rule of five
    # Cyclosporine
    "CC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N1)[C@@H]([C@H](C)C/C=C/C)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C", 
    
    # Paclitaxel / Taxol
    "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
    
    # Rifampicin
    r"C[C@H]1/C=C/C=C(\C(=O)NC2=C(C(=C3C(=C2O)C(=C(C4=C3C(=O)[C@](O4)(O/C=C/[C@@H]([C@H]([C@H]([C@@H]([C@@H]([C@@H]([C@H]1O)C)O)C)OC(=O)C)C)OC)C)C)O)O)/C=N/N5CCN(CC5)C)/C",
    
    # Erythromycin
    "CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O",

    "C[N+](C)(C)C",  # Tetramethylammonium (cation)

    # Cyclosporin A
    "CC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N1)[C@@H]([C@H](C)C/C=C/C)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C"
]

# Initialize collections
valid_mols = []
exclusion_log = []
table_rows = []
seen_canonical = set()

patterns = {
    "phenol": Chem.MolFromSmarts("c1ccccc1[OH]"),
}
all_highlights = []

print("Starting Cheminformatics Data Validation Pipeline...\n")

def truncate_smiles(smiles, max_len=20):
    """Truncates SMILES strings and adds dots if they exceed max_len."""
    return f"{smiles[:max_len]}..." if len(smiles) > max_len else smiles

for raw_s in smiles_list:
  # Sanitization and parsing
    mol = Chem.MolFromSmiles(raw_s)

    if mol is None:
        exclusion_log.append({
            "Identifier": truncate_smiles(raw_s, 30),
            "Pipeline_Stage": "SMILES Sanitization",
            "Rejection_Reason": "Invalid SMILES Syntax / Unparsable Graph"
        })
        continue

    # Deduplication
    canonical_s = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True) # Allows enantiomers to co-exist
    if canonical_s in seen_canonical:
        exclusion_log.append({
            "Identifier": truncate_smiles(raw_s, 30),
            "Pipeline_Stage": "Deduplication",
            "Rejection_Reason": f"Duplicate structure normalized to: {truncate_smiles(canonical_s)}"
        })
        continue
        
    seen_canonical.add(canonical_s)

    # Detect phenol
    has_phenol = mol.HasSubstructMatch(patterns["phenol"])

    # Detect nitrogen atoms
    n_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            n_count += 1
            charge = atom.GetFormalCharge()
            if charge != 0:
                print(f"[AUDIT LOG]: Ionized Nitrogen detected. Formal charge: {charge}.")

    # Lipinski Ro5 properties
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_donors = Lipinski.NumHDonors(mol)
    h_acceptors = Lipinski.NumHAcceptors(mol)

    violations = sum([
        mw >= 500,
        logp > 5,
        h_donors > 5,
        h_acceptors > 10
    ])

    if violations > 1:
        exclusion_log.append({
            "Identifier": truncate_smiles(canonical_s, 30),
            "Pipeline_Stage": "Physicochemical Filtering",
            "Rejection_Reason": f"Failed Lipinski's Rule of 5 ({violations} Violations)"
        })
    else:
        valid_mols.append(mol)

        # Highlight aromatic rings
        rdDepictor.Compute2DCoords(mol)
        ring_info = mol.GetRingInfo()
        highlight_atoms = set()
        
        for ring in ring_info.AtomRings():
            if all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring):
                highlight_atoms.update(ring)
        
        all_highlights.append(list(highlight_atoms))
        
        table_rows.append({
            "SMILES": canonical_s,
            "MolWt": round(mw, 2),
            "LogP": round(logp, 2),
            "H-Donors": h_donors,
            "H-Acceptors": h_acceptors,
            "Lipinski_Violations": violations,
            "Contains_Phenol": "YES" if has_phenol else "NO",
            "Nitrogen_Count": n_count
        })

df_clean = pd.DataFrame(table_rows)
df_exclusion_log = pd.DataFrame(exclusion_log)

grid_img = Draw.MolsToGridImage(
    valid_mols,
    molsPerRow=4,
    subImgSize=(250, 250),
    highlightAtomLists=all_highlights,
    legends=[f"MW: {row['MolWt']} | Phenol: {row['Contains_Phenol']} | N: {row['Nitrogen_Count']}" for _, row in df_clean.iterrows()],
    returnPNG=False
)
grid_img.save("molecules_in_grid.png")

# Output
print(f"\nPipeline Processing Complete.")
print(f"Total Source Rows: {len(smiles_list)} | Valid Compounds: {len(df_clean)} | Total Rejected Elements: {len(df_exclusion_log)}\n")

print("MOLECULE OUTPUT")
print(df_clean.to_string(index=False))

print("EXCLUSION LOG")
print(df_exclusion_log.to_string(index=False))

print("\n[SUCCESS]: Image of molecules saved as 'molecules_in_grid.png'.")
