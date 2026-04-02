from collections import defaultdict

# TODO keep phospho information
SINGLE_AA_CODE = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
    "SEP": "S",  # phospho
    "TPO": "T",  # phospho
    "PTR": "Y",  # phospho
    "MSE": "M",  # selenomethionine
}
TRIPLE_AA_CODE = defaultdict(list)
for k, v in SINGLE_AA_CODE.items():
    TRIPLE_AA_CODE[v].append(k)
VDW_RADII = {"C": 1.88, "N": 1.64, "O": 1.46, "S": 1.77, "P": 1.87, "H": 1.0}
# https://www.cgl.ucsf.edu/chimerax/docs/user/commands/clashes.html
