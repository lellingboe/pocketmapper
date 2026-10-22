"""
Offset tables: the UniProt coordinates of a Foldseek database's entries.
"""


def read_offset_table(path):
    """
    Read an offset table: the UniProt coordinates of every entry in a Foldseek database.

    Args:
        path (str): Path to the tab-separated table. Its header row is discarded.

    Returns:
        dict: database entry name -> region spec, as `seq_to_uniprot_map` parses.
    """
    offsets = {}
    with open(path) as f:
        next(f, None)  # header: unique_id, uniprot_domain
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            entry, domain = line.split("\t")
            offsets[entry] = domain
    return offsets


def seq_to_uniprot_map(domain):
    """
    Map each 0-indexed position of a domain onto its 1-indexed UniProt position.

    Args:
        domain (str): start1-stop1_start2-stop2_... in UniProt coords, 1-indexed and inclusive. A
            domain need not be contiguous, hence a list of regions rather than one offset.

    Returns:
        dict: 0-indexed domain position -> 1-indexed UniProt position. Increasing and injective as long
            as the regions are increasing and non-overlapping.
    """
    regions = domain.split("_")
    seq_pos_to_uniprot_pos = {}
    pos = 0
    for region in regions:
        start, stop = region.split("-")
        for i in range(int(start), int(stop) + 1):
            seq_pos_to_uniprot_pos[pos] = i
            pos += 1
    return seq_pos_to_uniprot_pos
