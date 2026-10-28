"""
Unit tests for `pocketmapper.compare.compare`: the check that the pockets file covers the alignment.
"""

import re

import pandas as pd
import pytest

from pocketmapper.compare.compare import check_coverage
from pocketmapper.lib.exceptions import PocketMapperError

# A human-domains hit name: not PDB-style, so its pocket would be synthesised
DB_HIT = "A0A024RBG1_chainsaw+pfam_001"


@pytest.fixture
def names(processor):
    """The preprocess name of each 4Q5J chain the cases use, keyed by chain, as parse computes them."""
    return {chain: processor.parse_individual_qt(f"4Q5J:{chain}", "auto")[0].preprocess_name for chain in "ABE"}


def alignment(pairs):
    """An alignment table holding only the query and target columns, a row per (query, target) pair."""
    return pd.DataFrame(pairs, columns=["query", "target"])


# (pairs as chains or DB_HIT, chains listed as chain -> pockets, synthesise, chain the message names or None)
COVERAGE_TABLE = [
    pytest.param([("A", "B")], {"A": ["4Q5J:A_E"], "B": ["4Q5J:B_F"]}, False, None, id="covered"),
    pytest.param([("A", "B")], {"A": ["4Q5J:A_E"]}, False, "B", id="target-unlisted"),
    pytest.param([("A", "B")], {"B": ["4Q5J:B_F"]}, False, "A", id="query-unlisted"),
    pytest.param([("A", "B"), ("A", "E")], {"A": ["4Q5J:A_E"], "E": []}, False, "B", id="one-of-two-unlisted"),
    pytest.param([("A", "B")], {"A": ["4Q5J:A_E"], "B": []}, False, None, id="listed-without-pocket"),
    pytest.param([("A", DB_HIT)], {"A": ["4Q5J:A_E"]}, True, None, id="synthesised-target-exempt"),
    pytest.param([("A", DB_HIT)], {"A": ["4Q5J:A_E"]}, False, DB_HIT, id="unsynthesised-target-checked"),
    pytest.param([("B", DB_HIT)], {"A": ["4Q5J:A_E"]}, True, "B", id="synthesised-query-checked"),
]


@pytest.mark.parametrize("pairs, listed, synthesise, missing", COVERAGE_TABLE)
def test_check_coverage(names, pairs, listed, synthesise, missing):
    """Every aligned name must be a key of `chains`, but a synthesised target; a listed chain may have no pocket."""
    resolved = {**names, DB_HIT: DB_HIT}
    alignment_df = alignment([(resolved[query], resolved[target]) for query, target in pairs])
    chain_pockets = {names[chain]: pocket_ids for chain, pocket_ids in listed.items()}
    if missing is None:
        assert check_coverage(alignment_df, chain_pockets, synthesise, "pockets.json") is None
    else:
        message = f"pockets.json does not cover the alignment: no entry for {resolved[missing]}. Run pockets"
        with pytest.raises(PocketMapperError, match=re.escape(message)):
            check_coverage(alignment_df, chain_pockets, synthesise, "pockets.json")
