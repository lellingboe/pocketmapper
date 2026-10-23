"""
The pockets step: building pockets from parsed query/target records.

`step.py` holds the step itself. `pocket_fetcher.fetch_pockets` is the entry point for building: it
dispatches each record to the builder of its pocket method and returns one pocket_id -> Pocket dict of
the pockets built.
Each builder lives beside the primitive it wraps (`structure`, `pisa`, `vdw`), and `pocket.py` declares
the shape they all return. `pisa_downloader` fills the PISA cache that `pisa` reads.

Nothing is re-exported -- import the modules directly, as
`from pocketmapper.pockets.<module> import <name>`.
"""
