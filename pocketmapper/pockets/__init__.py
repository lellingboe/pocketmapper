"""
Building pockets from parsed query/target records.

`pocket_fetcher.PocketFetcher` is the entry point: it dispatches each record to the builder of its
pocket method and returns one pocket_id -> Pocket dict. Each builder lives beside the primitive it
wraps (`pocket_parser`, `pisa_parser`, `pocket_calculator`), and `pocket.py` declares the shape they
all return.

Nothing is re-exported -- import the modules directly, as
`from pocketmapper.pockets.<module> import <name>`.
"""
