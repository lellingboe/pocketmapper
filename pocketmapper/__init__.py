"""
PocketMapper: compare binding pockets across structures via sequence or structural alignment.

Only `main` and `__version__` are exported. The submodules are reachable as `pocketmapper.lib`
etc. only as a side effect of `pocketmapper.cli` importing them, so always import them
explicitly -- `from pocketmapper.<module> import <name>`, or
`from pocketmapper.downloads.<module> import <name>` for the download components.
"""

from pocketmapper.cli import main

__version__ = "0.1.0"
