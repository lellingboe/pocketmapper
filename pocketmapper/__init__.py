"""
PocketMapper: compare binding pockets across structures via sequence or structural alignment.

Only `main` and `__version__` are exported. The submodules are reachable as `pocketmapper.lib`
etc. only as a side effect of `pocketmapper.pocketmapper` importing them, so always import them
explicitly -- `from pocketmapper.<module> import <name>`.
"""

from pocketmapper.pocketmapper import main

__version__ = "0.1.0"
