"""
PocketMapper: compare binding pockets across structures via sequence or structural alignment.

Only `PocketMapper` and `__version__` are exported. The submodules are reachable as
`pocketmapper.lib` etc. only as a side effect of `pocketmapper.pocketmapper` importing them, so
always import anything else explicitly -- `from pocketmapper.<module> import <name>`, or
`from pocketmapper.downloads.<module> import <name>` for the download components.
"""

import logging

from pocketmapper.constants import PACKAGE_LOGGER
from pocketmapper.pocketmapper import PocketMapper

# Library convention: records are dropped here unless the host application or the CLI adds a handler.
logging.getLogger(PACKAGE_LOGGER).addHandler(logging.NullHandler())

__all__ = ["PocketMapper", "__version__"]

__version__ = "0.2.2"
