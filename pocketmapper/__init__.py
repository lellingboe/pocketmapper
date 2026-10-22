"""
PocketMapper: compare binding pockets across structures via sequence or structural alignment.

Only `PocketMapper` and `__version__` are exported, and no subpackage re-exports anything. Importing the
package imports `pocketmapper.search`, which imports every other module but `cli`, so submodules such as
`pocketmapper.lib.constants` resolve as attributes only as a side effect of that chain. Always import
anything else explicitly from the module that defines it: `from pocketmapper.<module> import <name>`, or
`from pocketmapper.<package>.<module> import <name>` for a module in a subpackage
(`from pocketmapper.align.step import align`).
"""

import logging

from pocketmapper.lib.constants import PACKAGE_LOGGER
from pocketmapper.search import PocketMapper

# Library convention: records are dropped here unless the host application or the CLI adds a handler.
logging.getLogger(PACKAGE_LOGGER).addHandler(logging.NullHandler())

__all__ = ["PocketMapper", "__version__"]

__version__ = "0.2.6"
