"""
Fetching files from external services over HTTP.

Holds the urllib-based downloaders and the retry and atomic-write helpers they share. Fetching by
shelling out to an external tool is not in scope, so the Foldseek database download is not here.

Nothing is re-exported -- import the modules directly, as
`from pocketmapper.downloads.<module> import <name>`.
"""
