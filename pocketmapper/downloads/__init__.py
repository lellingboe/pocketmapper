"""
Fetching files from external services over HTTP.

Holds the PISA downloader. The retry and atomic-write helpers it shares with the structure downloader are
in `lib.download`.

Nothing is re-exported -- import the modules directly, as
`from pocketmapper.downloads.<module> import <name>`.
"""
