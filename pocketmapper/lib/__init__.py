"""
Generic helpers shared across PocketMapper.

Nothing here knows about the pipeline, Settings, or the Pocket shape -- each function takes plain values
and returns plain values, or does one thing to the filesystem or the package logger. Workflow logic
belongs in the component modules rather than here, and nothing here imports a step package.

Nothing is re-exported -- import the modules directly, as `from pocketmapper.lib.<module> import <name>`.
"""
