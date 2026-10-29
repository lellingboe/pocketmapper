"""
Helpers shared across PocketMapper.

`settings`, `entries` and `scopes` are pipeline-level: the `Settings` record and job-file layering, the
`entries.tsv` table every step reads and marks, and the step and run scopes. The other modules take plain
values and return plain values, or do one thing to the filesystem or the package logger. Workflow logic
belongs in the component modules rather than here, and nothing here imports a step package.

Nothing is re-exported -- import the modules directly, as `from pocketmapper.lib.<module> import <name>`.
"""
