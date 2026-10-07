"""
The pipeline's steps, one module each, in the order `search` runs them.

Each step reads its input files, does its work and writes its output files: it takes explicit
paths and values, never a `Settings`, and keeps no state between calls. The files they hand each
other are described in `pocketmapper.records`.
"""
