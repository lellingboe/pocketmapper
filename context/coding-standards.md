# Coding standards

A rule earns a line here only when breaking it fails silently rather than loudly. What a docstring
or comment at the code site already says appears as a pointer to that site, never as a second copy.

- **Comments explain why, and what breaks otherwise.** Density is deliberately high — see
  `constants.py`, where nearly every constant carries a rationale above it. A comment that restates
  the code is noise; one recording a measurement, a deliberate omission or a load-bearing ordering
  earns its place.
- **Docstrings are Google-style** — `Args:` / `Returns:` / `Raises:`. Match that rather than the
  reST outlier in `structure_preprocessor.py` or the section-less prose in `lib.py`, and keep
  `Returns:` honest — several already claim `None` and return a value.
- **No type annotations on functions.** No function in the package carries them. Annotate
  dataclass and `NamedTuple` fields only, in `str | None` form (`Settings`, `QTRecord`).
- **`os.path`, never `pathlib`; f-strings, never `%` or `.format()`** — inside log calls too.
  Currently no exceptions to either.
- **Reuse before adding.** `lib.py` already holds the generic helpers, `pocket_parser.py` the
  pocket-dict construction.
- **Formatting is not a standard here.** black and flake8 run on pre-commit; their settings live in
  `pyproject.toml` and `.flake8`.
