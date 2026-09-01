# Coding standards

A rule earns a line here only when breaking it fails silently rather than loudly. What a docstring
or comment at the code site already says appears as a pointer to that site, never as a second copy.

- **Comments explain why, and what breaks otherwise.** Density is deliberately high — see
  `constants.py`, where nearly every constant carries a rationale above it. A comment that restates
  the code is noise; one recording a measurement, a deliberate omission or a load-bearing ordering
  earns its place.
- **Docstrings are Google-style** — a one-line summary, any rationale paragraphs, then `Args:` /
  `Returns:` / `Raises:`. Every module, class and function carries one, including trivial private
  helpers. Keep `Returns:` honest: it is wrong far more quietly than a wrong `Args:`, since a
  caller only finds out by using the value the docstring said wasn't there.
- **No type annotations on functions.** No function in the package carries them. Annotate
  dataclass and `NamedTuple` fields only, in `str | None` form (`Settings`, `QTRecord`).
- **`os.path`, never `pathlib`; f-strings, never `%` or `.format()`** — inside log calls too.
  Currently no exceptions to either.
- **Reuse before adding.** `lib.py` already holds the generic helpers, `pocket.py` the declared
  pocket shape and `pocket_parser.py` the construction of one from a structure.
- **Formatting is not a standard here.** black and flake8 run on pre-commit; their settings live in
  `pyproject.toml` and `.flake8`.
