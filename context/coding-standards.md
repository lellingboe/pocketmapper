# Coding standards
Style, quality, and feature rules for code.

- **Comments explain why, and what breaks otherwise.** A comment that restates the code is noise; one recording a measurement, a deliberate omission or a load-bearing ordering earns its place.
- **Docstrings are Google-style** — a one-line summary, any rationale paragraphs, then `Args:` / `Returns:` / `Raises:`. Every module, class and function carries one. Keep docstrings up to date.
- **No type annotations on functions.** No function in the package carries them. Annotate
  dataclass and `NamedTuple` fields only, in `str | None` form (`Settings`, `QTRecord`).
- **`os.path`, never `pathlib`**
- **`f-strings`, never `%` or `.format()`**
- **Reuse before adding.** `lib.py` already holds the generic helpers.
- **Formatting is not a standard here.** black and flake8 run on pre-commit; their settings live in
  `pyproject.toml` and `.flake8`.
