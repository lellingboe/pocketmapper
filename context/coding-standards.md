# Coding standards
Style, quality, and feature rules for code.

- **Comments explain why, and what breaks otherwise.** A comment that restates the code is noise.
- **Docstrings are Google-style** — a one-line summary, any rationale paragraphs, then `Args:` / `Returns:` / `Raises:`. Every module, class and function carries one. Keep docstrings up to date.
- **No type annotations on functions.** No function in the package carries them. Annotate
  dataclass and `NamedTuple` fields only, in `str | None` form (`Settings`, `QTRecord`).
- **`os.path`, never `pathlib`**
- **`f-strings`, never `%` or `.format()`**
- **Formatting is managed by black and flake8.** Their settings live in `pyproject.toml` and `.flake8`.
- **Choose low codebase entropy solutions.** Reuse code when possible. Avoid sprawling code that will be hard to maintain.
