# Coding standards

## General
- **Choose the solution that is easiest to maintain.** Fewer special cases and fewer places to change
  beat cleverness.
- **Update docs, docstrings and comments in the same change as the code.** Stale documentation causes
  errors.

## Comments and docstrings
- **Local.** Say what this code does and what it depends on or must precede. Never say who calls it or
  how it is used elsewhere; that goes stale.
- **Concise.** A comment should not get in the way of reading the code.
- **Every module, class and function has a docstring.**
- **Google style**, with the summary on the line after the opening `"""`:
  - one-line summary;
  - paragraphs only for side effects or behaviour the caller would not expect;
  - `Args:` as `name (type): Description.`, adding `, optional` when it has a default; alternatives
    as `str or None`, never `str | None`;
  - `Returns:` / `Yields:` as `type: Description.`; `Raises:` as `ExceptionType: When.`

## Python
- **Type annotations only on dataclasses and NamedTuples.** No function carries them.
- **`os.path`, never `pathlib`.**
- **f-strings, never `%` or `.format()`**, including in logger calls.
- **Absolute imports** (`from pocketmapper.x import y`).
- **No leading-underscore names**; everything is public and documented.
- **`except Exception` only to clean up and re-raise, or to `logger.exception` and continue a
  per-item loop.** Never a bare `except:`.

## Logging and errors
Mechanism: `project-overview.md` "Logging and errors".
- **`logger = logging.getLogger(__name__)`**; never `logging.info(...)` etc. on the root logger.
- **`extra` is always named `log_extra`.** A single-stage class sets `self.log_extra` once in `__init__`
  and never mutates it; elsewhere build it per function. A declared `stage` is the Title Case step name.
- **Fail with `logger.critical(...)` then `raise PocketMapperError(...)`.** No `exit()`/`sys.exit()`
  outside `cli()`.
- **A writer into a configurable path makes its own directory**: `lib.make_dir` in the pipeline,
  `os.makedirs` in components.

## Formatting
- **ruff, configured in `pyproject.toml`.** pre-commit runs `ruff check --fix` and `ruff format`; CI
  checks both.
- **120 columns, comments and docstrings included.** `E501` is off and ruff does not wrap comments, so
  wrap them by hand.
