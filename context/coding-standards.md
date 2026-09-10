# Coding standards

## General
- **Write quality code!** Take the extra time to do things the right way.
- **Find solutions which are easy to maintain.** Choose solutions which minimize codebase entropy.
- **Keep documentation, docstrings, and comments up-to-date.** Stale documentation only leads to more errors.

## Comments
- **Concise**
  - To avoid getting in the way of just reading the code
- **Local**
  - Reference what a constant is or what a block of code does
  - Do not reference how it is being used elsewhere, this is likely to become stale.

## Docstrings
- **Google-style**
  - A one-line summary
  - Paragraphs with any details important when calling the function
  - `Args:` / `Returns:` / `Raises:`
- **Function docstrings are local**
  - Limit information to what the function does
  - Do not reference how it is used elsewhere, this is likely to become stale.

## Features
- **Limit type annotations to dataclasses.** No function in the package carries them.
- **`os.path`, never `pathlib`**
- **`f-strings`, never `%` or `.format()`**
- **Formatting is managed by black and flake8.** Their settings live in `pyproject.toml` and `.flake8`.
