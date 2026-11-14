# Publishing `gsua-csb` to PyPI

The package name **`gsua-csb` is free on PyPI** (checked 2026-11-07), so the first upload claims it.

## Already done (v0.5.0)

- Built `python/dist/gsua_csb-0.5.0.tar.gz` (sdist) and `python/dist/gsua_csb-0.5.0-py3-none-any.whl`
  (2026-11-14). 0.4.0 was never uploaded; 0.5.0 supersedes it (it is reproducible from commit
  `81276f0` if ever needed).
- Smoke test of the built wheel (imported from the archive itself): `import gsua_csb` reports `0.5.0`
  and `sensitivity_analysis(..., method="shapley")` reproduces the Ishigami Shapley effects.
- **Not yet run for 0.5.0:** `twine check dist/*` (twine is not installed in the base env — see below).
  Run it before uploading.

These are the exact files to upload. A rebuild (`python -m build`) is reproducible, but **do not change
the version** — PyPI refuses to overwrite an already-published version.

## ⚠ Install twine in isolation, NOT in the conda base env

`twine` (7.x) requires `packaging >= 26.1`, but this machine's `streamlit` pins `packaging < 25`.
Installing twine into the base env breaks streamlit (and vice versa). Use an isolated tool:

```bash
pipx install twine          # recommended; fully isolated
# or a throwaway venv:
python -m venv ~/twine-venv && ~/twine-venv/Scripts/pip install twine
```

(If you ran `pip install twine` in base earlier, undo it with `pip uninstall twine` and
`pip install "packaging<25"` to restore streamlit — already done in this session.)

## One-time setup

1. Create a PyPI account at https://pypi.org/account/register/ and enable 2FA.
2. Create an API token at https://pypi.org/manage/account/token/ — scope "Entire account" for the
   first upload (you can create a project-scoped token afterwards). It looks like `pypi-AgEI...`.

## (Recommended for a first publish) TestPyPI dry run

1. Register at https://test.pypi.org/ and make a separate token there.
2. Upload and install from TestPyPI in a clean venv:

```bash
cd python
twine upload --repository testpypi dist/gsua_csb-0.5.0*
python -m venv /tmp/t && /tmp/t/Scripts/pip install \
  -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ gsua-csb
/tmp/t/Scripts/python -c "import gsua_csb; print(gsua_csb.__version__)"
```

(The `--extra-index-url` is needed because numpy/scipy/sklearn aren't on TestPyPI.)

## Real upload

```bash
cd python
twine upload dist/gsua_csb-0.5.0*
# username: __token__
# password: <your pypi-... token>
```

Non-interactive alternative:

```bash
TWINE_USERNAME=__token__ TWINE_PASSWORD=pypi-AgEI... twine upload dist/gsua_csb-0.5.0*
```

## Verify

```bash
python -m venv /tmp/v && /tmp/v/Scripts/pip install gsua-csb
/tmp/v/Scripts/python -c "import gsua_csb; print(gsua_csb.__version__)"   # -> 0.5.0
```

## Future releases

Bump the version in `pyproject.toml` and `src/gsua_csb/__init__.py`, then:

```bash
cd python && rm -rf dist && python -m build && twine check dist/* && twine upload dist/*
```
