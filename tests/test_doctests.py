"""Execute every docstring example in the package."""

import doctest
import importlib
import pkgutil

import fastatacular

# Includes the sequence-deduplication module; no `>>>` examples yet.
# Raise these when you add modules or examples; never lower them without saying why in
# the commit message. MIN_EXAMPLES is inert until the first example lands, so MIN_MODULES
# is what catches a walk that imports nothing.
MIN_MODULES = 10
MIN_EXAMPLES = 0


def _all_modules():
    mods = [fastatacular]
    for _finder, name, _ispkg in pkgutil.walk_packages(fastatacular.__path__, prefix=f"{fastatacular.__name__}."):
        mods.append(importlib.import_module(name))
    return mods


def test_docstring_examples():
    mods = _all_modules()
    assert len(mods) >= MIN_MODULES, f"module walk found only {len(mods)}"
    attempted = failed = 0
    for mod in mods:
        result = doctest.testmod(mod, verbose=False)
        attempted += result.attempted
        failed += result.failed
    assert failed == 0, f"{failed} of {attempted} docstring examples failed"
    assert attempted >= MIN_EXAMPLES, f"only {attempted} examples ran, expected >= {MIN_EXAMPLES}"
