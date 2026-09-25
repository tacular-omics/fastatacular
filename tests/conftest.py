"""Hypothesis profiles, selected with ``HYPOTHESIS_PROFILE``.

- ``default`` (30 examples): the fast local run, ``uv run pytest tests``.
- ``ci`` (100 examples, plus the larger per-test counts in ``test_fast_paths.py``):
  what CI runs, ``HYPOTHESIS_PROFILE=ci uv run pytest tests``.
- ``thorough`` (2000 examples): local soak runs, ``HYPOTHESIS_PROFILE=thorough uv run pytest tests``.
"""

import os

from hypothesis import HealthCheck, settings

settings.register_profile("default", max_examples=30, deadline=None)
settings.register_profile("ci", max_examples=100, deadline=None)
settings.register_profile("thorough", max_examples=2000, deadline=None, suppress_health_check=[HealthCheck.too_slow])
settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "default"))
