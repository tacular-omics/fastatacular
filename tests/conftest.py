"""Hypothesis profiles: ``default`` keeps CI fast, ``thorough`` is for local soak runs.

Select with ``HYPOTHESIS_PROFILE=thorough uv run pytest tests``.
"""

import os

from hypothesis import HealthCheck, settings

settings.register_profile("default", max_examples=100, deadline=None)
settings.register_profile("thorough", max_examples=2000, deadline=None, suppress_health_check=[HealthCheck.too_slow])
settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "default"))
