UV := uv run --locked

.PHONY: check
check:
	$(UV) tools/checks/check-public.py
	$(UV) --group docs sphinx-build -W -q docs docs/_build/html
	$(UV) python -m unittest discover -s tests
