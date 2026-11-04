UV := uv run --locked

.PHONY: check
check:
	$(UV) tools/checks/check-public.py
	$(UV) tools/docs/api_reference.py --check
	$(UV) python -m unittest discover -s tests
