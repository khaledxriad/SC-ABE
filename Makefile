.PHONY: all test results verify clean
PY ?= python3
all: test results verify
test:
	$(PY) -m pytest -q tests
results:
	$(PY) experiments/run_experiments.py
verify:
	$(PY) scripts/verify.py
clean:
	rm -rf .pytest_cache
