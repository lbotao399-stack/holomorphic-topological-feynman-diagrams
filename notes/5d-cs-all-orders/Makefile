PYTHON ?= python3
LATEXMK ?= latexmk
ROOT := ../..
export OPENBLAS_NUM_THREADS := 1

.PHONY: all check check-new check-inherited clean

all: paper.pdf

paper.pdf: paper.tex $(wildcard sections/*.tex)
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -file-line-error -latexoption=-no-shell-escape paper.tex

check: check-new check-inherited

check-new:
	mkdir -p verification
	$(PYTHON) verify_lie_classification.py > lie_classification_checks.json
	$(PYTHON) checks/central_extension_exact.py > verification/central_extension_exact.json
	$(PYTHON) checks/cyclic_polynomial_referee.py > verification/cyclic_polynomial_referee.txt
	$(PYTHON) checks/paper_coefficients.py
	$(PYTHON) checks/four_loop_census.py

check-inherited:
	mkdir -p verification
	cd $(ROOT) && $(PYTHON) deepen/h1/verify/code/kernel_covariance.py > notes/5d-cs-all-orders/verification/kernel_covariance.txt
	cd $(ROOT) && $(PYTHON) notes/flagship/arith_check.py > notes/5d-cs-all-orders/verification/theta_constants.txt
	cd $(ROOT) && $(PYTHON) notes/flagship/referee2/r2_sl3_o3.py 7 notes/5d-cs-all-orders/verification/sl3_through_three_loops.txt
	cd $(ROOT) && $(PYTHON) deepen/o3/higher_orders/probe.py > notes/5d-cs-all-orders/verification/coefficient_detector.txt
	cd $(ROOT) && $(PYTHON) deepen/o3/higher_orders/four_loop_reduction/verify_four_point.py > notes/5d-cs-all-orders/verification/four_loop_full_cochain.txt

clean:
	$(LATEXMK) -c paper.tex
