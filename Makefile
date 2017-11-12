default: test

.PHONY: test
test: test-julia-0.6 test-julia-0.7

.PHONY: test-julia-0.6
test-julia-0.6:
	julia --color=yes -e 'Pkg.test("Granular")'

.PHONY: test-julia-0.6
test-julia-0.7:
	julia-0.7 --color=yes -e 'Pkg.test("Granular")'

.PHONY: docs
docs:
	cd docs && julia --color=yes make.jl
	open docs/build/index.html

