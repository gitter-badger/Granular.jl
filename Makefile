test: test-julia-0.6 test-julia-0.7

test-julia-0.6:
	julia --color=yes -e 'Pkg.test("Granular")'

test-julia-0.7:
	julia-0.7 --color=yes -e 'Pkg.test("Granular")'
