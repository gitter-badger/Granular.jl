REPONAME=seaice

default: docs

.PHONY: docs
docs:
	-mkdir -p ../$(REPONAME)-doc && \
		cd ../$(REPONAME)-doc && \
		git clone -b gh-pages git@github.com:anders-dc/$(REPONAME) html
	-make html -C doc/ && cd ../$(REPONAME)-doc/html && \
		pwd && \
		git add . && \
		git commit -m 'rebuilt docs' && \
		git push origin gh-pages
	-make latexpdf -C doc/ && \
		cp ../$(REPONAME)-doc/latex/seaice.pdf . && \
		git add seaice.pdf && \
		git commit -m 'rebuild docs' && \
		git push

.PHONY: coverage
coverage:
	@make coverage -C doc/ >/dev/null && \
		cat ../$(REPONAME)-doc/coverage/{python,c}.txt

.PHONY: doctest
doctest:
	@make coverage -C doc/ >/dev/null && \
		cat ../$(REPONAME)-doc/coverage/{python,c}.txt

clean:
	make clean -C doc/
