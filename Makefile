REPONAME=seaice

.PHONY: doc
doc:
	@make html -C doc/ && \
		cd .. && \
		git clone -b gh-pages git@github.com:anders-dc/$(REPONAME)
	@cd ../$(REPONAME)-doc/html && \
		pwd &&\
		git add . && \
		git commit -m 'rebuilt docs' && \
		git push origin gh-pages

.PHONY: coverage
coverage:
	@make coverage -C doc/ >/dev/null && \
		cat ../$(REPONAME)-doc/coverage/{python,c}.txt

.PHONY: doctest
doctest:
	@make coverage -C doc/ >/dev/null && \
		cat ../$(REPONAME)-doc/coverage/{python,c}.txt

manual.pdf:
	make latexpdf -C doc/
