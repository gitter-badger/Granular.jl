.PHONY: doc

doc:
	@make html -C doc/ && \
		pwd &&\
		cd ../seaice-doc/html && \
		pwd &&\
		git add . && \
		git commit -m 'rebuilt docs' && \
		git push origin gh-pages

manual.pdf:
	make latexpdf -C doc/
