.PHONY: release licenses
release: software-release/npr7150.docx

software-release/npr7150.docx: npr7150.md
	pandoc $< -o $@

licenses:
	conda run -n smap-download \
		pip-licenses \
		--format=csv \
		--with-urls \
		--output-file=software-release/licenses.csv
	conda run -n smap-download pip-licenses \
		--format=plain-vertical \
		--with-urls \
		--with-license-file \
		--output-file=software-release/licenses-full.txt
