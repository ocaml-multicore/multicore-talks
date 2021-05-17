all: eio.pdf

%.pdf: %.tex
	pdflatex $<

docker:
	docker build -t latex .
	docker run --rm -v "${PWD}:/mnt" latex make -C /mnt
