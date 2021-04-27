# ------------------------------------------------------------
# type "make" command in the Linux terminal create report.pdf
#
# type "make clean" to delete all generated files
# ------------------------------------------------------------

all: python
	
	make -C mat
	make -C sim
	make -C doc
	cp doc/report.pdf .
python:
	python2 t2_datagen.py 95804 > mat/data.txt

clean:
	make -C mat clean
	make -C sim clean
	make -C doc clean
	

.PHONY: all clean


