AGR=metae2e.agr ptmetae2e.agr

PDF=$(AGR:.agr=.pdf)

all: $(PNG) $(PDF)

.SUFFIXES : .agr .pdf

.agr.pdf:
	gracebat -printfile $*.pdf -hdevice PDF $*.agr


