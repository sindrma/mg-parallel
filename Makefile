SHELL=/bin/sh
BENCHMARK=mg
PROGRAM=$(BENCHMARK)

include make.def
#TODO: need to include cuda library -L$CUDA_HOME/lib64 -lcudart
OBJS = mg.o functions/setup.o functions/results.o timer.o random.o utility.o

ifdef class
    CLASS=${class}
else
    CLASS=A
endif

CFLAGS += -g

mg: ${OBJS}

${OBJS} : npbparams.h

sys/npbparams.h:
	@echo "Generating npbparams.h.."
	make -C ./sys
	cd sys && ./setparams mg ${CLASS}

npbparams.h: sys/npbparams.h
	mv sys/npbparams.h ./

clean:
	- rm -f *.o *~
	- rm -f npbparams.h
	- rm -f mg
	- rm -f functions/*.o
	- if [ -d rii_files ]; then rm -r rii_files; fi

run: mg
	./mg
