SOURCES = parkind1.F90 yomlun_ifsaux.F90 abor1.F90 yomhook.F90 \
    yomlun.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBIFSAUX = ../lib/libifsaux.a

all: $(LIBIFSAUX)

$(LIBIFSAUX): $(OBJECTS)
	ar r $(LIBIFSAUX) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(LIBIFSAUX)

yomlun_ifsaux.o, yomhook.o: parkind1.o
yomlun.o : yomlun_ifsaux.o
