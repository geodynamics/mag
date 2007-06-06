
all:
	cd src && make F77="$(F77)" all

clean:
	cd src && make F77="$(F77)" clean

