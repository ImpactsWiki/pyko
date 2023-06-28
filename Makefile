CC= gfortran
FLAGS1 =

all: 
	$(CC) $(FLAGS1) -o kov11e KOv11e.f 

clean:
	rm -f *.o kov11e
