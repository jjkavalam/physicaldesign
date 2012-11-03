
final : ./bin/floorplan.out	

./bin/floorplan.out : ./lib/libanalfp.a ./lib/libtextutils.a ./lib/libparser.a ./src/floorplan.c
	gcc -g -static ./src/floorplan.c -L ./lib -lanalfp -lparser -ltextutils -lgsl -lgslcblas -lm -I ./include -o ./bin/floorplan.out
    
./lib/libanalfp.a : ./include/analfp.h ./src/analfp.c
	echo "building analytical floorplanning library"
	gcc -g -lgsl -lgslcblas -lm -I ./include -c ./src/analfp.c -o ./install/analfp.o
	ar rcs ./lib/libanalfp.a ./install/analfp.o

./lib/libtextutils.a : ./include/text_utils.h ./src/text_utils.c
	@echo "building text_utils library"
	gcc -g -c ./src/text_utils.c -o ./install/text_utils.o 
	ar  rcs ./lib/libtextutils.a ./install/text_utils.o

./lib/libparser.a : ./include/parser.h ./src/parser.c ./lib/libtextutils.a
	@echo "building parser library"
	gcc -g -c ./src/parser.c  -L ./lib -ltextutils -I ./include -o ./install/parser.o
	ar  rcs ./lib/libparser.a ./install/parser.o

clean : 
			 rm bin/* install/* 
