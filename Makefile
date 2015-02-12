match : main.cpp lap.cpp lap.h gnrl.h
	g++ main.cpp lap.cpp -lopencv_core -lopencv_imgproc -lopencv_highgui -o match

clean :
	rm match
