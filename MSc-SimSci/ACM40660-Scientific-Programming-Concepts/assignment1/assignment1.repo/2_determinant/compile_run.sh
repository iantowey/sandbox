echo "_______________________________________________"
echo "_______________________________________________"
echo "_______________________________________________"
echo "Matrix Calculations" 
echo "_______________________________________________"
echo "_______________________________________________"
echo "***********************************************"
echo "Make"
echo "***********************************************"
make clean && make
echo "***********************************************"
echo "Run application"
echo "***********************************************"
./det_app
echo "***********************************************"
echo "Check for memory leak"
echo "***********************************************"
valgrind --leak-check=full ./det_app 
