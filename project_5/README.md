This is a  repository for our program trying to simulate the motion of a wuantum wavepacket through a double slit experiment. 


The main program will initialize the wave packet and simulate its motion through a box, with or without a barrier and slits. It is run by :

 - g++ p5.cpp -larmadillo -fopenmp -o p5 
 - ./p5


The visualization tool used is based of the program "animate.py" it is run by in-built terminal in VSCode. It is necessary to have all files  produced by the main program iot. produce sufficient data for the plotting. The animation function is not made to create all plots, so it will need editing to get desired plots/animations.

The .gifs are the motion of the wavepacket given the different conditons on the system.
