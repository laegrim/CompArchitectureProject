/*
   Implementation of ColoringSCC

   */

#include <cstdlib>
#include <cstdio>
#include <cuda.h>
#include <sys/time.h>

int main(int argc, char*argv[]){

	printf("Coloring SCC v1.0\n");
	
	//check command line for arguments; we need a graph for input
	if (argc != 2) {fprintf(stderr, "usage: %s graph_file\n", argv[0]); exit(-1);}
	
	//Do some preprocessing and error checking here on the incoming graph
	
	
	//We want to time both the kernel and the all of the relevant processes
	timeval allstart, allend, kstart, kend;
	gettimeofday(&allstart, NULL);

	//Allocate the necessary memory on the device


	//Start the kernel timer
	gettimeofday(&kstart, NULL);

	//Call the kernel

	//end the timers
	gettimeofday(&kend, NULL);
	gettimeofday(&allend, NULL);

	//Print the Times


	//Verify the Results


	//Free the memory


	return 0;

}
	
	
