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
	
	//read the graph from the file, creating an adjacency matrix as representation
	//Graph file should be formatted as a multiline adjancency list, and the first
	//line should be the number of nodes in the graph
	int num_nodes;

	fp = fopen(argv[1], "r");
	if (!fp) {fprintf(stderr, "Error Reading Graph File"); exit(-1);}

	fscanf(fp, "%d", &num_nodes);
	
	int ** h_graph = (int **) malloc(sizeof(int *) * num_nodes);

	int start, degree;
	for (int i = 0; i < num_nodes; i++){
		fscanf(fp, "%d %d", &start, &degree);
		h_graph[start] = (int *) malloc(sizeof(int) * degree);
		for (int i = 0; i < degree; i++) {fscanf(fp, "%d", &h_graph[start][i]);}
	}
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
	double all_runtime = allend.tv_sec + allend.tv_usec / 1000000.0 - allstart.tv_sec - allstart.tv_usec / 1000000.0;
	double kruntime = kend.tv_sec + kend.tv_usec / 1000000.0 - kstart.tv_sec - kstart.tv_usec / 1000000.0;

	//Verify the Results
	

	//Free the memory


	return 0;

}
	
	
