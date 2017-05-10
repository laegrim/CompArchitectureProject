/*
   *REALLY naive Implementation of FWBW

   */

#include <cstdlib>
#include <cstdio>
#include <stack>
#include <cuda.h>
#include <sys/time.h>
#include <unistd.h>
#include <ctype.h>

typedef struct Node {
	unsigned int order;
	unsigned int out_degree;
	unsigned int in_degree;
	unsigned int * successors;
	unsigned int * predecessors;
	int scc;
	unsigned int fw_reachable;
	unsigned int bw_reachable;
	int subgraph;
} Node;

__constant__ unsigned int d_num_nodes;

__global__ void bfs(Node * graph, unsigned int * finished, unsigned int * mode, int * subgraph){

	if (threadIdx.x == 0) {atomicCAS(finished, 0, 1);}
	//syncthreads
	__threadfence();

	unsigned int old;

	unsigned int v = (blockIdx.x * blockDim.x) + threadIdx.x;

	Node curr_node, neighbor;
	
	for (unsigned int i = v; i < d_num_nodes; i++) {
		
		//get the current node
		curr_node = graph[i];

		//if the current node is in the subgraph
		if (curr_node.subgraph == *subgraph) {
			//if mode == 0, it's a forward search
			if (*mode == 0) {
				//for each of it's neighbors
				for(int i = 0; i < curr_node.out_degree; i++){
					//if the current node is forward reachable and the neighbor is in the subgraph
					neighbor = graph[curr_node.successors[i]];
					if (curr_node.fw_reachable && neighbor.subgraph == *subgraph) {
						//tell it's neighbors they are forward reachable
						old = atomicCAS(&(graph[curr_node.successors[i]].fw_reachable), 0, 1);
						//if something changed, we aren't at stasis
						if (!old) {atomicCAS(finished, 1, 0);}
					}
				}
			}
			//mode is 1, and it's a backward search on the transpose graph
			else {
				for(int i = 0; i < curr_node.in_degree; i++){
					if (curr_node.bw_reachable && neighbor.subgraph == *subgraph) {
						//tell the neighbors they are backwards reachable
						old = atomicCAS(&(graph[curr_node.predecessors[i]].bw_reachable), 0, 1);
						//if something changed, we aren't at stasis
						if (!old) {atomicCAS(finished, 1, 0);}
					}
				}
			}
		}
	}
}

__global__ void reset_bwfw_reachability(Node * graph){
	//reset node.fw_reachable & node.bw_reachable	
	 
	unsigned int v = (blockIdx.x * blockDim.x) + threadIdx.x;
	Node curr_node;
	
	for (unsigned int i = v; i < d_num_nodes; i++){	
		//get current node
		curr_node = graph[i];
		//reset it
		curr_node.fw_reachable = 0;
		curr_node.bw_reachable = 0;
	}
}

__global__ void trim_kernel(Node * graph, unsigned int * finished, int * subgraph) {

	if (threadIdx.x == 0) {atomicCAS(finished, 0, 1);}

	unsigned int v = (blockIdx.x * blockDim.x) + threadIdx.x;
	Node curr_node, neighbor;
	
	for (unsigned int i = v; i < d_num_nodes; i++){
		
		//get the current node
		curr_node = graph[i];
		
		//if the current node is in the subgraph
		if (curr_node.subgraph == *subgraph) {
			//visit each neighbor in the subgraph, telling them they are reachable
			for(int i = 0; i < curr_node.out_degree; i++){
				neighbor = graph[curr_node.successors[i]];
				if (neighbor.subgraph == *subgraph) {
					atomicCAS(&(graph[curr_node.successors[i]].fw_reachable), 0, 1);
				}
			}
		}
	}
		
	//sync threads, may need threadfence?
	__threadfence();

	
	for (unsigned int i = v; i < d_num_nodes; i++){
		//if you are not reachable or do not reach
		if (curr_node.fw_reachable == 0 || curr_node.out_degree == 0) {
			//then you are trimmable
			curr_node.scc = curr_node.order;
			curr_node.subgraph = -1;
			atomicCAS(finished, 1, 0);
		}
	}
} 

__global__ void assign_scc(Node * graph, int * scc) {
	//We want to find the smallest order of a node in the subgraph both reachable from and by the pivot
	//That will be our SCC label
	
	//This implementation is just an atomic min; a reduction could yeild much better throughput

	unsigned int v = (blockIdx.x * blockDim.x) + threadIdx.x;
	Node curr_node;
	
	for (unsigned int i = v; i < d_num_nodes; i++){
		
		//get the current node
		curr_node = graph[i];

		if (curr_node.fw_reachable && curr_node.bw_reachable) {atomicMin(scc, curr_node.order);}
		
	}
	
	//make sure all atomic writes are done
	__threadfence();
	
	//scc should be the lowest now
	curr_node.scc = *scc;
	curr_node.subgraph = -1;
	
}

__global__ void ancestor_partition(Node * graph, int * subgraph, unsigned int * empty){
	
	//We need to assign the anscestors who are not in the scc to a new subgraph
	unsigned int v = (blockIdx.x * blockDim.x) + threadIdx.x;
	Node curr_node;

	if (threadIdx.x == 0) {atomicCAS(empty, 0, 1);}
	
	__threadfence();
	
	for (unsigned int i = v; i < d_num_nodes; i++){
		
		//get the current node
		curr_node = graph[i];

		if (curr_node.bw_reachable && !curr_node.fw_reachable){curr_node.subgraph = *subgraph; atomicCAS(empty, 1, 0);}
	}
}

__global__ void descendent_partition(Node * graph, int * subgraph, unsigned int * empty){

	//We need to assign the descendents who are not in the scc to a new subgraph
	unsigned int v = (blockIdx.x * blockDim.x) + threadIdx.x;
	Node curr_node;

	if (threadIdx.x == 0) {atomicCAS(empty, 0, 1);}
	
	__threadfence();
	
	for (unsigned int i = v; i < d_num_nodes; i++){
		
		//get the current node
		curr_node = graph[i];

		if (!curr_node.bw_reachable && curr_node.fw_reachable){curr_node.subgraph = *subgraph; atomicCAS(empty, 1, 0);}
	}
}

__global__ void remainder_partition(Node * graph, int * subgraph, unsigned int * empty){
	//We need to assign the descendents who are not in the scc to a new subgraph
	unsigned int v = (blockIdx.x * blockDim.x) + threadIdx.x;
	Node curr_node;

	if (threadIdx.x == 0) {atomicCAS(empty, 0, 1);}
	
	__threadfence();
	
	for (unsigned int i = v; i < d_num_nodes; i++){
		
		//get the current node
		curr_node = graph[i];

		if (!curr_node.bw_reachable && !curr_node.fw_reachable){curr_node.subgraph = *subgraph; atomicCAS(empty, 1, 0);}
	}
}

__global__ void pivot(Node * graph, int * subgraph, unsigned int * pivot){
	//here we just assign the lowest order node in a subgraph to be the pivot
	//this is suboptimal

	//pivot should be initialized to num_nodes

	unsigned int v = (blockIdx.x * blockDim.x) + threadIdx.x;
	unsigned int votes;
	unsigned int lowest;
	unsigned int lane_id = threadIdx.x % 32;
	Node curr_node;

	for (unsigned int i = v; i < d_num_nodes; i++){
		//we don't want to hammer away with atomics, so we reduce the number of
		//threads using bandwidth

		//this can be further optimized		

		//get the current node
		curr_node = graph[i];
		//see how many nodes in the warp are in the subgraph
		votes = __ballot(curr_node.subgraph == *subgraph);
		//find the first node, and see if it's this node
		lowest = __ffs(votes) - 1;
		if (lane_id == lowest) atomicMin(pivot, curr_node.order);
	}

	__threadfence();
	//now that pivot should be finished, we need to initialize the corresponding node
	for (unsigned int i = v; i < d_num_nodes; i++){
		curr_node = graph[i];
		if (lane_id == lowest && *pivot == curr_node.order) {
			curr_node.fw_reachable = 1;
			curr_node.bw_reachable = 1;
		}
	}
}
		

int main(int argc, char ** argv){

	printf("Basic FWBW SCC v1.0\n");
	int option;
	bool output = false;
	char * out_file, * in_file;
	bool trim = false;
	unsigned int blocks = 0;
	unsigned int threads = 0;		
	
	while ((option = getopt(argc, argv, "o:tb:x:")) != -1) {
		switch (option) {
			case 'o':
				output = true;
				out_file = optarg;
				break;
			case 't':
				trim = true;
				break;
			case 'b':
				blocks = (unsigned int) atoi(optarg);
				break;
			case 'x':
				threads = (unsigned int) atoi(optarg);
				break;
			case '?':
				if (optopt == 'o' || optopt == 'b' || optopt == 'x') {fprintf(stderr, "Option -%c requires an argument\n", optopt); exit(-1);}
				else if (isprint(optopt)) {fprintf(stderr, "Unknown option -%c\n", optopt); exit(-1);}
				else {fprintf(stderr, "Unknown option \\x%x\n", optopt); exit(-1);}
			default:
				exit(-1);
		}
	}
	
	if (optind >= argc || optind < argc - 1) {fprintf(stderr, "Input Graph File Required\n"); exit(-1);}
	
	in_file = argv[optind];

	FILE * fp = fopen(in_file, "r");
	if (!fp) {fprintf(stderr, "Error Reading Graph File"); exit(-1);}
	
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
	
	Node * graph = NULL;
	
	//read the graph from the file
	//in this implementation, the graph representation is AOS
	int num_nodes, num_edges;
	int source, dest;
	
	//initialize based on header in graph file and get the degree of each node
	while((read=getline(&line, &len, fp)) != -1){
		//ignore header lines
		if (line[0] == 'p'){
			//this line has the total number of nodes and edges, we can make an empty graph
			sscanf(line, "p sp %d %d", &num_nodes, &num_edges);
			graph = (Node *) malloc(sizeof(Node) * num_nodes);
			for (int i = 0; i < num_nodes; i++) {
				graph[i].out_degree = 0;
				graph[i].successors = NULL;
				}
		}
		else if (line[0] == 'a'){
			//now we know how much memory to allocate for each node
			sscanf(line, "a %d %d %*d", &source, &dest);
			graph[source - 1].out_degree++;
			graph[dest - 1].in_degree++;
			}
	}	
	if (graph == NULL) {fprintf(stderr, "Error Reading Graph File, No Graph Allocated"); exit(-1);}
	rewind(fp);
	for (int i = 0; i < num_nodes; i++){
		//allocate the memory for each node
		graph[i].successors = (unsigned int *) malloc(sizeof(unsigned int) * graph[i].out_degree);
		graph[i].predecessors = (unsigned int *) malloc(sizeof(unsigned int) * graph[i].in_degree);
	}
	int s = 0;
	int c = 0;
	unsigned int * predec_tracker = (unsigned int *) calloc(num_nodes, sizeof(unsigned int));
	//record edges in the graph structure
	//the format is sorted by source node
	while((read=getline(&line, &len, fp)) != -1){
		if (line[0] == 'a'){
			sscanf(line, "a %d %d %*d", &source, &dest);
			//record the transverse edge
			graph[dest - 1].predecessors[predec_tracker[dest - 1]] = (source - 1);
			predec_tracker[dest - 1]++;
			//record the edge, done this way we don't need another large array
			if (source - 1 == s){
				graph[s].successors[c] = (dest - 1);
				c++;
			}
			else {
				s = source - 1;
				c = 0;
				graph[s].successors[c] = (dest - 1);
				c++;
			}
		}
	}

	fclose(fp);
	free(line);
	free(predec_tracker);

	//Do some preprocessing and error checking here on the incoming graph
	for (unsigned int i = 0; i < num_nodes; i++){
		graph[i].order = i;
		graph[i].scc = -1;
		graph[i].subgraph = 0;
	}
			
	//if blocks and threads aren't pre-assigned, assign them	
	if (threads == 0) {threads = 512;}
	if (blocks == 0) {blocks = (num_nodes + (threads - 1)) / threads;}

	//set up the stack
	std::stack<unsigned int> subgraphs;
	subgraphs.push(0);
		
	//We want to time both the kernel and the all of the relevant processes
	timeval allstart, allend, kstart, kend;
	gettimeofday(&allstart, NULL);

	//Allocate the necessary memory on the device, and copy over any initial data
	Node * d_graph;
	unsigned int h_finished, h_mode, h_subgraph, h_empty;
	unsigned int * d_finished, * d_mode, * d_empty, * d_pivot;
	int * d_subgraph, * d_scc;
	
	if (cudaSuccess != cudaMalloc((void **) &d_graph, sizeof(Node) * num_nodes)) {fprintf(stderr, "Couldn't allocate d_graph\n"); exit(-1);}
	if (cudaSuccess != cudaMalloc((void **) &d_finished, sizeof(unsigned int))) {fprintf(stderr, "Couldn't allocate d_finished\n"); exit(-1);}
	if (cudaSuccess != cudaMalloc((void **) &d_mode, sizeof(unsigned int))) {fprintf(stderr, "Couldn't allocate d_mode\n"); exit(-1);}
	if (cudaSuccess != cudaMalloc((void **) &d_subgraph, sizeof(int))) {fprintf(stderr, "Couldn't allocate d_subgraph\n"); exit(-1);}
	if (cudaSuccess != cudaMalloc((void **) &d_scc, sizeof(int))) {fprintf(stderr, "Couldn't allocate d_scc\n"); exit(-1);}
	if (cudaSuccess != cudaMalloc((void **) &d_empty, sizeof(unsigned int))) {fprintf(stderr, "Couldn't allocate d_empty\n"); exit(-1);}
	if (cudaSuccess != cudaMalloc((void **) &d_pivot, sizeof(unsigned int))) {fprintf(stderr, "Couldn't allocate d_pivot\n"); exit(-1);}

	cudaMemcpyToSymbol(d_num_nodes, &num_nodes, sizeof(unsigned int));
	
	//we need to copy over the graph struct, and allocate/copy the predecessor and successor arrays
	if (cudaSuccess != cudaMemcpy(d_graph, graph, sizeof(Node) * num_nodes, cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy graph to device\n"); exit(-1);}
	unsigned int * d_successors, * d_predecessors;
	for (int i = 0; i < num_nodes; i++){
		if (cudaSuccess != cudaMalloc((void **) &d_successors, graph[i].out_degree * sizeof(unsigned int))) {fprintf(stderr, "Couldn't allocate a successors array for node %d\n", i); exit(-1);}
		if (cudaSuccess != cudaMalloc((void **) &d_predecessors, graph[i].in_degree * sizeof(unsigned int))) {fprintf(stderr, "Couldn't allocate a predecessors array for node %d\n", i); exit(-1);}
		if (cudaSuccess != cudaMemcpy(d_successors, graph[i].successors, graph[i].out_degree * sizeof(unsigned int), cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy successors array for node %d\n", i); exit(-1);}
		if (cudaSuccess != cudaMemcpy(d_predecessors, graph[i].predecessors, graph[i].in_degree * sizeof(unsigned int), cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy predeccessors array for node %d\n", i); exit(-1);}
		if (cudaSuccess != cudaMemcpy(&(d_graph->successors), &d_successors, sizeof(unsigned int *), cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy successors pointer for node %d\n", i); exit(-1);}
		if (cudaSuccess != cudaMemcpy(&(d_graph->predecessors), &d_predecessors, sizeof(unsigned int *), cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy predecessors pointer for node %d\n", i); exit(-1);}
	}

	//Start the kernel timer
	gettimeofday(&kstart, NULL);

	//while there's still work to do (subgraphs on the stack)
	while (!subgraphs.empty()){
	
		//take a subgraph off the stack	
		h_subgraph = subgraphs.top();
		subgraphs.pop();
		if (cudaSuccess != cudaMemcpy(d_subgraph, &h_subgraph, sizeof(unsigned int), cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy d_subgraph to device\n"); exit(-1);}

		//reset it's forward and backward reachable attributes
		reset_bwfw_reachability<<<blocks, threads>>>(d_graph);

		//if trim, then trim till you can trim no more
		if (trim){
			h_finished = 0;
			if (cudaSuccess != cudaMemcpy(d_finished, &h_finished, sizeof(unsigned int), cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy d_finished to device\n"); exit(-1);}
			while (!h_finished){
				//trim
				trim_kernel<<<blocks, threads>>>(d_graph, d_finished, d_subgraph);
				//reset the reachability
				reset_bwfw_reachability<<<blocks, threads>>>(d_graph);
				//check if finished
				if (cudaSuccess != cudaMemcpy(&h_finished, d_finished, sizeof(unsigned int), cudaMemcpyDeviceToHost)) {fprintf(stderr, "Couldn't copy d_finished to host\n"); exit(-1);}
			}
		}
						
		//choose a pivot
		pivot<<<blocks, threads>>>(d_graph, d_subgraph, d_pivot);

		//do a forwards reachability search
		h_finished = 0;
		h_mode = 0;
		if (cudaSuccess != cudaMemcpy(d_finished, &h_finished, sizeof(unsigned int), cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy d_finished to device\n"); exit(-1);}
		if (cudaSuccess != cudaMemcpy(d_mode, &h_mode, sizeof(unsigned int), cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy d_mode to device\n"); exit(-1);}
		while (!h_finished){
			bfs<<<blocks, threads>>>(d_graph, d_finished, d_mode, d_subgraph);
			if (cudaSuccess != cudaMemcpy(&h_finished, d_finished, sizeof(unsigned int), cudaMemcpyDeviceToHost)) {fprintf(stderr, "Couldn't copy d_finished to host\n"); exit(-1);}
		}

		//do a backwards reachability search
		h_finished = 0;
		h_mode = 1;
		if (cudaSuccess != cudaMemcpy(d_finished, &h_finished, sizeof(unsigned int), cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy d_finished to device\n"); exit(-1);}
		if (cudaSuccess != cudaMemcpy(d_mode, &h_mode, sizeof(unsigned int), cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy d_mode to device\n"); exit(-1);}
		while (!h_finished){
			bfs<<<blocks, threads>>>(d_graph, d_finished, d_mode, d_subgraph);
			if (cudaSuccess != cudaMemcpy(&h_finished, d_finished, sizeof(unsigned int), cudaMemcpyDeviceToHost)) {fprintf(stderr, "Couldn't copy d_finished to host\n"); exit(-1);}
		}

		//Add the nodes in the union of the FW and BW reachable sets to an SCC
		assign_scc<<<blocks, threads>>>(d_graph, d_scc);

		//Add the FW reachable nodes not in the SCC to the subgraph list
		h_subgraph++;
		if (cudaSuccess != cudaMemcpy(d_subgraph, &h_subgraph, sizeof(unsigned int), cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy d_finished to device\n"); exit(-1);}
		descendent_partition<<<blocks, threads>>>(d_graph, d_subgraph, d_empty);
		if (cudaSuccess != cudaMemcpy(&h_empty, d_empty, sizeof(unsigned int), cudaMemcpyDeviceToHost)) {fprintf(stderr, "Couldn't copy d_finished to host\n"); exit(-1);}
		if (!h_empty) {subgraphs.push(h_subgraph);}

		//Add the BW reachable nodes not in the SCC to the subgraph list
		h_subgraph++;
		if (cudaSuccess != cudaMemcpy(d_subgraph, &h_subgraph, sizeof(unsigned int), cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy d_finished to device\n"); exit(-1);}
		ancestor_partition<<<blocks, threads>>>(d_graph, d_subgraph, d_empty);
		if (cudaSuccess != cudaMemcpy(&h_empty, d_empty, sizeof(unsigned int), cudaMemcpyDeviceToHost)) {fprintf(stderr, "Couldn't copy d_finished to host\n"); exit(-1);}
		if (!h_empty) {subgraphs.push(h_subgraph);}

		//Add the nodes not reachable from the pivot to to the subgraph list
		h_subgraph++;
		if (cudaSuccess != cudaMemcpy(d_subgraph, &h_subgraph, sizeof(unsigned int), cudaMemcpyHostToDevice)) {fprintf(stderr, "Couldn't copy d_finished to device\n"); exit(-1);}
		remainder_partition<<<blocks, threads>>>(d_graph, d_subgraph, d_empty);
		if (cudaSuccess != cudaMemcpy(&h_empty, d_empty, sizeof(unsigned int), cudaMemcpyDeviceToHost)) {fprintf(stderr, "Couldn't copy d_finished to host\n"); exit(-1);}
		if (!h_empty) {subgraphs.push(h_subgraph);}
	}
	//end the timers
	gettimeofday(&kend, NULL);
	gettimeofday(&allend, NULL);

	//Print the Times
	double all_runtime = allend.tv_sec + allend.tv_usec / 1000000.0 - allstart.tv_sec - allstart.tv_usec / 1000000.0;
	double kruntime = kend.tv_sec + kend.tv_usec / 1000000.0 - kstart.tv_sec - kstart.tv_usec / 1000000.0;

	for (int i = 0; i < num_nodes; i++) {
		free(graph[i].successors);
		free(graph[i].predecessors);

	//Print the Results, if necessary
	if (cudaSuccess != cudaMemcpy(graph, d_graph, sizeof(Node) * num_nodes, cudaMemcpyDeviceToHost)) {fprintf(stderr, "Couldn't copy graph to host\n"); exit(-1);}
	
	if (output == true) {
		fp = fopen(out_file, "w+");
		int scc = 0;	
		int count = 0;
		bool first;
		while (count < num_nodes) {
			first = true;
			for (int i = 0; i < num_nodes; i++) {
				if (graph[i].scc == scc) {
					count++;
					if (first) {
						fprintf(fp, "\nSCC %d: ", scc);
						first = false;
					}
					fprintf(fp, "%d ", graph[i].order);
				}

			}
			scc++;
		}
		fprintf(fp, "\n");
		fclose(fp);
	}

	//Free the memory
	free(graph);
	cudaFree(d_finished);
	cudaFree(d_mode);
	cudaFree(d_scc);
	cudaFree(d_pivot);
	cudaFree(d_empty);
	cudaFree(d_subgraph);
	for (int i = 0; i < num_nodes; i++) {
		cudaFree(d_graph->successors);
		cudaFree(d_graph->predecessors);
	}
	cudaFree(d_graph);

	return 0;

}
	
	
