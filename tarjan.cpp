/*
   Implementation of Tarjan's SCC Algorithm

   */

#include <cstdlib>
#include <cstdio>
#include <sys/time.h>
#include <stack>
#include <algorithm>

struct Node {
	int order;
	int out_degree;
	int * successors;
	int lowlink;
	int index;
	bool onstack;
	int scc;
};

void strongconnect(Node * graph, std::stack<Node*> &S, Node &node, int &index) {
	//from wikipedia pseudocode
	
	Node * comp;
	node.index = index;
	node.lowlink = index;
	node.onstack = true;
	index += 1;
	S.push(&node);
	
	for (int i = 0; i < node.out_degree; i++){
		if (graph[node.successors[i]].index == -1) {
			strongconnect(graph, S, graph[node.successors[i]], index);
			node.lowlink = std::min(node.lowlink, graph[node.successors[i]].lowlink);
		}
		else if (graph[node.successors[i]].onstack) {
			node.lowlink = std::min(node.lowlink, graph[node.successors[i]].index);
		}
	}
	
	if (node.lowlink == node.index) {
		do {
			comp = S.top();
			S.pop();
			comp->scc = node.order;
			comp->onstack = false;
		}
		while (comp != &node);
	}			 
}

int main(int argc, char*argv[]){

	printf("Tarjan SCC v1.0\n");
	
	//check command line for arguments; we need a graph for input
	if (argc != 2 and argc != 3) {fprintf(stderr, "usage: %s <graph_file> <output_file>(optional)\n", argv[0]); exit(-1);}
	
	//read the graph from the file, creating an adjacency matrix as representation
	//Graph file should be formatted as a multiline adjancency list, and the first
	//line should be the number of nodes in the graph
	int num_nodes;
	int num_edges;

	FILE * fp = fopen(argv[1], "r");
	if (!fp) {fprintf(stderr, "Error Reading Graph File"); exit(-1);}
	
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
	
		
	Node * graph;
	//get the number of degrees on each edge
	int source, dest;
	
	//initialize based on header in graph file and get the degree of each node
	while((read=getline(&line, &len, fp))){
		//ignore header lines
		if (line[0] == 'p'){
			//this line has the total number of nodes and edges, we can make an empty graph
			sscanf(line, "p sp %d %d", &num_nodes, &num_edges);
			graph = (Node *) malloc(sizeof(Node) * num_nodes);
			for (int i = 0; i < num_nodes; i++) {graph[i].out_degree = 0;}
		}
		else if (line[0] == 'a'){
			//now we know how much memory to allocate for each node
			sscanf(line, "a %d %*d %*d", &source);
			graph[source].out_degree++;
			}
	}		
	rewind(fp);
	for (int i = 0; i < num_nodes; i++){
		//allocate the memory for each node
		graph[i].successors = (int *) malloc(sizeof(int) * graph[i].out_degree);
	}
	int c, s = 0;
	//record edges in the graph structure
	while((read=getline(&line, &len, fp))){
		if (line[0] == 'a'){
			sscanf(line, "a %d %d %*d", &source, &dest);
			if (source == s){
				graph[source].successors[c] = dest;
				c ++;
			}
			else {
				s = source;
				c = 0;
				graph[source].successors[c] = dest;
				c++;
			}
		}
	}

	fclose(fp);

	for (int i = 0; i < num_nodes; i++){
		graph[i].onstack = false;
		graph[i].index = -1;
		graph[i].order = i;
	}

	//We want to time the relevant processes
	timeval start, end;
	gettimeofday(&start, NULL);

	int index = 0;
	std::stack<Node *> S;
	for (int i = 0; i < num_nodes; i++){
		if (graph[i].index == -1) {strongconnect(graph, S, graph[i], index);}
	}
		
	//end the timers
	gettimeofday(&end, NULL);

	//Print the Time and the SCC's
	double runtime = end.tv_sec + end.tv_usec / 1000000.0 - start.tv_sec - start.tv_usec / 1000000.0;
	fprintf(stdout, "Time: %f\n", runtime);
	
	if (argc == 3) {
		fp = fopen(argv[2], "w+");
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
	for (int i = 0; i < num_nodes; i++) {free(graph[i].successors);}
	free(graph);

	return 0;

}
	
