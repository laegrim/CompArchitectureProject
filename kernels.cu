#ifndef _KERNELS_H_
#define _KERNELS_H_

struct Node {

}

__global__ void FWBW(Node * d_Graph, Node &pivot) {

	//If d_Graph

		//Get FW Set

		//Get BW Set

		//Find Intersection

		//Reorder Global Data

		//Call Recursively on Predecessors not in SCC
	
		//Call Recursively on Sucessors not in SCC

		//Call Recursively on Remainder

}

__global__ void FWBWTrim(Node * d_Graph, Node &pivot) {
	
	//if d_Graph

		//Trim Graph		
		
		//Get FW Set

		//Get BW Set

		//Find Intersection

		//Reorder Global Data

		//Call Recursively on Predecessors not in SCC

		//Call Recursively on Successors not in SCC

		//Call Recursively on Remainder

}

__global__ void Trim() {

	//Remove nodes with either no predecessor or no anscestor

	//Reorder Global Data

	//Call Recursively
}

//__global__ void ColorProp() {

//}
