
/**
 * @Description: This program is to implement a weighted graph to store and analyze crime records.
 * @author Evan Jiang
 * @Andrew_ID: mingyuaj
 * @Last_Modified_Date: 3/28/2016
 */

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class Graph {
	private List<Integer> vertexList;
	private List<LinkedList<Integer>> neighborList;
	private List<MinHeap> EdgeHeapList;
	private List<MinHeap> DuplicatedEdgeHeapList;
	//private List<LinkedList<Edge>> EdgeList;

	/**
	 * Constructor: Add the all involved vertices and edges into the graph
     * Pre-condition: A valid vertex array and a valid edges array are successfully passed in
     * Post-condition: A graph is constructed.
	 * @param vertexArray Array of the crime record
	 * @param edges Array of every location record and its adjacent location record
	 */
    protected Graph(int[] vertexArray, double[][] edges) {
     this.vertexList = new LinkedList<>();
     this.neighborList = new LinkedList<LinkedList<Integer>>();
     this.EdgeHeapList = new LinkedList<MinHeap>();
     this.DuplicatedEdgeHeapList = new LinkedList<MinHeap>();
     //this.EdgeList = new LinkedList<LinkedList<Edge>>();
      
      //Adding vertex
      for (int i = 0; i < vertexArray.length; i++){
    	  this.vertexList.add(vertexArray[i]);
    	  this.neighborList.add(new LinkedList<Integer>());
      }
      //Adding Edges
      for (int i = 0; i < edges.length; i++) {
    	  for (int j = 0; j < edges[i].length; j++){
    		  int TheVertex = i;
    		  int TheNeighborOfTheVertex = j;
    		  this.neighborList.get(TheVertex).add(TheNeighborOfTheVertex);
    	  }
      }
      //Adding distances
      for (int i = 0; i < this.vertexList.size(); i++) {
    	  this.EdgeHeapList.add(new MinHeap());
    	  this.DuplicatedEdgeHeapList.add(new MinHeap());
    	  //this.EdgeList.add(new LinkedList<Edge>());
      }
      for (int i = 0; i < edges.length; i++) {
    	  for (int j = 0; j < edges[i].length; j++){
    		  int TheVertex = i;
    		  int TheNeighborOfTheVertex = j;
    		  double Distance = edges[i][j]; 
    		  
    		  //Add Edge and distance to the MinHeap
    		  this.EdgeHeapList.get(TheVertex).add(new Edge(TheVertex, TheNeighborOfTheVertex, Distance));
    		  this.DuplicatedEdgeHeapList.get(TheVertex).add(new Edge(TheVertex, TheNeighborOfTheVertex, Distance));
    		 // this.EdgeList.get(TheVertex).add(new Edge(TheVertex, TheNeighborOfTheVertex, Distance));
    	  }
      }
      //Delete the edge with each vertex itself
      for(int i = 0; i < EdgeHeapList.size(); i++){
    	  this.EdgeHeapList.get(i).deleteMin();
    	  this.DuplicatedEdgeHeapList.get(i).deleteMin();
      }
      
      
    }
    
    /**
     * Print all the edges including distance information for all vertices
     * Pre-condition: A graph is constructed
     * Post-condition: The EdgeHeapList has been printed.
     * Best Case: Theta 1
     * Worst Case: Theta |V|^2
     */
    public void printDuplicatedEdgeHeapList(){
    	for(int i = 0; i < this.DuplicatedEdgeHeapList.size(); i++){
    		this.DuplicatedEdgeHeapList.get(i);
    		System.out.print("Vertex " + i + ": ");
    		for(int j = 0; j < this.DuplicatedEdgeHeapList.get(i).getHeapSize(); j++){
    			System.out.print(this.DuplicatedEdgeHeapList.get(i).getHeapList().get(j).getNeighbor() + ": "
    					+ this.DuplicatedEdgeHeapList.get(i).getHeapList().get(j).getDistance() + "  ");
    		}
    		System.out.println();
    	}
    }
    
    /**
     * Print all the edges including distance information for all vertices
     * Pre-condition: A graph is constructed
     * Post-condition: The EdgeHeapList has been printed.
     * Best Case: Theta 1
     * Worst Case: Theta |V|^2
     */
    /*
    public void printEdgeList(){
    	for(int i = 0; i < this.EdgeList.size(); i++){
    		this.EdgeList.get(i);
    		System.out.print("Vertex " + i + ": ");
    		for(int j = 0; j < this.EdgeList.get(i).size(); j++){
    			System.out.print(this.EdgeList.get(i).get(j).getNeighbor() + ": "
    					+ this.EdgeList.get(i).get(j).getDistance() + "  ");
    		}
    		System.out.println();
    	}
    }
    */
    
    
    /**
     * Print all the edges including distance information for all vertices
     * Pre-condition: A graph is constructed
     * Post-condition: The EdgeHeapList has been printed.
     * Best Case: Theta 1
     * Worst Case: Theta |V|*(|V|-1)
     */
    public void printEdgeHeapList(){
    	for(int i = 0; i < this.EdgeHeapList.size(); i++){
    		this.EdgeHeapList.get(i);
    		System.out.print("Vertex " + i + ": ");
    		for(int j = 0; j < this.EdgeHeapList.get(i).getHeapSize(); j++){
    			System.out.print(this.EdgeHeapList.get(i).getHeapList().get(j).getNeighbor() + ": "
    					+ this.EdgeHeapList.get(i).getHeapList().get(j).getDistance() + "  ");
    		}
    		System.out.println();
    	}
    }
    
    /**
     * Obtain the vertex based on index
     * Pre-condition: A graph is constructed
     * Post-condition: A vertex has been returned based on user input index.
     * Best Case: Theta 1
     * Worst Case: Theta n
     * @param VertexIndex User input index
     * @return The requested vertex to be returned
     */
    public int getVertex(int VertexIndex){
    	return vertexList.get(VertexIndex);
    }

    /**
     * Obtain the vertex list of the graph.
     * Pre-condition: A graph is constructed and it is not empty.
     * Post-condition: A vertex list has been returned.
     * Best Case: Theta 1
     * Worst Case: Theta 1
     * @return the vertex list of List Type
     */
    public List<Integer> getVertexList(){
    	return this.vertexList;
    }
    
    /**
     * Obtain the size of the vertex list
     * Pre-condition: A graph is constructed
     * Post-condition: A vertex list's size is returned.
     * Best Case: Theta 1
     * Worst Case: Theta 1
     * @return integer value of the vertex list's size.
     */
    public int getSize() {
      return this.vertexList.size();
    }
//*********************************************************Inner Edge Class************************************************************************    
 
    /**
     * Edge class for converting edge array to edge object and analyzing objects conveniently
     * @author Evan Jiang
     *
     */
    public static class Edge implements Comparable<Edge>{
      private int TheVertex; // Starting vertex of the edge
      private int TheNeighborOfTheVertex; // Ending vertex of the edge
      private double Distance;

      /**
       * Construct an Edge object
       * @param TheVertex One side of the edge
       * @param TheNeighborOfTheVertex The other side of the edge
       * @param Distance Distance of the edge
       */
      public Edge(int TheVertex, int TheNeighborOfTheVertex, double Distance) {
        this.TheVertex = TheVertex;
        this.TheNeighborOfTheVertex = TheNeighborOfTheVertex;
        this.Distance = Distance;
      }

      /**
       * Obtain the vertex based on index
       * Pre-condition: A graph is constructed
       * Post-condition: A vertex has been returned based on user input index.
       * Best Case: Theta 1
       * Worst Case: Theta 1
       * @return The requested vertex to be returned
       */
      public int getVertex(){
      	return this.TheVertex;
      }

      /**
       * Obtain the distance of the specific Edge
       * Pre-condition: A MinHeap is constructed.
       * Post-condition: The neighbor of the edge is returned
       * Best Case: Theta 1
       * Worst Case: Theta 1
       * @return double value of the neighbor
       */
      public int getNeighbor(){
    	  return this.TheNeighborOfTheVertex;
      }
      
      /**
       * Obtain the distance of the specific Edge
       * Pre-condition: A MinHeap is constructed.
       * Post-condition: The distance of the edge is returned
       * Best Case: Theta 1
       * Worst Case: Theta 1
       * @return double value of the distance
       */
      public double getDistance(){
    	  return this.Distance;
      }
      
      /**
       * Compare weights for Heap sort
       * Pre-condition: A MinHeap is constructed.
       * Post-condition: Edges are compared and an integer value is returned
       * Best Case: Theta 1
       * Worst Case: Theta 1
       * @param edge Another edge for comparison
       * @return integer variable to indicate comparison result
       */
      @Override
      public int compareTo(Edge edge) {
    	  if (this.Distance > edge.Distance)
    		  return 1;
    	  else if (this.Distance == edge.Distance)
    		  return 0;
    	  else
    		  return -1;
      }
    }
//************************************************Inner Edge Class Ends***************************************************************************

    
    
//************************************************Brute Force Methods****************************************************************
    /**
     * Find all the combinations of routes
     * Pre-condition: A graph is constructed, and all the vertices and edges are added.
     * Post-condition: all the combinations of routes have been found
     * Best Case: Theta |V|! 
     * Worst Case: Theta |V|!
     * @param afterMigration Linked List for adding new elements in the permutation
     * @param beforeMigration Linked List for removing new elements in the permutation
     * @param CombinationStorage Linked List for storing all combinations in the permutation
     */
    public void Permutation(LinkedList<Integer> afterMigration, LinkedList<Integer> beforeMigration, LinkedList<LinkedList<Integer>> CombinationStorage){
        
    	if(beforeMigration.size() > 0){
            for(int i = 0; i < beforeMigration.size(); i++){
            	//Remove the first element in the LinkedList (Queue)
                Integer firstOnList = beforeMigration.poll();
                
                //Declare a new list to store the copy of migrated elements
                LinkedList<Integer> afterMigrationDuplicated = new LinkedList<Integer>();
                afterMigrationDuplicated.addAll(afterMigration);
                //Add the migrated element to the list which is originally empty
                afterMigrationDuplicated.add(firstOnList);
                
                //Use recursion to generate combinations
                Permutation(afterMigrationDuplicated, beforeMigration, CombinationStorage);
                
                //Add back the first element when recursion returns to swap order of elements
                beforeMigration.add(firstOnList);
	        }
        }else{
        	//After generating one combination, store it to the list which stores all the combination
        	LinkedList<Integer> addToCombination = new LinkedList<Integer>();
        	
        	for(int i = 0; i < afterMigration.size(); i++){
        		addToCombination.add(afterMigration.get(i));
        	}
        	//Add the first one back to form a Hamiltonian Cycle
        	addToCombination.add(afterMigration.get(0));
        	
        	//System.out.println(afterMigration.toString());
        	
        	CombinationStorage.add(addToCombination);
        	
        	
        }
    }
    
    /**
     * Calculate distance for each combination.
     * Pre-condition: A list contains all combinations of routes, and a list to store calculated distance
     * Post-condition: Distance is calculated and stored in the list.
     * Best Case: Theta |V|^3
     * Worst Case: Theta |V|^3
     * @param CombinationStorage list contains all combinations of routes
     * @param DistanceList list to store calculated distance
     */
    public void calculateMinDistance(LinkedList<LinkedList<Integer>> CombinationStorage, LinkedList<Double> DistanceList){
    	for(int i = 0; i < CombinationStorage.size(); i++){
    		LinkedList<Integer> CombinationList = CombinationStorage.get(i);
    		
    		double SumOfDistance = 0;
    		
    		for(int j = 0; j < CombinationList.size() - 1; j++){
    			MinHeap EdgeSearch = this.DuplicatedEdgeHeapList.get(CombinationList.get(j));
    			for(int k = 0; k < EdgeSearch.getHeapSize(); k++){
    				if(EdgeSearch.getEdge(k).getNeighbor() == CombinationList.get(j + 1)){
    					SumOfDistance += EdgeSearch.getEdge(k).getDistance();
    				}
    			}
    		}
    		DistanceList.add(SumOfDistance);
    	}
    }

//************************************************Brute Force Methods End***********************************************************

    

//*****************************************************Prim's Algorithm***************************************************************************
    /**
     * Compare weights for Heap sort
     * Pre-condition: A MinHeap is constructed.
     * Post-condition: Edges are compared and an integer value is returned
     * Best Case: Theata |V|^3
     * Worst Case: Theta |V|^3
     * @param startingVertex Root of the minimum spanning tree
     * @param CycleStorageForKML Linked List to store TSP path vertices
     */
    public void getMST(int startingVertex, LinkedList<Integer> CycleStorageForKML){
    	//List to store all the visited vertices
        List<Integer> VisitedVertices = new LinkedList<Integer>();
        //Add the first vertex to the list
        VisitedVertices.add(startingVertex);

        //Check whether a vertex is visited
        boolean[] isVisited = new boolean[this.vertexList.size()];
        //The added vertex is set to be visited
        isVisited[startingVertex] = true;
        
        //Create an array for parents of vertices
        int vertexCount = this.vertexList.size();
        int[] parent = new int[vertexCount];
        //Initialize all the parents to -1
        for (int i = 0; i < parent.length; i++){
          parent[i] = -1;
        }
        
        //Start adding other vertices to the list and expanding the visited vertex list
        while (VisitedVertices.size() < vertexCount){
          //Initialize all neighbors to -1 and wait for the call
          int Neighbor = -1;
          double minDistance = Double.POSITIVE_INFINITY;
          
          for(int i = 0; i < VisitedVertices.size(); i++){
        	  //Delete all the visited neighbors for each vertex
        	  while(this.EdgeHeapList.get(VisitedVertices.get(i)).getHeapSize() > 0 && 
              		isVisited[this.EdgeHeapList.get(VisitedVertices.get(i)).PeekNearestNeighbor()]){
        		  this.EdgeHeapList.get(VisitedVertices.get(i)).deleteMin();
        	  }
        	  //Start finding the nearest neighbor of each vertex on the visited list
        	  if(this.EdgeHeapList.get(VisitedVertices.get(i)).getHeapSize() > 0){
            	  //Obtain the edge between nearest neighbor and the newly added element
                  Edge TheBestEdge = this.EdgeHeapList.get(VisitedVertices.get(i)).getBestEdge();
                
                  //Adjust distance value
                  if (TheBestEdge.getDistance() < minDistance) {
                	  Neighbor = TheBestEdge.getNeighbor();
                	  minDistance = TheBestEdge.getDistance();
                      //Add parent
                      parent[Neighbor] = VisitedVertices.get(i);
                  }
        	  }
          }//End of for loop

          //If the newly added element on the VisitedVertices has nearest neighbor, 
          //and at the same time its nearest neighbor is not visited, continue operation
          if (Neighbor != -1) {
        	  VisitedVertices.add(Neighbor); // Add a new vertex to the tree
        	  isVisited[Neighbor] = true;
        	  //SumOfDistance += minDistance;
          }
        } // End of while
 
        
        /*make an array to simulate the tree structure and store parent and children*/
        int[] MST = new int[(int) (Math.pow(2, parent.length - 1) - 1)];
        
        //Initialize all elements to be -2
        for(int i = 0; i < MST.length; i++){
        	MST[i] = -2;
        }
        
        //Call the tree construction method
        ConstructAndPrintMST(0, MST, parent);

        //Add the root back to make a Hamiltonian Cycle
        MST[MST.length - 1] = 0;
        
        //Print the Tree in Pre-Order format
        System.out.print("Hamiltonan Cycle (not necessarily optimum): ");
        LinkedList<Integer> MSTPreOrderList = new LinkedList<Integer>();
        MSTPreOrderTraversal(0, MST, MSTPreOrderList);
        for(int i = 0; i < MSTPreOrderList.size(); i++){
        	System.out.print(MSTPreOrderList.get(i) + " ");
        	//Store all the vertices in the MSTPreOrderList
        	CycleStorageForKML.add(MSTPreOrderList.get(i));
        }
        System.out.println();
        
        //Initialize the sum of distance of the best route
        double SumOfDistance = 0;
        
        //Find all the edges based on MST
        for(int i = 0; i < MSTPreOrderList.size() - 1; i++){
        	MinHeap EdgeSearch = this.DuplicatedEdgeHeapList.get(MSTPreOrderList.get(i));
        	for(int j = 0; j < EdgeSearch.getHeapSize(); j++){
        		if(EdgeSearch.getEdge(j).getNeighbor() == MSTPreOrderList.get(i + 1)){
        			SumOfDistance += EdgeSearch.getEdge(j).getDistance();
        		}
        	}
        }
        System.out.print("Length of Cycle: ");
        System.out.printf("%.2f miles\n", SumOfDistance);
    }

  //*****************************************************Prim's Algorithm Ends**************************************************************

    /**
     * Constructing the minimum spanning tree
     * Pre-condition: Prim's algorithm is executed and a tree like structure is constructed.
     * Post-condition: The minimum spanning tree has been returned.
     * Best Case: Theta 1
     * Worst Case: Theta |V|*Log|V|
     * @param MST integer array type of Minimum Spanning Tree which is waited to be constructed.
     * @param parentArray Parent Array to indicate the parent-children relationship
     * @return minimum spanning tree
     */
    private void ConstructAndPrintMST(int index, int[] MST, int[] parentArray){
    	
    	//Base case: return if the index reach to the end.
    	if(index > MST.length - 1) return;
    	
        //Inserting children
        for(int i = 0; i < parentArray.length; i++){
        	//Only the parent of the root is -1. Others are set to be -2
        	if(parentArray[i] == -1){
        		MST[0] = 0;
        	}else if(parentArray[i] == MST[index]){
        		if(((2 * index + 1) <= MST.length - 1) && (MST[2 * index + 1] == -2)) {
        			MST[2 * index + 1] = i;
        			
        		}else if (((2 * index + 2) <= MST.length - 1)
        				&& (MST[2 * index + 1] != -2)
        				&& (MST[2 * index + 2] == -2)){
        			MST[2 * index + 2] = i;
        		}else{
        			break;
        		}
        	}
        }
        //Recursion
        ConstructAndPrintMST(2 * index + 1, MST, parentArray);
        ConstructAndPrintMST(2 * index + 2, MST, parentArray);
    }
    
    /**
     * Print the Minimum Spanning Tree in Pre-Order format
     * Pre-condition: Prim's algorithm is executed and a tree like structure is constructed.
     * Post-condition: The minimum spanning tree has been printed in Pre-Order format.
     * Best Case: Theta |V|
     * Worst Case: Theta |V|
     * @param index The index to start the traversal
     * @param MST The Minimum Spanning Tree to be printed
     */
    private void MSTPreOrderTraversal(int index, int[] MST, LinkedList<Integer> MSTPreOrderList){
    	 if (index > MST.length - 1) return;
    	 if(MST[index] != -2) MSTPreOrderList.add(MST[index]);
    	 MSTPreOrderTraversal((2 * index + 1), MST, MSTPreOrderList); //left
    	 MSTPreOrderTraversal((2 * index + 2), MST, MSTPreOrderList); //right
    }
        
//******************************************************************MinHeap Class***********************************************************
    /**
     * MinHeap class for storing and sorting the Edges for a specific vertex
     * @author Evan
     *
     */
    public class MinHeap{
    	private ArrayList<Edge> HeapList;
    	
    	/**Constructor to create MinHeap*/
    	public MinHeap(){
    		this.HeapList = new ArrayList<Edge>();
    	}
    	
    	/**
    	 * Add Edge to the MinHeap, and swap will occur if necessary.
    	 * Pre-condition: A MinHeap is constructed.
    	 * Post-condition: The Edge object has been added.
         * Best Case: Theta 1
         * Worst Case: Theta Log|V|
    	 * @param edge The edge to be added.
    	 */
    	public void add(Edge edge) {
    		 //Add object
    		 this.HeapList.add(edge);
    		 // The index of the last node
    		 int currentIndex = this.HeapList.size() - 1;
    		 while (currentIndex > 0) {
    		 //Find parent's index
    			 int parentIndex = (currentIndex - 1) / 2;
    			 // Swap if the object is smaller than its parent
    			 if (this.HeapList.get(currentIndex).compareTo(this.HeapList.get(parentIndex)) == -1) {
    					 
    				 Edge NewParent = this.HeapList.get(currentIndex);
    				 this.HeapList.set(currentIndex, this.HeapList.get(parentIndex));
    				 this.HeapList.set(parentIndex, NewParent);
    					 
    			 }else break;
    			 
    			 currentIndex = parentIndex;
    			 }
    	 }
    	
    	
    	/**
    	 * Retrieve the specific edge based on index
    	 * Pre-condition: A MinHeap is constructed and it is not empty.
    	 * Post-condition: The Edge object has been returned.
         * Best Case: Theta 1
         * Worst Case: Theta 1
    	 * @param index specific index for the edge
    	 * @return The specific Edge
    	 */
    	public Edge getEdge(int index){
    		return this.HeapList.get(index);
    	}
    	
    	/**
    	 * Take a look at the nearest neighbor of the specific vertex
    	 * Pre-condition: A MinHeap is constructed.
    	 * Post-condition: Integer value of the nearest neighbor of the specific vertex has been returned.
         * Best Case: Theta 1
         * Worst Case: Theta 1
    	 * @return Integer value of the nearest neighbor of the specific vertex
    	 */
    	public int PeekNearestNeighbor(){
    		return this.HeapList.get(0).getNeighbor();
    	}

    	/**
    	 * Obtain the size of the MinHeap
    	 * Pre-condition: A MinHeap is constructed.
    	 * Post-condition: the size of the MinHeap has been returned
         * Best Case: Theta 1
         * Worst Case: Theta 1
    	 * @return the size of the MinHeap
    	 */
    	public int getHeapSize(){
    		return this.HeapList.size();
    	}
    	
    	/**
    	 * Obtain the Edge with the smallest distance in the MinHeap
    	 * Pre-condition: A MinHeap is constructed and it is not empty
    	 * Post-condition: The Edge with the smallest distance has been returned
         * Best Case: Theta 1
         * Worst Case: Theta 1
    	 * @return The Edge with the smallest value
    	 */
    	public Edge getBestEdge(){
    		return this.HeapList.get(0);
    	}
    	
    	/**
    	 * Delete the first element, which is the smallest element from the MinHeap. Swap occurs.
    	 * Pre-condition: A MinHeap is constructed.
    	 * Post-condition: The smallest element (the first element) is deleted, and corresponding swap occurs
         * Best Case: Theta 1
         * Worst Case: Theta logN
    	 * @return The smallest element in the MinHeap
    	 */
    	public Edge deleteMin() {
    		 if (this.HeapList.size() == 0) return null;
    		 
    		 //Obtain the first element in the heap
    		 Edge minEdge = this.HeapList.get(0);
    		 //Set the last element to be the temporal root
    		 this.HeapList.set(0, this.HeapList.get(this.HeapList.size() - 1));
    		 //Remove the last element since it is now the first one.
    		 this.HeapList.remove(this.HeapList.size() - 1);
 
    		 int currentIndex = 0;
    		 while (currentIndex < this.HeapList.size()) {
    			  int leftChildIndex = 2 * currentIndex + 1;
    			  int rightChildIndex = 2 * currentIndex + 2;
    			 
    			  //Find the minimum between two children
    			  //If all the elements are checked, break the loop
    			  if (leftChildIndex >= this.HeapList.size()) break;
    			  
    			  int minIndex = leftChildIndex;
    			  //If Left child is larger than the right child, set right child to be the new minimum value
    			  if (rightChildIndex < this.HeapList.size()) {
    				  if (this.HeapList.get(minIndex).compareTo(
    						  this.HeapList.get(rightChildIndex)) > 0) {
    					  			minIndex = rightChildIndex;
    				  }
    			  }
    			  //Swap if necessary
    			  if (this.HeapList.get(currentIndex).compareTo(this.HeapList.get(minIndex)) > 0) {
    			  Edge NewParent = this.HeapList.get(minIndex);
    			  this.HeapList.set(minIndex, this.HeapList.get(currentIndex));
    			  this.HeapList.set(currentIndex, NewParent);
    			  currentIndex = minIndex;
    			  }else break; //Heap sort complete
    			  } 
    		return minEdge;
    	}
    	
    	/**
    	 * Retrieve the list containing all egde pairs.
    	 * Pre-condition: A MinHeap is constructed.
    	 * Post-condition: The array list in the heap has been returned.
         * Best Case: Theta 1
         * Worst Case: Theta 1
    	 * @return The array list in the heap
    	 */
    	public ArrayList<Edge> getHeapList(){
    		return this.HeapList;
    	}

    }
//******************************************************************MinHeap Class Ends***********************************************************
	
	
	public static void main(String[] args) {

	}

}
