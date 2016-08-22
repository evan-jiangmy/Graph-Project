/**
 * @Description: This program is to implement a weighted graph to store and analyze crime records.
 * @author Evan Jiang
 * @Andrew_ID: mingyuaj
 * @Last_Modified_Date: 3/28/2016
 */
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.LinkedList;
import java.util.Scanner;
public class MainClass {
	public static void main(String[] args) {
//Read CSV file and construct vertex array and edge array********************************************************
		//User input to indicate start and end index
		int startIndex = 0;
		int endIndex = 0;
		
		//User input************************************************************
		try(Scanner Userinput = new Scanner(System.in)){
			System.out.println("Enter start index (the first index is 0): ");
			startIndex = Userinput.nextInt();
			
			System.out.println("Enter end index (the first index is 0): ");
			endIndex = Userinput.nextInt();
			
			//If the indexes are both 0, program terminates
		    if((startIndex == 0 && endIndex == 0) || (startIndex > endIndex)){
		    	System.out.println("No record is read or incorrect index combo. Program terminates.");
		    	System.exit(0);
		    }
		}
		//User input Ends***********************************************************
		
		
		//Import records from CSV File***************************************************
		int LineIndex = 0;
		//Create container to store user selected crime records.
    	LinkedList<String[]> CrimeRecordList = new LinkedList<String[]>();
    	
		try {
			Scanner input = new Scanner(new File("CrimeLatLonXY1990.csv"));
			input.useDelimiter(",");
			//Discard the first line (this line contains column names)
			input.nextLine();
		    
		    if(startIndex >= 0 && endIndex > 0){
		    	//Skip all lines until reaching the start index
		    	while(LineIndex < startIndex && input.hasNextLine()){
		    		input.nextLine();
		    		LineIndex++;
		    	}
		    	//Read every line and split them into string array
		    	while(LineIndex <= endIndex && input.hasNextLine()){
		    		String[] Capture = input.nextLine().split(",");
		    		CrimeRecordList.add(Capture);
		    		LineIndex++;
		    	}
		    }
		    input.close();
		}catch (FileNotFoundException e) {System.out.println(e.getMessage());}
		//Import records from CSV File Ends**********************************************
		
//Read CSV file and construct vertex array and edge array********************************************************
		
//Calculating vertex array and edge array************************************************************************
		int[] vertexArray = new int[CrimeRecordList.size()];
		double[][] edgeArray = new double [CrimeRecordList.size()][CrimeRecordList.size()];
		
		//Add Vertex
		int count = 0;
		for(int i = 0; i < vertexArray.length; i++){
			vertexArray[i] = count++;
		}
		
		//Add Edges and Distance
		try{
			//Calculate distance and add to edge array
			for(int i = 0; i < CrimeRecordList.size(); i++){
				double x = Double.valueOf(CrimeRecordList.get(i)[0]);
				double y = Double.valueOf(CrimeRecordList.get(i)[1]);
				for(int j = 0; j < CrimeRecordList.size(); j++){
					double x_Neighbor = Double.valueOf(CrimeRecordList.get(j)[0]);
					double y_Neighbor = Double.valueOf(CrimeRecordList.get(j)[1]);
					
					double Distance = Math.sqrt(Math.pow((x_Neighbor - x), 2) + Math.pow(y_Neighbor - y, 2)) * 0.00018939;
					edgeArray[i][j] = Distance;
				}
			}
		}catch(Exception e){System.out.println(e.getMessage());}
		
//Calculating vertex array and edge array Ends*******************************************************************

		
//Create Graph and execute Graph methods*************************************************************************
		Graph graph = new Graph(vertexArray, edgeArray);
		
		System.out.println("\nCrime Records Processed: ");
		for(int i = 0; i < CrimeRecordList.size(); i++){
			for(int j = 0; j < CrimeRecordList.get(i).length - 1; j++){
				System.out.print(CrimeRecordList.get(i)[j] + ", ");
			}
			System.out.print(CrimeRecordList.get(i)[CrimeRecordList.get(i).length - 1] + "\n");
		}
		System.out.println();
		
		//Calculate minimum distance using brute force method
		LinkedList<Integer> allElements = new LinkedList<Integer>();
		
		//Add all vertices to the list for permutation
		for (int i = 0; i < (endIndex - startIndex + 1); i++) allElements.add(i);
		
		//Empty LinkedList to facilitate permutation
	    LinkedList<Integer> afterMigration = new LinkedList<Integer>();
	    
	    //A list of list to store Combinations
	    LinkedList<LinkedList<Integer>> CombinationStorage = new LinkedList<LinkedList<Integer>>();
	    
	    //Execute permutation
	    graph.Permutation(afterMigration, allElements, CombinationStorage);
	    
	    //Calculate distances for all combinations
	    LinkedList<Double> DistanceList = new LinkedList<Double>();
	    graph.calculateMinDistance(CombinationStorage, DistanceList);
	    
	    /*
	    //*****TESTING*********************************
	    for(int i = 0; i < CombinationStorage.size(); i++){
	    	for(int j = 0; j < CombinationStorage.get(i).size(); j++){
	    		System.out.print(CombinationStorage.get(i).get(j) + " ");
	    	}
	    	System.out.println();
	    }
	    //TESINT***************************************
	    */
	    
	    //Find the smallest distance
	    LinkedList<Integer> MatchingIndexSet = new LinkedList<Integer>();
	    LinkedList<LinkedList<Integer>> minDistanceSets = new LinkedList<LinkedList<Integer>>();
	    double minDistance = DistanceList.get(0);
	    for(int i = 0; i < DistanceList.size(); i++){
	    	if(minDistance > DistanceList.get(i)){
	    		minDistance = DistanceList.get(i);
	    	}
	    }
	    
	    for(int i = 0; i < DistanceList.size(); i++){
	    	if(minDistance == DistanceList.get(i)){
	    		MatchingIndexSet.add(i);
	    		minDistanceSets.add(CombinationStorage.get(i));
	    	}
	    }
	    
	    //Remove duplicates
	    int RemoveThreshold = 0;
	    //Loop for set to be kept
	    for(int i = 0; i < minDistanceSets.size() - 1; i++){
	    	//Loop for the set to be tested and possibly deleted
	    	for(int k = i + 1; k < minDistanceSets.size(); k++){
	    		RemoveThreshold = 0;
	    		//Scan each element and compare 
		    	for(int j = 0; j < minDistanceSets.get(i).size() && j < minDistanceSets.get(k).size(); j++){
		    		//If head equals to tail, increase count
		    		if(minDistanceSets.get(i).get(j) == minDistanceSets.get(k).get(minDistanceSets.get(k).size() - 1 - j)){
		    			RemoveThreshold++;
		    		}
		    	}
		    	//If match count reaches the threshold, remove the duplicated list
		    	if(RemoveThreshold == minDistanceSets.get(k).size() && RemoveThreshold == minDistanceSets.get(i).size()){
		    		minDistanceSets.remove(k);
		    	}
	    	}
	    }
	    /*
	    //************TESTING****************************
	    for(int i = 0; i < minDistanceSets.size(); i++){
			System.out.println(minDistanceSets.get(i).toString() + " " + "Distance: " + DistanceList.get(MatchingIndexSet.get(i)));
		}
	    System.out.println();
	    //************TESTING****************************
	   */
	    
	    
	    
	    //TSP Path
	    LinkedList<Integer> TSPPathForKML = new LinkedList<Integer>();
	    graph.getMST(0, TSPPathForKML);
	    System.out.println();
	    
	    
	    
	    //Optimal Path
		System.out.print("Hamiltonan Cycle (minimal): ");
		
		if(minDistanceSets.size() > 1){
			for(int i = 0; i < minDistanceSets.get(1).size(); i++){
				System.out.print(minDistanceSets.get(1).get(i) + " ");
			}
		}else{
			for(int i = 0; i < minDistanceSets.get(0).size(); i++){
				System.out.print(minDistanceSets.get(0).get(i) + " ");
			}
		}

		System.out.println();
		
		System.out.print("Length of Cycle: ");
		System.out.print(minDistance + " miles.\n");
		
		//Processing KML
		if(minDistanceSets.size() > 1){
			toKML(CrimeRecordList, TSPPathForKML, minDistanceSets.get(1));
		}else{
			toKML(CrimeRecordList, TSPPathForKML, minDistanceSets.get(0));
		}
		
		
//Create Graph and execute Graph methods Ends********************************************************************
	}//End of main method
	/**
     * Organize location data in KML format and write the content to generate KML file
     * Pre-condition: The list contains locations of all related vertices, the TSP Path, and the Optimal Path are passed in.
     * Post-condition: KML file is generated with content
     * Best Case: Theata |V|
     * Worst Case: Theta |V|
	 * @param CrimeRecordList List contains crime record location data
	 * @param TSPPathForKML List contains TSP path vertices
	 * @param OptimalPathForKML List contains optimal path of vertices
	 */
	public static void toKML(LinkedList<String[]> CrimeRecordList, LinkedList<Integer> TSPPathForKML, LinkedList<Integer> OptimalPathForKML){
		String Prefix = "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n" + 
						"<kml xmlns=\"http://earth.google.com/kml/2.2\">\n" + 
						"<Document>\n" + 
						"<name>Pittsburgh TSP</name>\n" + 
						"<description>TSP on Crime</description>\n" + 
						"<Style id=\"style6\">\n" + 
						"<LineStyle>\n" + 
						"<color>73FF0000</color>\n" + 
						"<width>5</width>\n" + 
						"</LineStyle>\n" + 
						"</Style>\n" + 
						"<Style id=\"style5\">\n" + 
						"<LineStyle>\n" + 
						"<color>507800F0</color>\n" + 
						"<width>5</width>\n" + 
						"</LineStyle>\n" + 
						"</Style>\n" + 
						"<Placemark>\n" + 
						"<name>TSP Path</name>\n" + 
						"<description>TSP Path</description>\n" + 
						"<styleUrl>#style6</styleUrl>\n" + 
						"<LineString>\n" + 
						"<tessellate>1</tessellate>\n" + 
						"<coordinates>\n" + 
						"";
		
		//Build the TSP Path
		StringBuilder TSPPath = new StringBuilder();
		for(int i = 0; i < TSPPathForKML.size(); i++){
			TSPPath.append(CrimeRecordList.get(TSPPathForKML.get(i))[8] 
					+ "," + CrimeRecordList.get(TSPPathForKML.get(i))[7] 
					+ ",0.000000\n");
		}
		TSPPath.append("");
		System.out.println();
		System.out.println();
		
		
		String Middle = "</coordinates>\n" + 
						"</LineString>\n" + 
						"</Placemark>\n" + 
						"<Placemark>\n" + 
						"<name>Optimal Path</name>\n" + 
						"<description>Optimal Path</description>\n" + 
						"<styleUrl>#style5</styleUrl>\n" + 
						"<LineString>\n" + 
						"<tessellate>1</tessellate>\n" + 
						"<coordinates>\n" + 
						"";
		
		//Build the TSP Path
		StringBuilder OptimalPath = new StringBuilder();
		for(int i = 0; i < TSPPathForKML.size(); i++){
			OptimalPath.append(CrimeRecordList.get(OptimalPathForKML.get(i))[8] 
					+ "," + CrimeRecordList.get(OptimalPathForKML.get(i))[7] 
					+ ",0.000000\n");

		}
		
		String Suffix = "</coordinates>\n" + 
						"</LineString>\n" + 
						"</Placemark>\n" + 
						"</Document>\n" + 
						"</kml>\n";
		
		String KMLForWriter = Prefix + TSPPath.toString() + Middle + OptimalPath.toString() + Suffix;
		try {
			PrintWriter out = new PrintWriter("PGHCrimes.kml");
			out.print(KMLForWriter);
			out.close();
			System.out.println("\nThe crime data has been written to PGHCrimes.KML.");
		} 
		catch (FileNotFoundException e) {
			System.out.println(e.getMessage());
		}
		

	}
}//End of main class
