/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;

/**
 *
 * @author claudia.plant
 */
public class DataUtilsMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        //Graph g = ea.matlabToGraph("minnesotaFull.mat", "minnesotaFull");
         //Graph g = ea.matlabToGraph("football.mat", "graph"); //random 70 unverfaltet; ist es durch Grid schlechter? 
         //Graph g = ea.matlabToGraph("meshes.mat", "eppstein"); //random 70 unverfaltet; ist es durch Grid schlechter? 
        // Graph g = ea.matlabToGraph("dti.mat", "graph"); //random 70 unverfaltet; ist es durch Grid schlechter? 
       //Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        //Graph g = ea.matlabToGraph("powerGrid.mat", "graph");
         //Graph g = ea.matlabToGraph("karate.mat", "graph");
         
         
         
        DataUtils du = new DataUtils();
       // Graph g = du.readAmazon();
       
        //Graph g = du.readAdjList("pubmedAdj");
//            ea.writeGraphToMatlab(g, "dblpLarge.mat");
        
        //du.writeAdjacencyList(g, "airflightsAdj");
        //du.graphToTulipAdj(g, "football.ascii");
       // du.graphToTulipAdj(g, "heli.ascii");
//        
        
        
        double[][] coord = du.txtToDouble();
       du.writeCSVFile(coord, "pubmedLineC.txt");
        //du.saveAsMatlab(coord, "coord", "dblp128d.mat");
        
        
        
         //Graph g = ea.matlabToGraph("dblp.mat", "graph"); //random 70 unverfaltet; ist es durch Grid schlechter? 
        
        // double[][] dd = ea.readMatlabMatrix("dtigraph.mat", "vertices");
        
      //double[][] dd = ea.readMatlabMatrix("1138_bus.mat", "graph");
       //Graph g = ea.matlabToGraph("fullGraph.mat", "g");
       // Graph g = ea.matlabToGraph("dblp_small.mat", "graph");
       
       //Graph g = ea.matlabToGraph("meshes.mat", "eppstein");
       //Graph g = ea.matlabToGraph("bus.mat", "graph");
        
      // double[][] coordGt = ea.readMatlabMatrix("fullGraph.mat", "gt");
//       double[][] coordLine = ea.readMatlabMatrix("twoMoonsExcel.mat", "Z_GraRep");
       
      //double[][] coord = ea.readMatlabMatrix( "coord2DafterTsne.mat", "minn2D");
      // double[][] labels = ea.readMatlabMatrix( "labelsDblp_small.mat", "labels");
          //Graph g = ea.matlabToGraph("epsilonGraphs.mat", "epsilonGraph_02");
            //Graph g = ea.matlabToGraph("netScienceGraph.mat", "graph");
            //Graph g = ea.matlabToGraph("epsilonGraphs.mat", "epsilonGraph_015");
//               Graph g = ea.matlabToGraph("polblogs.mat", "graph");
//               double[][] labels = ea.readMatlabMatrix("polblogs.mat", "labels");
             //    Graph g = ea.matlabToGraph("twoMoons.mat", "g2");
         //Graph g = ea.matlabToGraph("football.mat", "graph");
//         double[][] coord = ea.readMatlabMatrix("airflights.mat", "airflights");
//          double[][] labels = ea.readMatlabMatrix("labels.mat", "labels");
//          int[] l = new int[labels.length];
//          for(int i = 0; i < labels.length; i++)
//              l[i] = (int) labels[i][0];
//System.out.println(g.getVertexCount() + " " + g.getEdgeCount());
  //      DataUtils du = new DataUtils();
//        g = du.getLargestComponent(g);
//        System.out.println(g.getVertexCount() + " " + g.getEdgeCount());
        //IO io = new IO();
       // du.writeArff(coord, l, 7, "airflightsGraRep2D");
        
        //du.procrustesEvaluation(coordGt, coordLine); 
        //du.removeDegreeOneNodes(g, labels); 
       // du.order1Evaluation(g, coord);
//        du.writeLineInputFile(g, "pubmed");
//    Graph gg = du.getLargestComponent(g);
//    System.out.println("m");
    //  du.writeAdjacencyList(g, "dblp");
       //ea.writeGraphToMatlab(g, "bus1C");
       
      
    }
}
