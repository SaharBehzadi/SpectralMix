/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;

/**
 *
 * @author plantc59cs
 */
public class VisuMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
       //Graph g = ea.matlabToGraph("powerGrid.mat", "graph");
         //Graph g = ea.matlabToGraph("football.mat", "graph");
         //Graph g = ea.matlabToGraph("airflights.mat", "graph");
         
      
         
         //  Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
          // Graph g = ea.matlabToGraph("meshes.mat", "tapir");
           
//            Visualization v = new Visualization(g);
//                    double[][] coord = v.getCoordinatesISOM();
                    
        // double[][] coord = v.getCoordinatesKKL();
           
           
           
            //Graph g = ea.matlabToGraph("bus.mat", "graph");
            //Graph g = ea.matlabToGraph("dblp.mat", "graph");
            
             DataUtils du = new DataUtils();
//             du.saveAsMatlab(coord, "coord", "coordSOM.mat");
              //Graph g = du.readLuxemburg();
              //Graph g = du.readAdjList("dblpAdj");
              Graph g = du.readAdjList("pubmedAdj");
              //du.writeLineInputFile(g, "dblp");
         //double[][] coord = ea.readMatlabMatrix("result9.mat", "coord");
        // double[][] coord = ea.readMatlabMatrix("Airflight_coordinates.mat", "CoordMatrix2D");
        // double[][] coord = ea.readMatlabMatrix("resultWithBlowup.mat", "coord");
////        
int numObj = 19717;
int d = 2;
        //double[][] coord = du.txtToDouble("coord-before", d, numObj);
//         coord = du.removeLastCol(coord);
////         
         //du.saveAsMatlab(coord, "coord", "cNew.mat");
        double[][] coord = ea.readMatlabMatrix("coordLineTsne.mat", "coord");
         //double[][] coord = ea.readMatlabMatrix("initMinn100.mat", "init");
          //double[][] coord = ea.readMatlabMatrix("init.mat", "init");
        //double[][] labels = ea.readMatlabMatrix("football.mat", "labels");
         double[][] labels = ea.readMatlabMatrix("labelsDiabetes.mat", "labels");
         int[] id = new int[labels.length];
         for(int i = 0; i < id.length; i++)
             id[i] = (int)labels[i][0];
        Visualization v = new Visualization(g);
        v.displayCoordNew(coord, "bla", id);
        // v.displayCoordNew(coord, "bla");
         
//         EvaluationGraphDrawing ed = new EvaluationGraphDrawing(g, coord);
//         double stress = ed.stress();
//         double jE = ed.avgJaccardError(10);
//         System.out.println(stress + " " + jE);
//         
//         GraphCompression gk = new GraphCompression(g, coord);
//         gk.mdlFunctionSimpleSigmoidComparisonMethods();
         
    }
    
}
