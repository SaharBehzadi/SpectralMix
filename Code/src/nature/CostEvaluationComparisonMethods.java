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
public class CostEvaluationComparisonMethods {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        DataUtils du = new DataUtils();
        Graph g = ea.matlabToGraph("minnesotaL.mat", "graph");
        String[] path = new String[9];
        path[0] = "GEM.txt";
        path[1] = "multipole.txt";
        path[2] = "FM3.txt";
        path[3] = "GRIP.txt";
        path[4] = "LinLog.txt";
        path[5] = "SOM.txt";
        path[6] = "PivotMDS.txt";
        path[7] = "Stress.txt";
        path[8] = "minnesota.schulz";
        int numObj = 2640;
        int dim = 2;
        double[][] coord = new double[1][1];
        for (int i = 0; i < path.length; i++) {
            if (i != 5 && i != 8) {
                coord = du.txtToDouble(path[i], dim + 1, numObj);
                coord = du.removeLastCol(coord);
            } 
            if(i == 5){
                
                coord = du.readMatlabMatrix("coordSOM.mat", "coord");
                
            }
            if(i == 8){
                coord = du.txtToDouble(path[i], dim, numObj);
            }
            coord = du.scaleCoordinates(coord);
            GraphCompression gk = new GraphCompression(g, coord);
            System.out.println(path[i]);
            gk.mdlFunctionSimpleSigmoidComparisonMethods();
        }
        //Graph g = ea.matlabToGraph("twoMoons.mat", "g6");
//        double[][] coordSOM = ea.readMatlabMatrix("resultSOM.mat", "result");
//        GraphCompression gk = new GraphCompression(g, coordSOM);
//        System.out.println("costSOM: ");
//        gk.mdlFunctionSimpleSigmoidComparisonMethods();

//        //double[][] coordFR = ea.readMatlabMatrix("Airflights_coordinates.mat", "CoordMatrix2D");
//        double[][] coordFR = ea.readMatlabMatrix("twomoons_coordinates_2D.mat", "CoordMatrix2D");
//        GraphCompression gk = new GraphCompression(g, coordFR);
//        System.out.println("costFR: ");
//        gk.mdlFunctionSimpleSigmoidComparisonMethods();
//        double[][] coordOwn = ea.readMatlabMatrix("resultOwn.mat", "coord");
//        DataUtils du = new DataUtils();
//        double[][] coords = du.scaleLargestAxis(coordOwn);
//        gk = new GraphCompression(g, coords);
//        System.out.println("cost Own: ");
//        gk.mdlFunctionSimpleSigmoidComparisonMethods();
    }
}
