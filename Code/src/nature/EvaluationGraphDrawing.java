/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;
import java.util.List;
import java.util.Vector;

/**
 *
 * @author plantc59cs
 */
public class EvaluationGraphDrawing {

    Graph g;
    int[][] shortestPath;
    double[][] coord;
    int d;

    public EvaluationGraphDrawing(Graph g, double[][] coord) {
        this.g = g;
        //this.coord = coord;
        DataUtils du = new DataUtils();
        this.coord = du.scaleCoordinates(coord);
        this.d = coord[0].length;
    }

    public double stress() {
        DijkstraShortestPath<Integer, Integer> alg = new DijkstraShortestPath(g);
        shortestPath = new int[g.getVertexCount()][g.getVertexCount()];
        //double maxDist_mds = 0.0;
        for (int i = 0; i < g.getVertexCount(); i++) {
            for (int j = 0; j < g.getVertexCount(); j++) {
                if (i > j) {
                    List<Integer> l = alg.getPath(i, j);
                    if (l.size() > 0) {
                        shortestPath[i][j] = l.size();
                        shortestPath[j][i] = shortestPath[i][j];
//                        if (dist[i][j] > maxDist_mds) {
//                            maxDist_mds = dist[i][j];
//                        }
//                    } else {
//                        dist[i][j] = maxDist;
//                        dist[j][i] = maxDist;
//                    }
                    }
                }
            }
        }
        double sum = 0.0;
        for (int i = 0; i < g.getVertexCount(); i++) {
            for (int j = 0; j < i; j++) {
                double weight = 1.0 / shortestPath[i][j];
                double squaredDev = Math.pow(dist(i,j) - shortestPath[i][j], 2);
                sum += weight * squaredDev;
            }
        }
        return sum;
    }
    
    public double avgJaccardError(int numNeighbors){
        double avgErr = 0.0;
        for(int i = 0; i < g.getVertexCount(); i++){
            int[] pN = knnPaths(i, numNeighbors);
            int[] vN = knnEuclid(i, numNeighbors);
            Vector<Integer> union = new Vector<Integer>();
            Vector<Integer> intersection = new Vector<Integer>();
            for(int j = 0; j < pN.length; j++){
                if(!union.contains(pN[j]))
                    union.add(pN[j]);
                if(!union.contains(vN[j]))
                    union.add(vN[j]);
                boolean inIntersect = false;
                for(int k = 0; k < vN.length; k++)
                    if(pN[j] == vN[k])
                        inIntersect = true;
                if(inIntersect && !intersection.contains(pN[j]))
                    intersection.add(pN[j]);       
            }
            double jaccErr = 1.0 - (double)intersection.size()/(double)union.size();
            avgErr += jaccErr/(double)g.getVertexCount();
        }
        return avgErr;
    }
    
    private int[] knnPaths(int node, int k){
      double[] erg = new double[k];        //Array der Grösse MinPts anlegen
        int[] knn = new int[k];
        for (int i = 0; i < k; i++) {
            erg[i] = Double.MAX_VALUE;  //alle auf Max_Value setzen.
        }
        for (int i = 0; i < g.getVertexCount(); i++) {   //Durchlauf alle Datenobjekte
            if(node != i){
            int c = k - 1;
            double d = shortestPath[node][i]; //Distanz zum i-ten Datenobjekt
            while (c >= 0 && d < erg[c]) { //MinPts >= 0 && Distanz zum i-ten Datenobjekt < erg[MinPts-1]
                if (c < k - 1) {
                    erg[c + 1] = erg[c];
                    knn[c + 1] = knn[c];
                }
                erg[c] = d;
                knn[c] = i;
                c--;                  //MinPts um 1 reduzieren
            }
        }
        }
        return knn;
   
    }
    
    private int[] knnEuclid(int node, int k){
         double[] erg = new double[k];        //Array der Grösse MinPts anlegen
        //DataObject[] knn = new DataObject[k];  //für knn
        int[] knn = new int[k];
        for (int i = 0; i < k; i++) {
            erg[i] = Double.MAX_VALUE;  //alle auf Max_Value setzen.
        }
        for (int i = 0; i < g.getVertexCount(); i++) {   //Durchlauf alle Datenobjekte
if(node != i){
            int c = k - 1;
            double d = dist(node, i); //Distanz zum i-ten Datenobjekt
            while (c >= 0 && d < erg[c]) { //MinPts >= 0 && Distanz zum i-ten Datenobjekt < erg[MinPts-1]
                if (c < k - 1) {
                    erg[c + 1] = erg[c];
                    knn[c + 1] = knn[c];
                }
                erg[c] = d;
                knn[c] = i;
                c--;                  //MinPts um 1 reduzieren
            }
        }
        }
        return knn;

        
    }

    private double dist(int i, int j) {
        double res = 0.0;
        for (int l = 0; l < d; l++) {
            res += (coord[i][l] - coord[j][l]) * (coord[i][l] - coord[j][l]);
        }
        return Math.sqrt(res);

    }

}
