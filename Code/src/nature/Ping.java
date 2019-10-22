/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import java.util.Random;

/**
 *
 * @author claudia.plant
 */
public class Ping {

    double[][] adj;
    Random r = new Random(17);

    
    public Ping(double[][] adj) {
        this.adj  = adj;
    }

    public double[][] pingMatrix() {
        int n = adj.length;
        int numPings = 50;
        int walkLength = 1000;
        double[][] result = new double[numPings][];
        int[][] adjList = new int[n][];
        for (int i = 0; i < n; i++) {
            int count = 0;
            for (int j = 0; j < n; j++) {
                if (adj[i][j] != 0) {
                    count++;
                }
            }
            adjList[i] = new int[count];
            count = 0;
            for (int j = 0; j < n; j++) {
                if (adj[i][j] != 0) {
                    adjList[i][count++] = j;
                }
            }
        }
        for (int i = 0; i < numPings; i++) {
            result[i] = randomWalk(adjList, walkLength);
        }
        return result;
    }

    public double[] randomWalk(int[][] adjList, int walkLength) {
        int n = adjList.length;
        double[] result = new double[n];
        int k = r.nextInt(n);
        result[k] += 1.0 / walkLength;
        for (int i = 0; i < walkLength; i++) {
            k = adjList[k][r.nextInt(adjList[k].length)];
            result[k] += 1.0 / walkLength;
        }
        return result;
    }
}
