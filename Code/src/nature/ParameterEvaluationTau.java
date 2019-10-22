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
public class ParameterEvaluationTau {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        Graph g = ea.matlabToGraph("twoMoons.mat", "g6");
        double[] tau = new double[9];
        tau[0] = 1E-1;
        tau[1] = 0.7;
        tau[2] = 0.8;
        tau[3] = 0.9;
        tau[4] = 1.0;
        tau[5] = 0.6;
        tau[6] = 0.7;
        tau[7] = 0.8;
        tau[8] = 0.9;
        double[][] peRuntime = new double[9][2];

        for (int i = 0; i < 5; i++) {
            //for (int j = 0; j < 5; j++) {
                HierarchicalEmbedding he = new HierarchicalEmbedding(g, 0.8);
                long startTime = System.currentTimeMillis();
                he.run(tau[i]);
                long endTime = System.currentTimeMillis();
                peRuntime[i][0] += he.pe;
                peRuntime[i][1] += (double) (endTime - startTime) / 1000;
           // }
//            peRuntime[i][0] /= 5.0;
//            peRuntime[i][1] /= 5.0;
            System.out.println("Tau " + tau[i] + " pe: " + peRuntime[i][0] + " runtime: " + peRuntime[i][1]);
        }
        ea.writeDoubleToMatlab(peRuntime, "peRuntime.mat");
    }
}
