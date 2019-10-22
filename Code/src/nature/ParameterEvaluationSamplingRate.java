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
public class ParameterEvaluationSamplingRate {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        //Graph g = ea.matlabToGraph("epsilonGraphs.mat", "epsilonGraph_015");
           Graph g = ea.matlabToGraph("twoMoons.mat", "g6");
        double[] samplingRate = new double[9];
        samplingRate[0] = 0.1;
        samplingRate[1] = 0.2;
        samplingRate[2] = 0.3;
        samplingRate[3] = 0.4;
        samplingRate[4] = 0.5;
        samplingRate[5] = 0.6;
        samplingRate[6] = 0.7;
        samplingRate[7] = 0.8;
        samplingRate[8] = 0.9;
        double[][] peRuntime = new double[9][2];

       // for (int i = 5; i < samplingRate.length; i++) {
            for (int j = 0; j < 10; j++) {
                //HierarchicalEmbedding he = new HierarchicalEmbedding(g, samplingRate[i]);
                HierarchicalEmbeddingWeighted he = new HierarchicalEmbeddingWeighted(g, samplingRate[j]);
                long startTime = System.currentTimeMillis();
                he.run();
                long endTime = System.currentTimeMillis();
                peRuntime[j][0] += he.pe;
                peRuntime[j][1] += (double) (endTime - startTime) / 1000;
                 System.out.println("Sampling rate: " + samplingRate[j] + " pe: " + peRuntime[j][0] + " runtime: " + peRuntime[j ][1]);
            }
//            peRuntime[i][0] /= 10.0;
//            peRuntime[i][1] /= 10.0;
           // System.out.println("Sampling rate: " + samplingRate[i] + " pe: " + peRuntime[i][0] + " runtime: " + peRuntime[i][1]);
       // }
        ea.writeDoubleToMatlab(peRuntime, "peRuntime.mat");
    }
}
