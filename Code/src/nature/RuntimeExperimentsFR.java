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
public class RuntimeExperimentsFR {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        IO ea = new IO();
        Graph g50 = ea.matlabToGraph("g50.mat", "graph");
        Graph g100 = ea.matlabToGraph("g100.mat", "graph");
        Graph g150 = ea.matlabToGraph("g150.mat", "graph");
        Graph g200 = ea.matlabToGraph("g200.mat", "graph");
        Graph g300 = ea.matlabToGraph("g300.mat", "graph");
        Graph g400 = ea.matlabToGraph("g400.mat", "graph");
        Graph g500 = ea.matlabToGraph("g500.mat", "graph");
        Graph g600 = ea.matlabToGraph("g600.mat", "graph");
        Graph g700 = ea.matlabToGraph("g700.mat", "graph");
        Graph g800 = ea.matlabToGraph("g800.mat", "graph");
        Graph g900 = ea.matlabToGraph("g900.mat", "graph");
        Graph g1000 = ea.matlabToGraph("g1000.mat", "graph");

        double[][] rt = new double[12][1];

        long startTime = System.currentTimeMillis();
        Visualization v = new Visualization(g50);
        double[][] coord = v.getCoordinatesFRRuntime();
        long endTime = System.currentTimeMillis();
        double runtime = (endTime - startTime) / 1000;
        rt[0][0] = runtime;
        System.out.println("rt: " + runtime);

        startTime = System.currentTimeMillis();
        v = new Visualization(g100);
        coord = v.getCoordinatesFRRuntime();
        endTime = System.currentTimeMillis();
        runtime = (endTime - startTime) / 1000;
        rt[1][0] = runtime;
        System.out.println("rt: " + runtime);

        startTime = System.currentTimeMillis();
        v = new Visualization(g150);
        coord = v.getCoordinatesFRRuntime();
        endTime = System.currentTimeMillis();
        runtime = (endTime - startTime) / 1000;
        rt[2][0] = runtime;
        System.out.println("rt: " + runtime);

        startTime = System.currentTimeMillis();
        v = new Visualization(g200);
        coord = v.getCoordinatesFRRuntime();
        endTime = System.currentTimeMillis();
        runtime = (endTime - startTime) / 1000;
        rt[3][0] = runtime;
        System.out.println("rt: " + runtime);

        startTime = System.currentTimeMillis();
        v = new Visualization(g300);
        coord = v.getCoordinatesFRRuntime();
        endTime = System.currentTimeMillis();
        runtime = (endTime - startTime) / 1000;
        rt[4][0] = runtime;
        System.out.println("rt: " + runtime);

        startTime = System.currentTimeMillis();
        v = new Visualization(g400);
        coord = v.getCoordinatesFRRuntime();
        endTime = System.currentTimeMillis();
        runtime = (endTime - startTime) / 1000;
        rt[5][0] = runtime;
        System.out.println("rt: " + runtime);

        startTime = System.currentTimeMillis();
        v = new Visualization(g500);
        coord = v.getCoordinatesFRRuntime();
        endTime = System.currentTimeMillis();
        runtime = (endTime - startTime) / 1000;
        rt[6][0] = runtime;
        System.out.println("rt: " + runtime);

        startTime = System.currentTimeMillis();
        v = new Visualization(g600);
        coord = v.getCoordinatesFRRuntime();
        endTime = System.currentTimeMillis();
        runtime = (endTime - startTime) / 1000;
        rt[7][0] = runtime;
        System.out.println("rt: " + runtime);

        startTime = System.currentTimeMillis();
        v = new Visualization(g700);
        coord = v.getCoordinatesFRRuntime();
        endTime = System.currentTimeMillis();
        runtime = (endTime - startTime) / 1000;
        rt[8][0] = runtime;
        System.out.println("rt: " + runtime);

        startTime = System.currentTimeMillis();
        v = new Visualization(g800);
        coord = v.getCoordinatesFRRuntime();
        endTime = System.currentTimeMillis();
        runtime = (endTime - startTime) / 1000;
        rt[9][0] = runtime;
        System.out.println("rt: " + runtime);

        startTime = System.currentTimeMillis();
        v = new Visualization(g900);
        coord = v.getCoordinatesFRRuntime();
        endTime = System.currentTimeMillis();
        runtime = (endTime - startTime) / 1000;
        rt[10][0] = runtime;
        System.out.println("rt: " + runtime);

        startTime = System.currentTimeMillis();
        v = new Visualization(g1000);
        coord = v.getCoordinatesFRRuntime();
        endTime = System.currentTimeMillis();
        runtime = (endTime - startTime) / 1000;
        rt[11][0] = runtime;
        System.out.println("rt: " + runtime);

        DataUtils du = new DataUtils();
        du.saveAsMatlab(rt, "rt", "rt.mat");

    }

}
