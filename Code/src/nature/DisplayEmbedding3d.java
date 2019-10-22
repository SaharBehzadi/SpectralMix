/*
 * Copyright (c) 2005, the JUNG Project and the Regents of the University of
 * California All rights reserved.
 *
 * This software is open-source under the BSD license; see either "license.txt"
 * or http://jung.sourceforge.net/license.txt for a description.
 *
 * 
 */
package nature;

/**
 *
 */
import edu.uci.ics.jung.algorithms.layout3d.FRLayout;
import java.awt.BorderLayout;
import java.util.HashMap;
import java.util.Map;

import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JPanel;

import edu.uci.ics.jung.algorithms.layout3d.Layout;
import edu.uci.ics.jung.algorithms.layout3d.SpringLayout;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.TestGraphs;
import edu.uci.ics.jung.visualization.decorators.ToStringLabeller;
import edu.uci.ics.jung.visualization3d.VisualizationViewer;

import java.awt.Color;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.vecmath.Point3f;

/**
 *
 * @author Tom Nelson - tomnelson@dev.java.net
 *
 */
public class DisplayEmbedding3d extends JPanel {

//	Graph<String,Number> demoGraph = TestGraphs.getDemoGraph();
//	Graph<String,Number> oneComponentGraph = TestGraphs.getOneComponentGraph();
//	Map<String,Graph<String,Number>> graphMap = new HashMap<String,Graph<String,Number>>();
////	Map<String,Class> layoutMap = new HashMap<String,Class>();
//	JComboBox layoutBox, graphBox;
    Graph g;
    WeightedMajorizationDispl3d wm;
    double[][] coord;

    public DisplayEmbedding3d() {
        super(new BorderLayout());
        //setBackground(Color.yellow);
        //vv.setBackground(Color.yellow);
        //Graph<String,Number> graph = TestGraphs.getOneComponentGraph();
        IO ea = new IO();
        Graph<Integer, Integer> graph = ea.matlabToGraph("football.mat", "graph");
        FRLayout<Integer, Integer> layout = new FRLayout<Integer, Integer>(graph);

        //Layout<Integer, Integer> layout = new WeightedMajorizationDispl3d(graph);
        //   WeightedMajorizationDispl3d layout = new WeightedMajorizationDispl3d(graph);
        VisualizationViewer<Integer, Integer> vv = new VisualizationViewer<Integer, Integer>();
        double[][] coordRes = new double[graph.getVertexCount()][3];

        //vv.setBackground(Color.yellow);
        //Graph<String,Number> graph = TestGraphs.getOneComponentGraph();
        //TestGraphs.getDemoGraph();
        //vv.getRenderContext().setVertexStringer(new ToStringLabeller<String>());
        //Layout<String,Number> layout = new SpringLayout<String,Number>(graph);
        vv.setGraphLayout(layout);
        //vv.setBackground(Color.WHITE);
//        
//          try {
//            wait(100);
//        } catch (InterruptedException ex) {
//             Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
//        }

        add(vv);
        layout.setMaxIterations(5000);
        while (!layout.done()) {

        }

        System.out.println("done");
        for (Integer v : layout.getGraph().getVertices()) {
            Point3f p = layout.transform(v);
            coordRes[v][0] = p.x;
            coordRes[v][1] = p.y;
            coordRes[v][2] = p.z;

        }
        DataUtils du = new DataUtils();
        du.saveAsMatlab(coordRes, "coord", "FR3D.mat");

    }

    public DisplayEmbedding3d(Graph g, double[][] coordInit) {
        super(new BorderLayout());
        this.g = g;
        this.wm = new WeightedMajorizationDispl3d(g, coordInit);
        Layout<Integer, Integer> layout = wm;
        VisualizationViewer<Integer, Integer> vv = new VisualizationViewer<Integer, Integer>();
        vv.setGraphLayout(layout);
        add(vv);

    }

    public void display() {
        JFrame f = new JFrame();
        //f.setBackground(Color.yellow);
        f.add(this);
        f.setSize(600, 600);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setVisible(true);
        while (!wm.done()) {
//            String s = l.getStatus();
//            System.out.println(s);
//            double[][] cc = l.getCoordinates();
//            GraphCompression gcc = new GraphCompression(g, cc, 30);
//            System.out.println(gcc.mdlFunction());
//            try {
//               // wait(100);
//            } catch (InterruptedException ex) {
//                // Logger.getLogger(Visualization.class.getName()).log(Level.SEVERE, null, ex);
//            }
        }
        coord = wm.getCoordinates();

    }

    public double[][] getCoord() {
        return coord;
    }

    public static void main(String argv[]) {
        final DisplayEmbedding3d demo = new DisplayEmbedding3d();
        JFrame f = new JFrame();
        //f.setBackground(Color.yellow);
        f.add(demo);
        f.setSize(600, 600);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setVisible(true);
    }
}
