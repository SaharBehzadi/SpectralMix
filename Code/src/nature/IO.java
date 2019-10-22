/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import com.jmatio.io.MatFileReader;
import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLDouble;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Vector;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;

/**
 *
 * @author claudia
 */
public class IO {

    public Graph matlabToGraph(String dir, String variableName) {
        Graph<Integer, Integer> g = new UndirectedSparseGraph<Integer, Integer>();
        DataUtils du = new DataUtils();


        MatFileReader mfr = null;
        try {
            mfr = new MatFileReader(dir);
        } catch (IOException e) {
        }

        if (mfr != null) {
            double[][] data = ((MLDouble) mfr.getMLArray(variableName)).getArray();
            for (int i = 0; i < data.length; i++) {
                g.addVertex(i);

            }
            int counter = 0;
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < data.length; j++) {
                    if (i > j && data[i][j] == 1.0) {
                        //g.addEdge(du.getIndex(i, j, data.length), i, j);
                        g.addEdge(counter, i, j);
                        counter++;
                    }
                }
            }

        }
        return g;
    }
    
     public void writeGraphToMatlabSorted(Graph<Integer, Integer> g, int[] labels, String filename) {
        int numObj = g.getVertexCount();
        double[][] dist = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                if (g.isNeighbor(i, j)) {
                    dist[i][j] = 1.0;
                } else {
                    dist[i][j] = 0.0;
                }
            }
        }
        MLDouble q = new MLDouble("graph", dist);
        ArrayList ll = new ArrayList();
        ll.add(q);
        MatFileWriter mw = new MatFileWriter();
        try {
            String name = filename + ".mat";
            mw.write(name, ll);
        } catch (IOException ex) {
            //     Logger.getLogger(VI.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void writeGraphToMatlab(Graph<Integer, Integer> g, String filename) {
        int numObj = g.getVertexCount();
        double[][] dist = new double[numObj][numObj];
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numObj; j++) {
                if (g.isNeighbor(i, j)) {
                    dist[i][j] = 1.0;
                } else {
                    dist[i][j] = 0.0;
                }
            }
        }
        MLDouble q = new MLDouble("graph", dist);
        ArrayList ll = new ArrayList();
        ll.add(q);
        MatFileWriter mw = new MatFileWriter();
        try {
            String name = filename + ".mat";
            mw.write(name, ll);
        } catch (IOException ex) {
            //     Logger.getLogger(VI.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void writeFigure() {
//        double[][] gt = readMatlabMatrix("pr.mat", "gtCoordRed");
//        double[][] embedding = readMatlabMatrix("pr.mat", "Z_OwnBlowup");
        
        
         double[][] gt = readMatlabMatrix("twoMoonsProcrustes.mat", "groundTruth");
        double[][] embedding = readMatlabMatrix("twoMoonsProcrustes.mat", "speP");
        
        
        String fn = "resultSPE.txt";
        try {
            FileOutputStream fout = new FileOutputStream(new File(fn));
            for (int i = 0; i < gt.length; i++) {
                for (int j = 0; j < gt[i].length; j++) {
                    fout.write((gt[i][j] + " ").getBytes());
                }
                fout.write(("\n").getBytes());
            }
            fout.write(("\n").getBytes());
            fout.write(("\n").getBytes());
            for (int i = 0; i < embedding.length; i++) {
                for (int j = 0; j < embedding[i].length; j++) {
                    fout.write((embedding[i][j] + " ").getBytes());
                }
                fout.write(("\n").getBytes());
            }
            fout.write(("\n").getBytes());
            fout.write(("\n").getBytes());
            for (int i = 0; i < gt.length; i++) {
                for (int j = 0; j < gt[i].length; j++) {
                    fout.write((gt[i][j] + " ").getBytes());
                }
                fout.write(("\n").getBytes());
                for (int j = 0; j < embedding[i].length; j++) {
                    fout.write((embedding[i][j] + " ").getBytes());
                }
                fout.write(("\n").getBytes());
                fout.write(("\n").getBytes());

            }

            fout.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }


    }

    public void writeDoubleToMatlab(double[][] d, String variablename, String filename) {
        MLDouble q = new MLDouble(variablename, d);
        ArrayList ll = new ArrayList();
        ll.add(q);
        MatFileWriter mw = new MatFileWriter();
        try {
            String name = filename + ".mat";
            mw.write(name, ll);
        } catch (IOException ex) {
            //     Logger.getLogger(VI.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void writeDoubleToMatlab(double[][] d, String filename) {
        MLDouble q = new MLDouble(filename, d);
        ArrayList ll = new ArrayList();
        ll.add(q);
        MatFileWriter mw = new MatFileWriter();
        try {
            String name = filename + ".mat";
            mw.write(name, ll);
        } catch (IOException ex) {
            //     Logger.getLogger(VI.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public double[][] readMatlabMatrix(String dir, String variableName) {
        double[][] data = new double[1][1];

        MatFileReader mfr = null;
        try {
            mfr = new MatFileReader(dir);
        } catch (IOException e) {
        }

        if (mfr != null) {
            data = ((MLDouble) mfr.getMLArray(variableName)).getArray();
        }
        return data;
    }

    public Instances graphToArff(Graph g) {
        int numVar = g.getEdgeCount();
        FastVector fvWekaAttributes = new FastVector(numVar);
        for (int i = 0; i < numVar; i++) {
            FastVector my_nominal_values = new FastVector(2);
            my_nominal_values.addElement("true");
            my_nominal_values.addElement("false");
            Integer edgeIndex = (Integer) i;
            Pair endpoints = g.getEndpoints(edgeIndex);
            Attribute a = new Attribute(endpoints.toString(), my_nominal_values);
            fvWekaAttributes.addElement(a);
        }
        Instances dd = new Instances("Graph", fvWekaAttributes, g.getVertexCount());
        for (int i = 0; i < g.getVertexCount(); i++) {
            Instance inst = new Instance(dd.numAttributes());
            for (int j = 0; j < dd.numAttributes(); j++) {
                Integer EdgeIndex = (Integer) j;
                Integer NodeIndex = (Integer) i;
                if (g.isIncident(NodeIndex, EdgeIndex)) {
                    inst.setValue((Attribute) fvWekaAttributes.elementAt(j), "true");
                } else {
                    inst.setValue((Attribute) fvWekaAttributes.elementAt(j), "false");
                }

            }
            dd.add(inst);
        }
        ArffFileWriter af = new ArffFileWriter();
        af.saveFile("graph.arff", dd);
        return dd;
    }

    public void pingMatrixToArff(String dir, String variableName) {
        double[][] data = new double[1][1];
        MatFileReader mfr = null;
        try {
            mfr = new MatFileReader(dir);
        } catch (IOException e) {
        }

        if (mfr != null) {
            data = ((MLDouble) mfr.getMLArray(variableName)).getArray();
            int numPings = data.length;
            int numObj = data[0].length;
            FastVector fvWekaAttributes = new FastVector(numPings);
            for (int i = 0; i < numPings; i++) {
                //count values of Attribute i
                Vector<Integer> clusterIDs = new Vector<Integer>();
                for (int j = 0; j < numObj; j++) {
                    Integer ii = new Double(data[i][j]).intValue();
                    if (!clusterIDs.contains(ii)) {
                        clusterIDs.add(ii);
                    }
                }
                FastVector my_nominal_values = new FastVector(clusterIDs.size());
                for (int k = 0; k < clusterIDs.size(); k++) {
                    my_nominal_values.addElement(clusterIDs.elementAt(k).toString());
                }
                Attribute aa = new Attribute("clustering_" + i, my_nominal_values);
                fvWekaAttributes.addElement(aa);
            }
            Instances dd = new Instances("clusteringPings", fvWekaAttributes, numObj);
            for (int i = 0; i < numObj; i++) {
                Instance inst = new Instance(dd.numAttributes());
                for (int j = 0; j < numPings; j++) {
                    inst.setValue(j, data[j][i]);
                }
                dd.add(inst);
            }
            ArffFileWriter af = new ArffFileWriter();
            af.saveFileWithoutInstanceToString("clPings.arff", dd);

        }//mfr
    }

    public void matlabToArff(String dir, String variableName) {
        double[][] data = new double[1][1];

        MatFileReader mfr = null;
        try {
            mfr = new MatFileReader(dir);
        } catch (IOException e) {
        }

        if (mfr != null) {
            data = ((MLDouble) mfr.getMLArray(variableName)).getArray();



            FastVector fvWekaAttributes = new FastVector(data.length);
            for (int i = 0; i < data.length; i++) {
                fvWekaAttributes.addElement(new Attribute(i + "_numeric"));
            }
            Instances dd = new Instances("Rel", fvWekaAttributes, data.length);

            for (int i = 0; i < data.length; i++) {
                Instance inst = new Instance(dd.numAttributes());
                for (int j = 0; j < data[0].length; j++) {
                    inst.setValue((Attribute) fvWekaAttributes.elementAt(j), data[i][j]);
                }
                dd.add(inst);
            }
            ArffFileWriter af = new ArffFileWriter();
            af.saveFile("test.arff", dd);

        }

    }
}
