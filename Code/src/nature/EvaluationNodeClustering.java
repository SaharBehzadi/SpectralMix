/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.functions.SMO;
import weka.classifiers.lazy.IBk;
import weka.clusterers.ClusterEvaluation;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.clusterers.SimpleKMeans;

/**
 *
 * @author plantc59cs
 */
public class EvaluationNodeClustering {

    Graph g;
    double[][] coord;
    int[] labels;
    double[][][] confusionMatrices; //numExp x numClasses x numClasses (assigned)
    Instances nodes;
    int numClasses;
    int d;
    static int numExp = 100;
    boolean verbose = false;

    public EvaluationNodeClustering(Graph g, double[][] coord, int[] labels, int numClasses, int d) {
        this.g = g;
        this.coord = coord;
        this.labels = labels;
        this.numClasses = numClasses;
        this.d = d;
        int numVar = d;
        FastVector fvWekaAttributes = new FastVector(numVar);
        for (int i = 0; i < numVar - 1; i++) {
            Attribute a = new Attribute("d" + i);
            fvWekaAttributes.addElement(a);
        }
        FastVector my_nominal_values = new FastVector(numClasses);
        for (int i = 0; i < numClasses; i++) {
            int vv = i + 1;
            String s = Integer.toString(vv);
            my_nominal_values.addElement(s);
        }
        Attribute cl = new Attribute("class", my_nominal_values);
        //fvWekaAttributes.addElement(cl);
        nodes = new Instances("Graph", fvWekaAttributes, g.getVertexCount());
        //nodes.setClassIndex(d);
        for (int i = 0; i < g.getVertexCount(); i++) {
            Instance inst = new Instance(nodes.numAttributes());
            for (int j = 0; j < nodes.numAttributes() - 1; j++) {
                inst.setValue((Attribute) fvWekaAttributes.elementAt(j), coord[i][j]);
            }
            inst.setValue((Attribute) fvWekaAttributes.elementAt(fvWekaAttributes.size() - 1), labels[i]);
            nodes.add(inst);
        }

    }

    public void nodeClustering() {
        try {
            String[] options = new String[4];
            options[0] = "-init";                 // max. iterations
            options[1] = "1";
            options[2] = "-N";
            options[3] = "3";

            SimpleKMeans clusterer = new SimpleKMeans();   // new instance of clusterer
            clusterer.setOptions(options);     // set the options
            clusterer.setPreserveInstancesOrder(true);
            clusterer.buildClusterer(nodes);
            int[] ids = clusterer.getAssignments();

            double[][] idd = new double[ids.length][1];
            for (int i = 0; i < ids.length; i++) {
                idd[i][0] = ids[i];
            }
            DataUtils du = new DataUtils();
            du.saveAsMatlab(idd, "ids", "ids.mat");
        } catch (Exception ex) {
            Logger.getLogger(EvaluationNodeClustering.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    //as defined in http://rushdishams.blogspot.com/2011/08/micro-and-macro-average-of-precision.html
    //averaged over all classes
    public double[] microMacroF1Measure(int testSize) {
        double[] res = new double[2];
        double[] clSize = new double[numClasses];
        //double sumCl = 0.0;
        for (int i = 0; i < clSize.length; i++) {
            double sum = 0.0;
            for (int j = 0; j < clSize.length; j++) {
                sum += confusionMatrices[0][i][j] / (double) testSize;
            }
            clSize[i] = sum;
            //sumCl += sum;
        }
//        //Check--ok
//        System.out.println(sumCl);
//        //
        double[] sumTruePos = new double[numClasses];
        double[] sumTruePosFalsePos = new double[numClasses];
        double[] sumTruePosFalseNeg = new double[numClasses];
        double[] microF1 = new double[numClasses];
        double[] macroF1 = new double[numClasses];
        double[] avgPrecision = new double[numClasses];
        double[] avgRecall = new double[numClasses];
        for (int i = 0; i < numExp; i++) {
            for (int j = 0; j < numClasses; j++) {
                sumTruePos[j] += confusionMatrices[i][j][j];
                double truePosFalsePos = 0.0; //obj of this class and obj of other classes classified to it
                double truePosFalseNeg = 0.0; //all obj of this class
                for (int k = 0; k < numClasses; k++) {
                    sumTruePosFalsePos[j] += confusionMatrices[i][k][j];
                    truePosFalsePos += confusionMatrices[i][k][j];
                }
                for (int k = 0; k < numClasses; k++) {
                    sumTruePosFalseNeg[j] += confusionMatrices[i][j][k];
                    truePosFalseNeg += confusionMatrices[i][j][k];
                }
                if (confusionMatrices[i][j][j] > 0) {
                    avgPrecision[j] += confusionMatrices[i][j][j] / truePosFalsePos;
                    if (Double.isNaN(avgPrecision[j])) {
                        System.out.println("m");
                    }
                    avgRecall[j] += confusionMatrices[i][j][j] / truePosFalseNeg;
                }
            }
        }
        for (int j = 0; j < numClasses; j++) {

            double microPrec = 0.0;
            if (sumTruePosFalsePos[j] > 0) {
                microPrec = sumTruePos[j] / sumTruePosFalsePos[j];
            }
            double microRec = 0.0;
            if (sumTruePosFalseNeg[j] > 0) {
                microRec = sumTruePos[j] / sumTruePosFalseNeg[j];
            }
            double macroPrec = avgPrecision[j] / numExp;
            double macroRec = avgRecall[j] / numExp;
            if ((microPrec + microRec) > 0) {
                microF1[j] = 2.0 * (microPrec * microRec) / (microPrec + microRec);
            }

            //System.out.println(avgPrecision[j] + " " + avgRecall[j]);
            if ((macroPrec + macroRec) > 0) {
                macroF1[j] = 2.0 * (macroPrec * macroRec) / (macroPrec + macroRec);
            }

        }
        //return weighted average among classes
        for (int j = 0; j < numClasses; j++) {
            res[0] += clSize[j] * microF1[j];
            res[1] += clSize[j] * macroF1[j];
        }
        if (verbose) {
            System.out.println("Micro: " + res[0] + " " + " macro: " + res[1]);
        }
        return res;

    }

}
