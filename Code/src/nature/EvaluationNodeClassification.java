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
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;

/**
 *
 * @author plantc59cs
 */
public class EvaluationNodeClassification {

    Graph g;
    double[][] coord;
    int[] labels;
    double[][][] confusionMatrices; //numExp x numClasses x numClasses (assigned)
    Instances nodes;
    int numClasses;
    int d;
    static int numExp = 100;
    boolean verbose = false;

    public EvaluationNodeClassification(Graph g, double[][] coord, int[] labels, int numClasses, int d) {
        this.g = g;
        this.coord = coord;
        this.labels = labels;
        this.numClasses = numClasses;
        this.d = d;
        confusionMatrices = new double[numExp][numClasses][numClasses];
        int numVar = d + 1;
        FastVector fvWekaAttributes = new FastVector(numVar);
        for (int i = 0; i < numVar - 1; i++) {
            Attribute a = new Attribute("d" + i);
            fvWekaAttributes.addElement(a);
        }
        FastVector my_nominal_values = new FastVector(numClasses);
        for (int i = 0; i < numClasses; i++) {
            my_nominal_values.addElement("c" + i);
        }
        Attribute cl = new Attribute("class", my_nominal_values);
        fvWekaAttributes.addElement(cl);
        nodes = new Instances("Graph", fvWekaAttributes, g.getVertexCount());
        nodes.setClassIndex(d);
        for (int i = 0; i < g.getVertexCount(); i++) {
            Instance inst = new Instance(nodes.numAttributes());
            for (int j = 0; j < nodes.numAttributes() - 1; j++) {
                inst.setValue((Attribute) fvWekaAttributes.elementAt(j), coord[i][j]);
            }
            inst.setValue((Attribute) fvWekaAttributes.elementAt(fvWekaAttributes.size() - 1), "c" + labels[i]);
            nodes.add(inst);
        }

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

    public double[] nodeClassification(double percentageLabeled) {
        int numLabeled = (int) Math.ceil(percentageLabeled * g.getVertexCount());
        int testDataSize = g.getVertexCount() - numLabeled;
        for (int i = 0; i < numExp; i++) {
            Instances train = new Instances(nodes);
            train.delete();
            Instances test = new Instances(nodes);
            Random r = new Random(i);
            for (int j = 0; j < numLabeled; j++) {
                int index = r.nextInt(test.numInstances());
                train.add(test.instance(index));
                test.delete(index);
            }
            if (verbose) {
                System.out.println("training data " + train.numInstances());
                System.out.println("test data " + test.numInstances());
            }
            //Classifier svm = new SMO();
            Classifier svm = new IBk(5);
            try {
                svm.buildClassifier(train);
                Evaluation eval = new Evaluation(train);
                eval.evaluateModel(svm, test);
                confusionMatrices[i] = eval.confusionMatrix();
                if (verbose) {
                    System.out.println(eval.toSummaryString("\nResults\n======\n", false));
                }
            } catch (Exception ex) {
                Logger.getLogger(EvaluationNodeClassification.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return microMacroF1Measure(testDataSize);

    }
    
    
   
    

    public double[][] nodeClassificationExperiments() {
        double[] pL = new double[7];
        pL[0] = 0.005;
        pL[1] = 0.01;
        pL[2] = 0.02;
        pL[3] = 0.04;
        pL[4] = 0.08;
        pL[5] = 0.16;
        pL[6] = 0.32;

//        pL[0] = 0.08;
//        pL[1] = 0.16;
//        pL[2] = 0.32;
//        pL[3] = 0.04;
//        pL[4] = 0.08;
//        pL[5] = 0.16;
//        pL[6] = 0.32;

        double[][] res = new double[7][2];
        for (int i = 0; i < pL.length; i++) {
            res[i] = nodeClassification(pL[i]);
        }
        return res;
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
