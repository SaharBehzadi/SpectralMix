/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import Jama.Matrix;

/**
 *
 * @author claudia
 */
public class NumericCompression {

    double[][] orig; //original coordinates, numObj x dim
    double[][] proj; //projected coord., numObj x dim
    double[] individualCost; //cost per object
    int origDim;
    int projDim;
    int numObj;

    public NumericCompression(double[][] orig, double[][] proj) {
        this.orig = new Matrix(orig).transpose().getArrayCopy(); //since orig = numNumericalVariables x numObj
        this.proj = proj;
        this.origDim = this.orig[0].length;
        this.projDim = proj[0].length;
        this.numObj = proj.length;
        individualCost = new double[this.orig.length];
    }

    public double costWithoutEmbedding(){
        double[] mean = new double[origDim];
        for (int i = 0; i < numObj; i++) {
            for(int j = 0; j < origDim; j++)
            mean[j] += orig[i][j]/(double)numObj;
        }
        double[] var = new double[origDim];
         for (int i = 0; i < numObj; i++) {
            for(int j = 0; j < origDim; j++)
            var[j] += Math.pow((orig[i][j]-mean[j]), 2.0)/(double)numObj;
        }
        double sumEntropy = 0.0;
        for(int j = 0; j < origDim; j++)
            sumEntropy += lg2(Math.sqrt(2.0 * Math.PI * Math.E * var[j]));
        return numObj * sumEntropy;

    }


    public double getCost(int i, int j){
        return individualCost[i] * individualCost[j];
    }

    public double codingCost() {
        Matrix projM = new Matrix(proj);
        Matrix origM = new Matrix(orig);
        Matrix scProj = projM.transpose().times(projM);
        Matrix beta = scProj.inverse().times(projM.transpose()).times(origM);
        double[][] origEst = projM.times(beta).getArrayCopy();
        double[] err = new double[numObj];
        double meanErr = 0.0;
        for (int i = 0; i < numObj; i++) {
            err[i] = dist(origEst[i], orig[i]);
            meanErr += err[i] / (double) numObj;
        }
        double varErr = 0.0;
        for (int i = 0; i < numObj; i++) {
            varErr += Math.pow(err[i] - meanErr, 2) / (double) numObj;
        }
        double factor = 1.0/(Math.sqrt(2 * Math.PI * varErr));
        double sumCost = 0.0;
        double paramCost = paramCost();
        //distribute paramCost among all objects
        for(int i = 0; i < numObj; i++){
            double ll = factor * Math.exp(-Math.pow(err[i]-meanErr, 2)/(2 * varErr));
            individualCost[i] = lg2(1.0/ll);
            individualCost[i] += paramCost/(double)numObj;
            sumCost += individualCost[i];
        }
        double bla = numObj *  0.5 * lg2(2.0 * Math.PI * Math.E * varErr);
        double costProj = bla + paramCost;



        double costOrig = costWithoutEmbedding();
        return costProj - costOrig;
    }

    public double paramCost(){
        return (double)(numObj * projDim)/2.0 * lg2(numObj);
    }

    private double lg2(double a) {
        return Math.log(a) / Math.log(2);
    }

    private double dist(double[] a, double[] b) {
        double sum = 0.0;
        for (int i = 0; i < a.length; i++) {
            sum += (a[i] - b[i]) * (a[i] - b[i]);
        }
        return Math.sqrt(sum);
    }
}
