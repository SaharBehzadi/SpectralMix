/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import Jama.Matrix;
import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLDouble;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;
import java.util.Vector;
import weka.core.Instances;

/**
 *
 * @author claudia
 */
public class Princals {

    Instances data;
    int numObj;
    int numAttributes;
    double[][] fall;
    double[][] intercept;
    double[] meanL;
    double[] normL;
    double[] scalarL;
    int[] classLabel; //for visualization
    Matrix[] categoryCoord; //each matrix numCategories x dim
    Matrix objectCoord; //numObj x dim
    int projDim;
    int numNumeric;
    double objF;
    //static double cheat = 0.00001;
    static double convConst = 0.001;
    boolean verbose = true;
    boolean display = true;
    int iter;

    public Princals(Instances data, int projDim) {
        this.data = data;
        this.projDim = projDim;
        numObj = data.numInstances();
        numAttributes = data.numAttributes();
        numNumeric = 0;
        categoryCoord = new Matrix[numAttributes];
        objectCoord = new Matrix(numObj, projDim);
        //determine number of numerical attributes in the data set
        for (int i = 0; i < numAttributes; i++) {
            if (data.attribute(i).isNumeric()) {
                numNumeric++;
                categoryCoord[i] = new Matrix(numObj, projDim);
            } else {
                categoryCoord[i] = new Matrix(data.attribute(i).numValues(), projDim);
            }
        }
        fall = new double[numNumeric][projDim];
        intercept = new double[numNumeric][projDim];
        objF = Double.MAX_VALUE;
        iter = 0;
        classLabel = new int[data.numInstances()];

    }

    public void setClassLabel(double[] cl) {
        this.classLabel = new int[cl.length];
        for (int i = 0; i < cl.length; i++) {
            classLabel[i] = (int) cl[i];
        }
    }

//     double[][] fall;
//    double[][] intercept;
//    double[] meanL;
//    double[] normL;
//    double[] scalarL;
//    Matrix[] categoryCoord; //each matrix numCategories x dim
//    Matrix objectCoord; //numObj x dim
    private void zscaleResults() {
        for (int i = 0; i < fall.length; i++) {
            for (int j = 0; j < fall[i].length; j++) {
                fall[i][j] *= Math.sqrt(numObj - 1);
                intercept[i][j] *= Math.sqrt(numObj - 1);
            }
        }
        for (int i = 0; i < meanL.length; i++) {
            meanL[i] *= Math.sqrt(numObj - 1);
            normL[i] *= Math.sqrt(numObj - 1);

        }
        for (int i = 0; i < scalarL.length; i++) {
            scalarL[i] *= Math.sqrt(numObj - 1);
        }
        for (int i = 0; i < categoryCoord.length; i++) {
            categoryCoord[i] = categoryCoord[i].times(Math.sqrt(numObj - 1));
        }
        objectCoord = objectCoord.times(Math.sqrt(numObj - 1));


    }

    public Princals(Instances data, int projDim, double[][] initialObjCoord) {
        this.data = data;
        this.projDim = projDim;
        numObj = data.numInstances();
        numAttributes = data.numAttributes();
        numNumeric = 0;
        categoryCoord = new Matrix[numAttributes];
        objectCoord = new Matrix(initialObjCoord);
        //determine number of numerical attributes in the data set
        for (int i = 0; i < numAttributes; i++) {
            if (data.attribute(i).isNumeric()) {
                numNumeric++;
                categoryCoord[i] = new Matrix(numObj, projDim);
            } else {
                categoryCoord[i] = new Matrix(data.attribute(i).numValues(), projDim);
            }
        }
        fall = new double[numNumeric][projDim];
        intercept = new double[numNumeric][projDim];
        objF = Double.MAX_VALUE;
        iter = 0;
        classLabel = new int[data.numInstances()];
    }

//    //generate an additional orthogonal dim
//    public Princals(Instances data, double[] init, double[] scalar, int newDim) {
//        this.data = data;
//        this.projDim = 1;
//        scalarL = new double[(newDim * newDim - newDim) / 2];
//        for (int i = 0; i < scalar.length; i++) {
//            scalarL[i] = scalar[i];
//        }
//        numObj = data.numInstances();
//        numAttributes = data.numAttributes();
//        numNumeric = 0;
//        categoryCoord = new Matrix[numAttributes];
//        double[][] initm = new double[init.length][1];
//        for (int i = 0; i < init.length; i++) {
//            initm[i][0] = init[i];
//        }
//        objectCoord = new Matrix(initm);
//        //determine number of numerical attributes in the data set
//        for (int i = 0; i < numAttributes; i++) {
//            if (data.attribute(i).isNumeric()) {
//                numNumeric++;
//                categoryCoord[i] = new Matrix(numObj, projDim);
//            } else {
//                categoryCoord[i] = new Matrix(data.attribute(i).numValues(), projDim);
//            }
//        }
//        fall = new double[numNumeric][projDim];
//        intercept = new double[numNumeric][projDim];
//        objF = Double.MAX_VALUE;
//    }
    public void runAdditionalColumn(double[][] existingCoord) {
        initialize();
        boolean converged = false;
        while (!converged) {
            updateCategories();
            updateObjectsOrtEx(existingCoord);
            converged = checkConvergence();
            iter++;
            if (verbose) {
                System.out.println(iter + " " + objF);
                writeTest(objectCoord.getArrayCopy(), "obj.mat");
            }

        }
        //zscaleResults();
        if (display) {
            categoryPlot();
            objectPlot();
            combinedPlot();
            writeTest(objectCoord.getArrayCopy(), "obj.mat");
            writeTest(categoryCoord[0].getArrayCopy(), "bla1.mat");
            writeTest(categoryCoord[1].getArrayCopy(), "bla2.mat");
            writeTest(categoryCoord[2].getArrayCopy(), "bla3.mat");
            writeTest(categoryCoord[3].getArrayCopy(), "bla4.mat");
            writeTest(categoryCoord[4].getArrayCopy(), "bla5.mat");

        }

    }

    public Princals() {
        iter = 0;
    }

    public int getIter() {
        return iter;
    }

    public double[] getScalarL() {
        return scalarL;
    }

    public double[] getNormL() {
        return normL;
    }

    public Matrix[] getCategoryCoord() {
        return categoryCoord;
    }

    public double[][] getFall() {
        return fall;
    }

    public double[] getMeanL() {
        return meanL;
    }

    public double[][] getIntercept() {
        return intercept;
    }

    public Matrix getObjectCoord() {
        return objectCoord;
    }

    public boolean runWithInit(double[][] init) {
        boolean success = true;
        initializeAdditionalDim(init);
        boolean converged = false;
        while (!converged) {
            updateCategories();
            updateObjects();
            iter++;
            if (projDim > 1) {
                for (int i = 0; i < normL.length; i++) {
                    if (normL[i] < 0.0000000001) {
                        success = false;
                    }
                }
            }
            converged = checkConvergence();
            if (verbose) {
                System.out.println(iter + " " + objF);
                writeTest(objectCoord.getArrayCopy(), "obj.mat");
            }

        }
        //zscaleResults();
        if (display) {
            categoryPlot();
            objectPlot();
            combinedPlot();
            writeTest(objectCoord.getArrayCopy(), "obj.mat");
//            writeTest(categoryCoord[0].getArrayCopy(), "bla1.mat");
//            writeTest(categoryCoord[1].getArrayCopy(), "bla2.mat");
//            writeTest(categoryCoord[2].getArrayCopy(), "bla3.mat");
//            writeTest(categoryCoord[3].getArrayCopy(), "bla4.mat");
//            writeTest(categoryCoord[4].getArrayCopy(), "bla5.mat");

        }
        return success;
    }

    public boolean runWithoutInit() {
        boolean converged = false;
        boolean success = true;
        while (!converged) {
            updateCategories();
            updateObjects();
            iter++;
            for (int i = 0; i < normL.length; i++) {
                if (normL[i] < 0.0000000001) {
                    success = false;
                }
            }
            converged = checkConvergence();
            if (verbose) {
                System.out.println(iter + " " + objF);
                writeTest(objectCoord.getArrayCopy(), "obj.mat");
            }

        }
        //zscaleResults();
        if (display) {
            categoryPlot();
            objectPlot();
            combinedPlot();
            writeTest(objectCoord.getArrayCopy(), "obj.mat");
//            writeTest(categoryCoord[0].getArrayCopy(), "bla1.mat");
//            writeTest(categoryCoord[1].getArrayCopy(), "bla2.mat");
//            writeTest(categoryCoord[2].getArrayCopy(), "bla3.mat");
//            writeTest(categoryCoord[3].getArrayCopy(), "bla4.mat");
//            writeTest(categoryCoord[4].getArrayCopy(), "bla5.mat");

        }
        return success;

    }

    public void forPresentation() {
        initialize();
        updateCategories();
        updateObjects();
        updateCategories();
        categoryPlot();
        objectPlot();
        combinedPlot();
        writeTest(objectCoord.getArrayCopy(), "obj.mat");
//        writeTest(categoryCoord[0].getArrayCopy(), "bla1.mat");
//        writeTest(categoryCoord[1].getArrayCopy(), "bla2.mat");
//        writeTest(categoryCoord[2].getArrayCopy(), "bla3.mat");
//        writeTest(categoryCoord[3].getArrayCopy(), "bla4.mat");
//        writeTest(categoryCoord[4].getArrayCopy(), "bla5.mat");

    }

    public boolean run() {
        boolean success = true;
        initialize();
        boolean converged = false;
        while (!converged && success) {
            updateCategories();
            updateObjects();
            iter++;
            converged = checkConvergence();
            for (int i = 0; i < normL.length; i++) {
                if (normL[i] < 0.0000000001) {
                    updateCategories();
                    success = false;
                }
                if (Double.isNaN(objF)) {
                    success = false;
                }
            }
            if (verbose) {
                System.out.println(iter + " " + objF);
                writeTest(objectCoord.getArrayCopy(), "obj.mat");
            }

        }
        //zscaleResults();
        if (display) {
            categoryPlot();
            objectPlot();
            combinedPlot();

//            writeTest(categoryCoord[0].getArrayCopy(), "bla1.mat");
//            writeTest(categoryCoord[1].getArrayCopy(), "bla2.mat");
//            writeTest(categoryCoord[2].getArrayCopy(), "bla3.mat");
//            writeTest(categoryCoord[3].getArrayCopy(), "bla4.mat");
//            writeTest(categoryCoord[4].getArrayCopy(), "bla5.mat");
//            writeTest(categoryCoord[5].getArrayCopy(), "bla6.mat");
//            writeTest(categoryCoord[6].getArrayCopy(), "bla7.mat");
//            writeTest(categoryCoord[7].getArrayCopy(), "bla8.mat");


        }
        return success;

    }

    private void writeTest(double[][] d, String filename) {
        //write burt matrix to matlab
        MLDouble t = new MLDouble("test", d);
        ArrayList ll = new ArrayList();
        ll.add(t);
        MatFileWriter mw = new MatFileWriter();
        try {
            String name = filename;
            mw.write(name, ll);
        } catch (IOException ex) {
        }
    }

    //objectCoord (numObj x oldDim] : previous object locations
    public void initializeAdditionalDim(double[][] objCoord) {
        double[][] start = new double[numObj][projDim];
        Random r = new Random(1);
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < projDim - 1; j++) {
                start[i][j] = objCoord[i][j];
            }
            start[i][projDim - 1] = r.nextGaussian();
        }
        //center last column
        double colMean = 0.0;
        for (int j = 0; j < numObj; j++) {
            colMean += start[j][projDim - 1] / (double) numObj;
        }
        for (int j = 0; j < numObj; j++) {
            start[j][projDim - 1] -= colMean;
        }
        start = modifiedGramSchmidt(start);
        // writeTest(start, "startAddInit.mat");

        for (int j = 0; j < numObj; j++) {
            start[j][projDim - 1] *= Math.sqrt(numObj - 1);
        }

        objectCoord = new Matrix(start);

    }

    public void initialize() {
        double[][] start = new double[numObj][projDim];
        Random r = new Random(1);
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < projDim; j++) {
                start[i][j] = r.nextGaussian();
            }
        }
//standardize start
//        DataTransformer dt = new DataTransformer();
//        Instances d = dt.create(start);
//        d = dt.zScale(d);
//        start = dt.extract(d);
        //(start, "start1.mat");

        //center each column
        double[] colMean = new double[projDim];
        for (int i = 0; i < projDim; i++) {
            for (int j = 0; j < numObj; j++) {
                colMean[i] += start[j][i] / (double) numObj;
            }
        }
        for (int i = 0; i < projDim; i++) {
            for (int j = 0; j < numObj; j++) {
                start[j][i] -= colMean[i];
            }
        }
        //writeTest(start, "start2.mat");

        //colunm-orthogonalize
        start = modifiedGramSchmidt(start);

        for (int i = 0; i < projDim; i++) {
            for (int j = 0; j < numObj; j++) {
                start[j][i] *= Math.sqrt(numObj - 1);
            }
        }
        // writeTest(start, "start3.mat");

        if (verbose) {
            writeTest(start, "start.mat");
        }
        objectCoord = new Matrix(start);
    }

    private void updateCategories() {
        //set category coordinates to zero
        for (int i = 0; i < numAttributes; i++) {
            categoryCoord[i] = new Matrix(categoryCoord[i].getRowDimension(), categoryCoord[i].getColumnDimension());
        }
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numAttributes; j++) {
                if (data.attribute(j).isNumeric()) { //set category to the location of object i
                    for (int k = 0; k < projDim; k++) {
                        categoryCoord[j].set(i, k, objectCoord.get(i, k));
                    }
                } else { //set category location to the mean of all corresponding objects
                    int[] count = data.attributeStats(j).nominalCounts;
                    for (int k = 0; k < projDim; k++) {
                        int categoryIndex = (int) data.instance(i).value(j);
                        double oldValue = categoryCoord[j].get(categoryIndex, k);
                        double newValue = oldValue + (objectCoord.get(i, k) / (double) count[categoryIndex]);
                        categoryCoord[j].set(categoryIndex, k, newValue);
                    }
                }
            }
        }
        //correct the numerical category locations by regression
        int counter = 0;
        for (int i = 0; i < numAttributes; i++) {
            if (data.attribute(i).isNumeric()) {
                double[] originalValues = new double[numObj];
                double[] transformedCoord = new double[numObj];
                double mean_orig = 0.0;
                for (int j = 0; j < numObj; j++) {
                    originalValues[j] = data.instance(j).value(i);
                    mean_orig += originalValues[j] / (double) numObj;
                }
                double[][] transformedCoordNew = new double[numObj][projDim];
//                //DEBUG
//                double[][] transformedTest = new double[numObj][projDim];
//                        //DEBUG
                for (int j = 0; j < projDim; j++) {
                    double[][] d = categoryCoord[i].getArrayCopy();
                    double mean_transformed = 0;
                    for (int k = 0; k < numObj; k++) {
                        transformedCoord[k] = d[k][j];
                        mean_transformed += transformedCoord[k] / (double) numObj;
                    }
                    //perform a linear regression with transformedCoord as dependent variable y and originalValues as variable x
                    double upper = 0;
                    double lower = 0;
                    for (int k = 0; k < numObj; k++) {
                        upper += (originalValues[k] - mean_orig) * (transformedCoord[k] - mean_transformed);
                        lower += (originalValues[k] - mean_orig) * (originalValues[k] - mean_orig);
                    }
                    double b = upper / lower;
                    double a = mean_transformed - b * mean_orig;
                    fall[counter][j] = b;
                    intercept[counter][j] = a;

                    for (int k = 0; k < numObj; k++) {
                        transformedCoordNew[k][j] = a + b * originalValues[k];

//                        //test
//                        transformedTest[k][j] =  fall[counter][j] * originalValues[k] + intercept[counter][j] + regError[k][counter][j]; //ok: transformedTest = transformed.
//                        System.out.println("transformedTest: " + transformedTest[k][j] + " transformedCoord: " + transformedCoord[k]);
//                        //test
                    }

                }
                counter++;
                categoryCoord[i] = new Matrix(transformedCoordNew);
            }
        }
    }

    //writes scalarL and norm: scalarproducts of new vector to existing ones, norm of new vector.
    private void updateObjectsOrtEx(double[][] existingCoord) {
        double[][] newLoc = new double[numObj][projDim];
        for (int i = 0; i < numObj; i++) {
//            //DEBUG
//            if(i == 101)
//                System.out.println("m");
//            //DEBUG
            for (int j = 0; j < numAttributes; j++) {
                for (int k = 0; k < projDim; k++) {
                    if (data.attribute(j).isNumeric()) {
                        //DEBUG
                        double singleLoc = categoryCoord[j].get(i, k);
                        double under = (double) numAttributes;
                        //DEBUG
                        newLoc[i][k] += categoryCoord[j].get(i, k) / (double) numAttributes;
                    } else {
                        //DEBUG
                        double singleLoc = categoryCoord[j].get((int) data.instance(i).value(j), k);
                        double under = (double) numAttributes;
                        //DEBUG
                        newLoc[i][k] += categoryCoord[j].get((int) data.instance(i).value(j), k) / (double) numAttributes;
                    }
                } //k
            }
        }
        //center each column
        double[] colMean = new double[projDim];
        for (int i = 0; i < projDim; i++) {
            for (int j = 0; j < numObj; j++) {
                colMean[i] += newLoc[j][i] / (double) numObj;
            }
        }
        meanL = colMean;
        for (int i = 0; i < projDim; i++) {
            for (int j = 0; j < numObj; j++) {
                newLoc[j][i] -= colMean[i];
            }
        }
        // orthogonalize new vector to the old ones
        int oldDim = existingCoord[0].length;

        normL = new double[1];
        scalarL = new double[oldDim];
        int counter = 0;
        ///  for (int j = 0; j < k; j++) {
        for (int i = 0; i < oldDim; i++) {
            double skalarprod = 0.0;
            double self_i = 0.0;
            double proj_vi_vj = 0.0;
            for (int l = 0; l < numObj; l++) {
                skalarprod += existingCoord[l][i] * newLoc[l][0];
                self_i += existingCoord[l][i] * existingCoord[l][i];
            }
            scalarL[i] = skalarprod / self_i;
            counter++;
            for (int l = 0; l < numObj; l++) {
                proj_vi_vj = (skalarprod / self_i) * existingCoord[l][i];
                //scalarL[j] += proj_vi_vj;
                newLoc[l][0] = newLoc[l][0] - proj_vi_vj;
            }
        } //i
        double norm_j = 0.0;
        for (int l = 0; l < numObj; l++) {
            norm_j += newLoc[l][0] * newLoc[l][0];
        }
        norm_j = Math.sqrt(norm_j);
        normL[0] = norm_j;

        for (int l = 0; l < numObj; l++) {
            newLoc[l][0] = newLoc[l][0] / norm_j;
            newLoc[l][0] = newLoc[l][0] * Math.sqrt(numObj - 1);
        }
        objectCoord = new Matrix(newLoc);





    }

    //set objects to the mean of their category locations
    private void updateObjects() {
        double[][] newLoc = new double[numObj][projDim];
        for (int i = 0; i < numObj; i++) {
//            //DEBUG
//            if(i == 101)
//                System.out.println("m");
//            //DEBUG
            for (int j = 0; j < numAttributes; j++) {
                for (int k = 0; k < projDim; k++) {
                    if (data.attribute(j).isNumeric()) {
                        //DEBUG
                        double singleLoc = categoryCoord[j].get(i, k);
                        double under = (double) numAttributes;
                        //DEBUG
                        newLoc[i][k] += categoryCoord[j].get(i, k) / (double) numAttributes;
                    } else {
                        //DEBUG
                        double singleLoc = categoryCoord[j].get((int) data.instance(i).value(j), k);
                        double under = (double) numAttributes;
                        //DEBUG
                        newLoc[i][k] += categoryCoord[j].get((int) data.instance(i).value(j), k) / (double) numAttributes;
                    }
                } //k
            }
        }
        //center each column
        double[] colMean = new double[projDim];
        for (int i = 0; i < projDim; i++) {
            for (int j = 0; j < numObj; j++) {
                colMean[i] += newLoc[j][i] / (double) numObj;
            }
        }
        meanL = colMean;
        for (int i = 0; i < projDim; i++) {
            for (int j = 0; j < numObj; j++) {
                newLoc[j][i] -= colMean[i];
            }
        }
        //writeTest(newLoc, "beforeOrt.mat");
        //colunm-orthogonalize
        newLoc = modifiedGramSchmidt(newLoc);
        //writeTest(newLoc, "afterOrt.mat");
////        //z-scale
//        for (int i = 0; i < projDim; i++) {
//            for (int j = 0; j < numObj; j++) {
//                newLoc[j][i] *= Math.sqrt(numObj - 1);
//            }
//        }
        for (int i = 0; i < projDim; i++) {
            for (int j = 0; j < numObj; j++) {
                newLoc[j][i] *= Math.sqrt(numObj);
            }
        }
        objectCoord = new Matrix(newLoc);
    }

    //sum up squared difference between object and category locations.
    private boolean checkConvergence() {
        int size = numObj * projDim;
        double objNew = 0;
        for (int i = 0; i < numObj; i++) {
            for (int j = 0; j < numAttributes; j++) {
                for (int k = 0; k < projDim; k++) {
                    if (data.attribute(j).isNumeric()) {
                        objNew += Math.pow(objectCoord.get(i, k) - categoryCoord[j].get(i, k), 2);

                    } else {
                        objNew += Math.pow(objectCoord.get(i, k) - categoryCoord[j].get((int) data.instance(i).value(j), k), 2);

                    }
                }
            }

        }
        //objNew /= (double) numAttributes;
        objNew /= size;
        if (objF - objNew < -0.1) {
            if (verbose) {
                System.err.println("increase of objF: " + (objF - objNew));
            }
            //objF = objNew;
            return true;
        }
        if (objF - objNew < convConst) {
            objF = objNew;

            return true;
        } else {
            objF = objNew;
            return false;
        }

    }

    //orthogonalize the column vectors in v
    public double[][] modifiedGramSchmidt(double[][] v) {
        int k = v[0].length;
        int vLength = v.length;
        normL = new double[projDim];
        scalarL = new double[(projDim * projDim - projDim) / 2];
        int counter = 0;
        for (int j = 0; j < k; j++) {
            for (int i = 0; i < j; i++) {
                double skalarprod = 0.0;
                double self_i = 0.0;
                double proj_vi_vj = 0.0;
                for (int l = 0; l < vLength; l++) {
                    skalarprod += v[l][i] * v[l][j];
                    self_i += v[l][i] * v[l][i];
                }
                scalarL[counter] = skalarprod / self_i;
                counter++;
                for (int l = 0; l < vLength; l++) {
                    proj_vi_vj = (skalarprod / self_i) * v[l][i];
                    //scalarL[j] += proj_vi_vj;
                    v[l][j] = v[l][j] - proj_vi_vj;
                }
            } //i
            double norm_j = 0.0;
            for (int l = 0; l < vLength; l++) {
                norm_j += v[l][j] * v[l][j];
            }
            norm_j = Math.sqrt(norm_j);
            normL[j] = norm_j;

            for (int l = 0; l < vLength; l++) {
                v[l][j] = v[l][j] / norm_j;
            }
        }//j
        //writeTest(v, "orth.mat");
        return v;
    }

    private void objectPlot() {
        double[][] x = objectCoord.getArrayCopy();
        DataObject[] d = new DataObject[numObj];
        for (int i = 0; i < numObj; i++) {
            StringBuffer sb = new StringBuffer();
            for (int j = 0; j < numAttributes; j++) {
                if (data.attribute(j).isNumeric()) {
                    sb.append(Double.toString(data.instance(i).value(j)) + " ");
                }
                if (data.attribute(j).isNominal()) {
                    sb.append(data.instance(i).stringValue(j) + " ");
                }
            }
            for (int j = 0; j < x[i].length; j++) {
                sb.append(new Double(x[i][j]).toString() + " ");
            }
            d[i] = new DataObject(x[i], sb.toString(), i);
//            d[i].classID = (int) concepts.instance(i).value(0); //concept 1 :categorical
            d[i].numericClassID = 1.0;
            d[i].classID = classLabel[i];
            //d[i].numericClassID = concepts.instance(i).value(0);


        }
        for (int i = 0; i < x[0].length; i++) {
            for (int j = 0; j < x[0].length; j++) {
                if (i < j) {
                    String title = "object plot: " + i + " " + j;
                    VisuInfo vi = new VisuInfo(d, title, i, j);
                    vi.setSize(600, 600);
                    vi.setLocation(100, 0);
                    vi.setVisible(true);
                }
            }
        }
    }

    private void combinedPlot() {
        double[][] x = objectCoord.getArrayCopy();
        DataObject[] d = new DataObject[numObj];
        for (int i = 0; i < numObj; i++) {
            StringBuffer sb = new StringBuffer();
            for (int j = 0; j < numAttributes; j++) {
                if (data.attribute(j).isNumeric()) {
                    sb.append(Double.toString(data.instance(i).value(j)) + " ");
                }
                if (data.attribute(j).isNominal()) {
                    sb.append(data.instance(i).stringValue(j) + " ");
                }
            }
            for (int j = 0; j < x[i].length; j++) {
                sb.append(new Double(x[i][j]).toString() + " ");
            }
            d[i] = new DataObject(x[i], sb.toString(), i);
            d[i].classID = 100; //concept 1 :categorical
            d[i].numericClassID = 1.0;
            //d[i].numericClassID = concepts.instance(i).value(0);


        }

        Vector<DataObject> ddv = new Vector<DataObject>();
        //under construction

        int counter = 0;
        for (int i = 0; i < numAttributes; i++) {
            if (data.attribute(i).isNumeric()) {
                for (int j = 0; j < numObj; j++) {
                    double[] cc = new double[projDim];
                    for (int k = 0; k < projDim; k++) {
                        cc[k] = categoryCoord[i].get(j, k);
                    }

                    //public DataObject(double[] coord, String label, int clusterID, double numericClassID) {
                    //public DataObject(double[] coord, String label, int classID, double numericClassID, int number) {
                    StringBuffer sb = new StringBuffer();
                    sb.append(Double.toString(data.instance(j).value(i)) + " ");
                    for (int l = 0; l < cc.length; l++) {
                        sb.append(new Double(cc[l]).toString() + " ");
                    }
                    ddv.add(new DataObject(cc, sb.toString(), i, 1.0, counter));
                    counter++;
                }

            } else {
                for (int j = 0; j < data.attribute(i).numValues(); j++) {
                    double[] cc = new double[projDim];
                    for (int k = 0; k < projDim; k++) {
                        cc[k] = categoryCoord[i].get(j, k);
                    }
                    ddv.add(new DataObject(cc, data.attribute(i).toString() + " " + j, i, 1.0, counter));
                    counter++;
                }

            }

        }
        DataObject[] dd = new DataObject[ddv.size()];
        for (int i = 0; i < ddv.size(); i++) {
            dd[i] = ddv.elementAt(i);
        }


//        DataObject[] dd = new DataObject[categories.length];
//        int counter = 0;
//        for (int i = 0; i < numAttributes; i++) {
//            for (int j = 0; j < g[i].colNames.length; j++) {
//                dd[counter] = new DataObject(categories[counter], data.attribute(i).toString() + " " + g[i].colNames[j], counter);
//                // dd[counter].numericClassID = 1.0;
//                dd[counter].classID = 100;
//                counter++;
//            }
//        }
        DataObject[] combined = new DataObject[d.length + dd.length];
        for (int i = 0; i < d.length; i++) {
            combined[i] = new DataObject(d[i]);
            combined[i].numericClassID = 1.0;
        }
        for (int i = 0; i < dd.length; i++) {
            combined[d.length + i] = new DataObject(dd[i]);
            combined[d.length + i].numericClassID = 1.0;
        }

        for (int i = 0; i < x[0].length; i++) {
            for (int j = 0; j < x[0].length; j++) {
                if (i < j) {
                    String title = "combined plot: " + i + " " + j;
                    VisuInfo vi = new VisuInfo(combined, title, i, j);
                    vi.setSize(600, 600);
                    vi.setLocation(100, 0);
                    vi.setVisible(true);
                }
            }
        }
    }

    private void categoryPlot() {
        Vector<DataObject> dd = new Vector<DataObject>();
        //under construction

        int counter = 0;
        for (int i = 0; i < numAttributes; i++) {
            if (data.attribute(i).isNumeric()) {
                for (int j = 0; j < numObj; j++) {
                    double[] cc = new double[projDim];
                    for (int k = 0; k < projDim; k++) {
                        cc[k] = categoryCoord[i].get(j, k);
                    }

                    //public DataObject(double[] coord, String label, int clusterID, double numericClassID) {
                    //public DataObject(double[] coord, String label, int classID, double numericClassID, int number) {
                    dd.add(new DataObject(cc, data.attribute(i).toString() + " " + Double.toString(data.instance(j).value(i)), i, 1.0, counter));
                    counter++;
                }

            } else {
                for (int j = 0; j < data.attribute(i).numValues(); j++) {
                    double[] cc = new double[projDim];
                    for (int k = 0; k < projDim; k++) {
                        cc[k] = categoryCoord[i].get(j, k);
                    }
                    //public DataObject(double[] coord, String label, int classID, double numericClassID, int number) {
                    double size = 0.0;
                    if(j == 0)
                        size = 0.0;
                    else
                        size = 1.0;
                    dd.add(new DataObject(cc, data.attribute(i).toString() + " " + j, j, size, counter));
                    counter++;
                }

            }

        }
        DataObject[] d = new DataObject[dd.size()];
        for (int i = 0; i < dd.size(); i++) {
            d[i] = dd.elementAt(i);
        }
        for (int i = 0; i < projDim; i++) {
            for (int j = 0; j < projDim; j++) {
                if (i < j) {
                    String title = "category plot: " + i + " " + j;
                    VisuInfo vi = new VisuInfo(d, title, i, j);
                    vi.setSize(600, 600);
                    vi.setLocation(100, 0);
                    vi.setVisible(true);
                }
            }//j
        }
    }
}
