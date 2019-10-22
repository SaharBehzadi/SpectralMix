/*
 * Ica.java
 *
 * Created on 13. Oktober 2005, 15:46
 *
 */
package nature;

import Jama.*;
import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLDouble;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author boehm
 */
public class Ica {

    public double[][] ev;   //eigenvalues
    public double wscCost;  // Cost of White Space Compensation in Bits
    double[][] icaInput;
    double[][] centered;
    double[][] whitened;
    double[][] projPCA;
    double[][] ics;
    double[] mean;
    Matrix E;  //Eigenvectors
    Matrix D; //Eigenvalues
    Matrix reducedE; //dimensionality reduced matrices
    Matrix reducedD;
    Matrix whitening;
    Matrix dewhitening;
    Matrix A; //mixing matrix
    Matrix W; //inverse of mixing matrix
    boolean verbose = true;
    boolean minimalOutput = false;

    /**
     * Creates a new instance of Ica
     */
    public Ica() {
    }

    public Ica(double[][] icaInput){
        this.icaInput = icaInput;
    }
    
    public Matrix matrix(double[][] m) {
        return new Matrix(m);
    }

    public Matrix columnVector(double[] v) {
        return new Matrix(v, v.length);
    }

    public Matrix rowVector(double[] v) {
        return new Matrix(v, 1);
    }

    public double[][] getEV() {
        return this.ev;
    }
    
    public void doPca(int numDimensions){
      if (!minimalOutput) {
            dataCheck(icaInput, "icainput");
        }
        if (E == null) {
            pca(icaInput);
        }
        reduceDimensionality(numDimensions);   
        projPCA = new double[numDimensions][icaInput[0].length];
        projPCA = reducedE.transpose().times(new Matrix(icaInput)).getArrayCopy();
    }
          

    //first pca (if required), then dimensionality reduction to numdimensions then ica
    public void doIca(int numDimensions) {
        if (!minimalOutput) {
            dataCheck(icaInput, "icainput");
        }
        if (E == null) {
            pca(icaInput);
        }
        reduceDimensionality(numDimensions);
        whitening();
        //fasticaHv(numDimensions); FastICA nach Matlab/Hyverinnen
        double[][] w = optimizeEntropy(new Matrix(whitened).transpose().getArrayCopy());
        //to DO CHECK RESULT
//        OK
//        dataCheck(centered, "centered");
//        dataCheck(E.getArrayCopy(), "eigenvectors");
//        dataCheck(D.getArrayCopy(), "eigenvalues");
//        dataCheck(whitening.getArrayCopy(), "whitening");
//        dataCheck(dewhitening.getArrayCopy(), "dewhitening");
        if (!minimalOutput) {
            dataCheck(whitened, "whitened");
        }
        W = new Matrix(w).times(whitening); //4 x 400
        A = dewhitening.times(new Matrix(w).transpose()); //400 x 4
        //compute ICs. add mean back to data
        //icasig = W * mixedsig + (W * mixedmean) * ones(1, NumOfSampl);
        // Matrix icm = W.times(new Matrix(icaInput)); //4 x 1600 - ok
        double[][] ones = new double[1][icaInput[0].length];
        for (int i = 0; i < ones.length; i++) {
            for (int j = 0; j < ones[i].length; j++) {
                ones[i][j] = 1;
            }
        }
        Matrix tt = W.times(rowVector(mean).transpose()).times(new Matrix(ones)); //4 x 1600
        Matrix icm = W.times(new Matrix(icaInput)).plus(tt);
        ics = icm.getArrayCopy();
        //normalizeICs();
        if (!minimalOutput) {
            dataCheck(ics, "ics");
        }
        //clusterICs(numDimensions);
    }

    //writes whitening, dewhitening, whitened
    private void whitening() {
        //whitenig = sqrt(D)^-1 * E'
        double[][] d = reducedD.getArrayCopy();
        for (int i = 0; i < d.length; i++) {
            d[i][i] = 1.0 / Math.sqrt(d[i][i]);
        }
        Matrix dm = new Matrix(d);
        whitening = dm.times(reducedE.transpose());
        //dewhiteningMatrix = E * sqrt (D);
        d = reducedD.getArrayCopy();
        for (int i = 0; i < d.length; i++) {
            d[i][i] = Math.sqrt(d[i][i]);
        }
        dm = new Matrix(d);
        dewhitening = reducedE.times(dm);
        //whitenedData =  whitening * data;
        whitened = whitening.times(new Matrix(centered)).getArrayCopy();

    }

    public double[][] centerData(double[][] data) {
        int numObj = data.length;
        int d = data[0].length;
        double[][] result = new double[numObj][d];

        for (int i = 0; i < d; i++) {
            double sum = 0.0;
            for (int j = 0; j < numObj; j++) {
                sum += data[j][i];
            }
            sum /= numObj;
            for (int j = 0; j < numObj; j++) {
                result[j][i] = data[j][i] - sum;
            }
        }

        return result;
    }

    public double[][] centerData(double[][] data, double[] center) {
        int numObj = data.length;
        int d = data[0].length;
        if (center == null) {
            center = new double[d];
        }
        double[][] result = new double[numObj][d];

        for (int i = 0; i < d; i++) {
            center[i] = 0.0;
            for (int j = 0; j < numObj; j++) {
                center[i] += data[j][i];
            }
            center[i] /= numObj;
            for (int j = 0; j < numObj; j++) {
                result[j][i] = data[j][i] - center[i];
            }
        }

        return result;
    }

    public double[][] centerByMedian(double[][] data) {
        int numObj = data.length;
        int dim = data[0].length;
        double[] sortarray = new double[numObj];
        double[][] result = new double[numObj][dim];
        for (int i = 0; i < dim; i++) {
            double median = 0.0;
            for (int j = 0; j < numObj; j++) {
                sortarray[j] = data[j][i];
            }
            Arrays.sort(sortarray);
            if (numObj % 2 == 0) {
                median = (sortarray[numObj / 2 - 1] + sortarray[numObj / 2]) / 2;
            } else {
                median = (sortarray[numObj / 2 - 1] + sortarray[numObj / 2] + sortarray[numObj / 2 + 1]) / 3;
            }
            for (int j = 0; j < numObj; j++) {
                result[j][i] = data[j][i] - median;
            }
        }
        return result;
    }

    public Matrix covar(Matrix data) {
        return data.transpose().times(data).times(1.0 / data.getRowDimension());
    }

    public double[][] covarByMedian(double[][] data) {
        int numObj = data.length;
        int dim = data[0].length;
        double[] sortarray = new double[numObj];
        double[][] cov = new double[dim][dim];
        for (int l1 = 0; l1 < dim; l1++) {
            for (int l2 = 0; l2 < dim; l2++) {
                for (int j = 0; j < numObj; j++) {
                    sortarray[j] = data[j][l1] * data[j][l2];
                }

                Arrays.sort(sortarray);
                if (numObj % 2 == 0) {
                    cov[l1][l2] = (sortarray[numObj / 2 - 1] + sortarray[numObj / 2]) / 2;
                } else {
                    cov[l1][l2] = (sortarray[numObj / 2 - 1] + sortarray[numObj / 2] + sortarray[numObj / 2 + 1]) / 3;
                }
            }
        }
        double maxdiff = 0.0;
        for (int i = 0; i < dim; i++) {
            double sum = 0.0;
            for (int j = 0; j < dim; j++) {
                sum += Math.abs(cov[i][j]);
            }
            double diff = 1.1 * sum - 2.1 * cov[i][i];
            if (diff > maxdiff) {
                maxdiff = diff;
            }
        }
        if (maxdiff > 0.0) {
            for (int i = 0; i < dim; i++) {
                cov[i][i] += maxdiff;
            }
        }
        return cov;
    }

    public double[][] whitenData(double[][] data, double[][] transform, double[][] retransform) {
        int numObj = data.length;
        int dim = data[0].length;
        double[][] result = new double[numObj][dim];

        Matrix dm = new Matrix(data);
        // T E S T    T E S T    T E S T    T E S T    T E S T    T E S T
        Matrix c = covar(dm);
        //robust fit
        //Matrix c = new Matrix(covarByMedian(data)) ;
        EigenvalueDecomposition e = new EigenvalueDecomposition(c);
        Matrix evm = e.getV();
        //DEBUG
        double[] d = e.getRealEigenvalues();
        for (int i = 0; i < d.length; i++) {
            System.out.println(d[i]);
        }

        //DEBUG
        Matrix ewm = e.getD();
        if (verbose) {
            System.out.println("Whiten Data:");
            System.out.println("covariance matrix");
            // c.print(5,3) ;
            System.out.println("eigenvectors");
            // evm.print(5, 3);
            System.out.println("eigenvalues");
            //ewm.print(5, 3);
        }
        ev = ewm.getArrayCopy();
        double[] diag = new double[dim];
        for (int i = 0; i < dim; i++) {
            diag[i] = Math.sqrt(Math.min(Math.max(ewm.get(i, i), 1.0e-6), 1.0e+6));
            wscCost += Math.log(diag[i]) / Math.log(2.0);
        }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                transform[i][j] = evm.get(i, j) / diag[j];
                retransform[i][j] = evm.get(j, i) * diag[i];
            }
        }

        Matrix res = dm.times(new Matrix(transform));
        if (verbose) {
            System.out.println("Some sanity checks for data whitening:");
            System.out.println("transform &* retransform: This should be unity:");
            //(new Matrix(transform)).times(new Matrix(retransform)) .print(10,8) ;
        }
        return res.getArrayCopy();

    }

    public double g_func(double x) {
        return Math.tanh(x);
    }

    public double g_slope(double x) {
        return 1.0 - Math.pow(Math.tanh(x), 2);
    }

    public double curtosis(double[][] data, double[] w) {
        int numObj = data.length;
        int dim = data[0].length;
        double res = 0.0;
        double xsq = 0.0;
        for (int i = 0; i < numObj; i++) {
            double x = 0;
            for (int j = 0; j < dim; j++) {
                x += data[i][j] * w[j];
            }
            res += Math.pow(x, 4.0);
            xsq += x * x;
        }
        if (xsq / numObj > 1.1 || xsq / numObj < 0.9) {
            System.out.println("Error xsq=" + xsq / numObj);
        }
        return res / numObj - 3;
    }

//    private void rotateAndDisplay(double[][] data, double[][] w, int iter){
//        Matrix mw = new Matrix(w);
//        Matrix mwhitened = new Matrix(data);
//         double[][] projectedWhitened = mwhitened.times(mw.transpose()).getArrayCopy();
//         String title = new Integer(iter).toString();
//        displayData(projectedWhitened, title);
//        System.out.println("projectedWhitened:") ;
//            mwhitened.times(mw.transpose()).print(6,4) ;
//    }
    //data: Zeilenvektoren
    /**
     * Standard fastICA. Re-transformation to original space commented out.
     *
     * @param data data as row vetors.
     * @return w as row vectors.
     */
    public double[][] fastICA(double[][] data) {
        //verbose = true;
        // fastICA algorithm acc. to p. 194 with deflationary orthogonalization
        int numObj = data.length;
        int dim = data[0].length;
        Random r = new Random();

        // T E S T    T E S T    T E S T    T E S T    T E S T    T E S T
        //double[][] data1 = centerByMedian(data) ;
        double[][] data1 = centerData(data);
        //displayData(data1, "centered");
        double[][] transform = new double[dim][dim];
        double[][] retransform = new double[dim][dim];

        data1 = whitenData(data1, transform, retransform);
        //displayData(data1, "whitened");

        if (verbose) {
            Matrix dd = new Matrix(data1);
            System.out.println("Size of whitened data matrix: " + dd.getRowDimension() + " " + dd.getColumnDimension());
            System.out.println("For whitened data the covariance matrix must also be unity: ...");
            //covar(dd).print(10, 8) ;
            System.out.println("Daten:");
            //dd.print(6,4) ;
            System.out.println("-------");
        }

        double[][] w = new double[dim][dim];
        int alliter = 0;

        for (int n = 0; n < dim; n++) {

            for (int i = 0; i < dim; i++) {
                w[n][i] = r.nextDouble();
            }
            double sum = 0.0;
            for (int i = 0; i < dim; i++) {
                sum += w[n][i] * w[n][i];
            }
            for (int i = 0; i < dim; i++) {
                w[n][i] /= Math.sqrt(sum);
            }

            double[] old_w;
            double[] wUnorth;

            int iter = 0;
            do {
                if (verbose) {
                    System.out.println("w = [" + w[n][0] + " " + w[n][1] + "] curtosis = " + curtosis(data1, w[n]));
                }
                //if(iter == 1)
                //rotateAndDisplay(data1, w , iter);
                old_w = w[n].clone();
                double eLeft[] = new double[dim];
                double eRight = 0.0;
                for (int i = 0; i < numObj; i++) {
                    double scalar = 0.0;
                    for (int j = 0; j < dim; j++) {
                        scalar += w[n][j] * data1[i][j];
                    }
                    double g = g_func(scalar) / numObj;
                    for (int j = 0; j < dim; j++) {
                        eLeft[j] += g * data1[i][j];
                    }
                    eRight += g_slope(scalar) / numObj;

                }
                for (int j = 0; j < dim; j++) {
                    w[n][j] = eLeft[j] - eRight * w[n][j];
                }

                // deflationary orthogonalize w[n] w.r.t. vectors w[0]..w[n-1]

                wUnorth = w[n].clone();

                for (int m = 0; m < n; m++) {
                    double scalar = 0.0;
                    for (int j = 0; j < dim; j++) {
                        scalar += wUnorth[j] * w[m][j];
                    }
                    for (int j = 0; j < dim; j++) {
                        w[n][j] -= scalar * w[m][j];
                    }
                }
                sum = 0.0;
                for (int i = 0; i < dim; i++) {
                    sum += w[n][i] * w[n][i];
                }
                for (int i = 0; i < dim; i++) {
                    w[n][i] /= Math.sqrt(sum);
                }
                sum = 0.0;
                for (int i = 0; i < dim; i++) {
                    sum += w[n][i] * old_w[i];
                }
                iter++;
            } while (Math.abs(sum) < 0.9999 && iter < 100);
            if (verbose) {
                System.out.println("w = [" + w[n][0] + " " + w[n][1] + "] curtosis = " + curtosis(data1, w[n]));
                //if(iter == 1)
                //rotateAndDisplay(data1, w , iter);
                System.out.println("ITERATIONS: " + iter);
                Matrix ww = rowVector(w[n]).times(new Matrix(retransform));
                ww.print(6, 4);
//              (new Matrix(data)).print(6,4) ;
            }
            alliter += iter;
        }

        // retransform the matrix w to the original (not whitened) data space
        // w' = w' * retransform
        double[][] wn = new double[dim][dim];

        for (int m = 0; m < dim; m++) {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    wn[m][i] += w[m][j] * retransform[j][i];
                }
            }
        }

        // normalize the result (for the last time...)

        for (int m = 0; m < dim; m++) {
            double sum = 0.0;
            for (int i = 0; i < dim; i++) {
                sum += wn[m][i] * wn[m][i];
            }
            for (int i = 0; i < dim; i++) {
                wn[m][i] /= Math.sqrt(sum);
            }
        }


        if (verbose) {
            System.out.println("SUM OF ITERATIONS: " + alliter);
            (new Matrix(wn)).print(6, 4);
        }

        // To transform data such that independent components are located at
        // the coordinate axes, we have to multiply the data with wn^{-1} which
        // does not coincide with wn' (unlike the Eigenvector matrix of PCA)

        return wn;

//        if (verbose) {
//            System.out.println("SUM OF ITERATIONS: "+alliter) ;
//            (new Matrix(w)).print(6,4) ;
//        }
//        return w;

    }

//    private void displayData(double[][] d, String title){
//        int dim = 2;
//        DataObject[] data = new DataObject[d.length];
//        for(int i = 0; i < d.length; i++){
//            double[] coord = new double[dim];
//            for(int j = 0; j < dim; j++)
//                coord[j] = d[i][j];
//            data[i] = new DataObject(coord, i);
//        }
//        VisuInfo visu = new VisuInfo(data, title);
//        visu.setSize(600,600);
//        visu.setLocation(100,0);
//        visu.setVisible(true);
//    }
//    
    /**
     * Standard fastICA. Without whitening and re-transformation to original
     * space.
     *
     * @param data data[][] rows: objects, columns: dimensions
     * @return Matrix W: independent components as row vectors
     */
    public double[][] optimizeEntropy(double[][] data) {
        // fastICA algorithm acc. to p. 194 with deflationary orthogonalization
        int numObj = data.length;
        int dim = data[0].length;
        Random r = new Random();
        double[][] w = new double[dim][dim];
        int alliter = 0;
        for (int n = 0; n < dim; n++) {
            for (int i = 0; i < dim; i++) {
                w[n][i] = r.nextDouble();
            }
            double sum = 0.0;
            for (int i = 0; i < dim; i++) {
                sum += w[n][i] * w[n][i];
            }
            for (int i = 0; i < dim; i++) {
                w[n][i] /= Math.sqrt(sum);
            }
            double[] old_w;
            double[] wUnorth;
            int iter = 0;
            do {
                if (verbose) {
                    System.out.println("w = [" + w[n][0] + " " + w[n][1] + "] curtosis = " + curtosis(data, w[n]));
                }
                old_w = w[n].clone();
                double eLeft[] = new double[dim];
                double eRight = 0.0;
                for (int i = 0; i < numObj; i++) {
                    double scalar = 0.0;
                    for (int j = 0; j < dim; j++) {
                        scalar += w[n][j] * data[i][j];
                    }
                    double g = g_func(scalar) / numObj;
                    for (int j = 0; j < dim; j++) {
                        eLeft[j] += g * data[i][j];
                    }
                    eRight += g_slope(scalar) / numObj;
                }
                for (int j = 0; j < dim; j++) {
                    w[n][j] = eLeft[j] - eRight * w[n][j];
                }
                // deflationary orthogonalize w[n] w.r.t. vectors w[0]..w[n-1]
                wUnorth = w[n].clone();
                for (int m = 0; m < n; m++) {
                    double scalar = 0.0;
                    for (int j = 0; j < dim; j++) {
                        scalar += wUnorth[j] * w[m][j];
                    }
                    for (int j = 0; j < dim; j++) {
                        w[n][j] -= scalar * w[m][j];
                    }
                }
                sum = 0.0;
                for (int i = 0; i < dim; i++) {
                    sum += w[n][i] * w[n][i];
                }
                for (int i = 0; i < dim; i++) {
                    w[n][i] /= Math.sqrt(sum);
                }
                sum = 0.0;
                for (int i = 0; i < dim; i++) {
                    sum += w[n][i] * old_w[i];
                }
                iter++;
            } while (Math.abs(sum) < 0.99999999 && iter < 100); // <0.9999 && iter<100
            if (verbose) {
                System.out.println("w = [" + w[n][0] + " " + w[n][1] + "] curtosis = " + curtosis(data, w[n]));
                System.out.println("ITERATIONS: " + iter);
            }
            alliter += iter;
        }
        if (verbose) {
            System.out.println("SUM OF ITERATIONS: " + alliter);
            (new Matrix(w)).print(6, 4);
        }

        // To transform data such that independent components are located at
        // the coordinate axes, we have to multiply the data with w^{-1} which
        // does not coincide with wn' (unlike the Eigenvector matrix of PCA)

        //(new Matrix (data1)).times((new Matrix (w)).transpose()).print (6,4) ;
        return w;
    }

    private double[][] covar(double[][] d) {
        int rows = d.length;
        int cols = d[0].length;
        double[][] cov = new double[rows][cols];
        double[] colmeans = new double[cols];
        //remove the mean of each column
        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < rows; j++) {
                colmeans[i] += d[j][i] / rows;
            }
        }
        for (int j = 0; j < cols; j++) {
            for (int k = 0; k < rows; k++) {
                d[k][j] = d[k][j] - colmeans[j];
            }
        }
        Matrix dm = new Matrix(d);
        Matrix cm = dm.transpose().times(dm);
        int size = cm.getRowDimension();
        cov = cm.getArrayCopy();
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                cov[i][j] /= (double) rows;
            }
        }
        return cov;

    }

    private void dataCheck(double[][] d, String filename) {
        MLDouble q = new MLDouble("test", d);
        ArrayList ll = new ArrayList();
        ll.add(q);
        MatFileWriter mw = new MatFileWriter();
        try {
            String name = filename + ".mat";
            mw.write(name, ll);
        } catch (IOException ex) {
            //  Logger.getLogger(VI.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    //centers and performs eigenvalue decomposition. Writes centered, E and D
    private void pca(double[][] d) {
        if (!minimalOutput) {
            System.out.println("---pca---");
        }
//       //DEBUG: data test
        if (!minimalOutput) {
            dataCheck(d, "start");
        }
//        //DEBUG
        //center data
        //d = centerData(d);
        d = rowCenterData(d);
        centered = d;
        //DEBUG
        if (!minimalOutput) {
            dataCheck(d, "centered");
        }
        //DEBUG
        //compute covariance Matrix of data
//        double[][] cov = covar(new Matrix(d).transpose().getArrayCopy());
//        //DEBUG
//        dataCheck(cov, "covar");
//        //DEBUG
        Matrix dd = new Matrix(d);//transpose??
        Matrix covar = new Matrix(covar(dd.transpose().getArrayCopy()));
        EigenvalueDecomposition e = new EigenvalueDecomposition(covar);
        E = e.getV();
        D = e.getD();
    }

    private double[][] rowCenterData(double[][] data) {
        int numObj = data.length;
        int d = data[0].length;
        double[][] result = new double[numObj][d];
        mean = new double[numObj];
        for (int i = 0; i < numObj;
                i++) {
            double sum = 0.0;
            for (int j = 0; j < d;
                    j++) {
                sum += data[i][j];
            }
            sum /=
                    d;
            mean[i] = sum;
            for (int j = 0; j < d;
                    j++) {
                result[i][j] = data[i][j] - sum;
            }
        }

        return result;
    }

    //reduce dimensionality to numDimensions
    private void reduceDimensionality(int numDimensions) {
        //select numdimensions colums of ev and ew
        int remove = D.getRowDimension() - numDimensions;
        double[] discardedEigenvalues = new double[remove];
        double[] eigenvalues = new double[D.getRowDimension()];
        for (int i = 0; i < D.getRowDimension(); i++) {
            eigenvalues[i] = D.get(i, i);
        }
        double[][] redEW = new double[numDimensions][numDimensions];
        double[][] ww = D.getArrayCopy();
        int counter = 0;
        for (int i = 0; i < remove; i++) {
            discardedEigenvalues[i] = ww[i][i];
        }
        for (int i = remove; i < D.getRowDimension(); i++) {
            redEW[counter][counter] = ww[i][i];
            counter++;
        }
        double[][] redEV = new double[E.getRowDimension()][numDimensions];
        double[][] vv = E.getArrayCopy();
        for (int i = 0; i < E.getRowDimension(); i++) {
            counter = 0;
            for (int j = remove; j < E.getColumnDimension(); j++) {
                redEV[i][counter] = vv[i][j];
                counter++;
            }

        }

//        //DEBUG
//        double[] v = e.getRealEigenvalues();
//        for (int i = 0; i < v.length; i++) {
//            System.out.println(v[i]);
//        }
//    //DEBUG
        reducedE = new Matrix(redEV);
        reducedD = new Matrix(redEW);
        // double bic = mdlF(discardedEigenvalues, D.getRowDimension(), numDimensions);
//        //double bic = mdlNew(eigenvalues, D.getRowDimension(), numDimensions);
        //System.out.println("DIM: " + numDimensions + " MDL: " + bic);

    }
}
