/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

/**
 *
 * @author claudia.plant
 */
public class MetricCorrection {

    double[] distanceMatrix;
    double[] metric;
    int numObj;
    static double maxVio = 0.0000000001;
    boolean verbose = true;
    int iter;

    public MetricCorrection(double[] distanceMatrix, int numObj) {
        this.distanceMatrix = distanceMatrix;
        metric = new double[distanceMatrix.length];
        for(int i = 0; i < metric.length; i++)
            metric[i] = distanceMatrix[i];
        this.numObj = numObj;
        iter = 0;
    }
    
    public void run(){
        boolean metric = false;
        while(! metric){
            metric = metricCorrectionOneStep();
            if (verbose)
                    System.out.println(iter);
            iter ++;
        }
        
    }
    
    
    public boolean metricCorrectionOneStep() {
        DataUtils du = new DataUtils();
        boolean passed = true;
     
      
            for (int i = 0; i <= numObj - 3; i++) {
                for (int j = i + 1; j <= numObj - 2; j++) {
                    for (int kk = j + 1; kk <= numObj - 1; kk++) {
                        // int ijk = du.getIndex3(i, j, kk, numObj);
                        int ki = du.getIndex(kk, i, numObj);
                        int jk = du.getIndex(j, kk, numObj);
                        int ij = du.getIndex(i, j, numObj);

                        double diff = metric[ki] + metric[jk] - metric[ij];
                         if (diff < -maxVio) {
                        //if (diff < -minIntervallSize) {
                            diff = -diff / 3.0;
                            metric[ij] -= diff;
                            metric[ki] += diff;
                            metric[jk] += diff;
                            passed = false;
                        }
                        diff = metric[ki] + metric[ij] - metric[jk];
                        if (diff < -maxVio) {
                            diff = -diff / 3.0;
                            metric[ij] += diff;
                            metric[ki] += diff;
                            metric[jk] -= diff;
                            passed = false;
                        }
                        diff = metric[ij] + metric[jk] - metric[ki];
                        if (diff < -maxVio) {
                            diff = -diff / 3.0;
                            metric[ij] += diff;
                            metric[ki] -= diff;
                            metric[jk] += diff;
                            passed = false;
                        }
                        //System.out.println("m");

                    }
                }
            }
        
        return passed;

    }
    
    
}
