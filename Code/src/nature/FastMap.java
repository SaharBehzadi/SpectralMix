


package nature;

import javax.vecmath.GMatrix;
import javax.vecmath.GVector;
import java.util.Random;
import java.util.Vector;

/**
 * FastMap clustering algorithm (Java implementation)
 *
 * [FastMap: A Fast Algorithm for Indexing, Data-Mining and Visualization of
 *  Traditional and Multimedia Datasets;
 *  Christos Faloutsos, King-Ip (David) Lin;
 *  University of Maryland]
 *
 * @author    Kai Jauslin
 * @version   1.0
 */

public class FastMap {
    /* FastMap cached pivot objects */
    private int[][] pivots;
    private GVector[] result;
    private GMatrix distances;
    private int targetDimension;
    private int numberOfObjects;
    private int dim = 0;
    private double pivotDistance;
    
    
    public FastMap(GMatrix distances, int targetDimension) {
        this.distances = distances;
        if (distances.getNumCol() != distances.getNumRow()) {
            // error: distance matrix is not square
        }
        
        this.targetDimension = targetDimension;
        this.numberOfObjects = distances.getNumCol();
        
        // initialize result (can maybe use object cloning?)
        this.result = new GVector[numberOfObjects];
        for (int i=0; i<numberOfObjects; i++)
            result[i] = new GVector(targetDimension);
        
        // initialize pivot storage
        pivots = new int[targetDimension][2];
        
        // initialize random number generator
        Math.random();
    }
    
    public GVector[] getResult(){
        return this.result;
    }
    
    public int getDim(){
        return this.dim;
    }
    
    public int getNumObj(){
        return this.numberOfObjects;
    }
    
    private double getDistance(int p, int q) {
        // only lower half of distance square matrix is filled
        if (p == q) return 0.0;
        
        // get base distance in n-dimensional space
        double basedist;
        if (p > q) {
            basedist = distances.getElement(q, p);
        } else {
            basedist = distances.getElement(p, q);
        }
        
        // take into account dimensionality reduction sum(x_p - x_q)
        if (dim > 0) {
            double diff = 0.0;
            double sqdiff = 0.0;
            for (int i=0; i<dim; i++) {
                diff = result[p].getElement(i) - result[q].getElement(i);
                sqdiff += diff*diff;
            }
            
            double basedistsq = basedist*basedist;
            if (basedistsq - sqdiff < 0) {
                return Math.sqrt(sqdiff - basedistsq);
            } else {
                return Math.sqrt(basedistsq - sqdiff);
            }
        } else {
            return basedist;
        }
    }
    
    /**
     * Finds object with maximum distance to object o.
     * @dist distance matrix
     * @o source object
     */
    private int findMaxDistance(int o) {
        double maxdist = 0.0;
        double value;
        int targetObject = o;
        
        for (int j=0; j<numberOfObjects; j++) {
            value = getDistance(j, o);
            if (value > maxdist) {
                maxdist = value;
                targetObject = j;
            }
        }
        
        pivotDistance = maxdist;
        return targetObject;
    }
    
    /**
     * Heuristic to choose most distant objects.
     * @dist distance matrix of objects
     * @returns index of two pivot objects
     */
    private void choosePivots() {
        // choose random object O_b
        double maxDistance = 0.0;
        
        // number of trials to find pivot objects
        int numberOfTrials = 5;
        
        int O_a, O_b;
        for (int i=0; i<numberOfTrials; i++) {
            O_b = (int)(Math.random()*numberOfObjects);
            // "they" always start with the (second last) element?
            //O_b = numberOfObjects-2;
            if (O_b < 0) O_b = 0;
            O_a = findMaxDistance(O_b);
            O_b = findMaxDistance(O_a);
            
            // get distance
            if (pivotDistance > maxDistance) {
                maxDistance = pivotDistance;
                pivots[dim][0] = O_a;
                pivots[dim][1] = O_b;
            }
        }
        
        pivotDistance = maxDistance;
    }
    
    /**
     * Performs one iteration step of the fastmap mapping.
     * @dist distance matrix needed for this dimension
     * @dim dimension this step is mapping to
     */
    private void mapStep() {
        
        // choose two distant objects as pivots
        choosePivots();
        int pa = pivots[dim][0];
        int pb = pivots[dim][1];
        
        // project objects to line
        double x, dai, dbi;
        for (int i=0; i<numberOfObjects; i++) {
            dai = getDistance(pa, i);
            dbi = getDistance(pb, i);
            x = (dai*dai + pivotDistance*pivotDistance - dbi*dbi)/(2*pivotDistance);
            result[i].setElement(dim, x);
        }
    }
    
    public GVector[] map() {
        for (dim=0; dim<targetDimension; dim++) {
            //System.out.println("Step "+dim);
            mapStep();
            // write out object positions and pivots
            //System.out.println("Pivots = ("+pivots[dim][0]+", "+pivots[dim][1]+")");
            //for (int i=0; i<numberOfObjects; i++) {
            //  System.out.println("Object "+i+" = "+result[i].getElement(dim));
            //}
        }
        return result;
    }
    
   /* public static void main(String[] args) {
        GMatrix test = new GMatrix(5,5);
        double[] r1 = {0,1,1,100,100};
        test.setColumn(0, r1);
        double[] r2 = {1,0,1,100,100};
        test.setColumn(1, r2);
        double[] r3 = {1,1,0,100,100};
        test.setColumn(2, r3);
        double[] r4 = {100,100,100,0,1};
        test.setColumn(3, r4);
        double[] r5 = {100,100,100,1,0};
        test.setColumn(4, r5);
    
        FastMap f = new FastMap(test, 2);
        f.map();
        System.out.println("numberOfObjects " + f.getNumObj());
        System.out.println("dim " + f.getDim());
        for (int i=0; i< f.getNumObj(); i++) {
             System.out.println("Object "+i+" Koordinaten: ");
             for(int j=0; j<f.getDim(); j++){
                 System.out.print(f.getResult()[i].getElement(j) + " ");
             }
             System.out.println();
        }
    }*/
    
    
}
