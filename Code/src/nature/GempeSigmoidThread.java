/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

/**
 *
 * @author plantc59cs
 */
public class GempeSigmoidThread implements Runnable {   
    GempePar p;
//    int n;
//    int d;
    boolean verbose = true;

    public GempeSigmoidThread(GempePar p) {
        this.p = p;
    }

    public void run() {
        if (verbose) {
            System.out.println("Sigmoid calculation thread started.");
        }
        p.computeSigmoid();
        
//        p.computeSigmoid();
//        Visualization v = new Visualization(p.g);
//        v.displayCoordNew(p.coord, Integer.toString(seed));

    }

  
}
