/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import java.util.HashSet;

/**
 *
 * @author plantc59cs
 */
public class Cluster {
    HashSet<Integer> members;
    double[] sum;

    public Cluster(int d) {
        sum = new double[d];
        members = new HashSet<Integer>();
    }
    
}
