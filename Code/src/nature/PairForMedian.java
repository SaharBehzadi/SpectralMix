/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nature;

/**
 *
 * @author boehm
 */
public class PairForMedian implements Comparable {

    public double value;
    public double weight;

    public PairForMedian(double value, double weight) {
        this.value = value;
        this.weight = weight;
    }

    public int compareTo(Object o) {
        PairForMedian other = (PairForMedian) o;
        if (this.value < other.value) {
            return -1;
        }
        if (this.value > other.value) {
            return 1;
        }
        return 0;
    }
}

