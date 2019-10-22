/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;


import java.util.*;

/**
 *
 * @author boehm
 */
public class SortPack implements Comparable {
    
    public double value ;
    public int index ;
    
    /** Creates a new instance of sortPack */
    public SortPack() {
    }
    
    public SortPack(double v, int i) {
        value=v ;
        index=i ;
    }
    
    public int compareTo (Object obj) {
        SortPack o=(SortPack)obj ;
        if (value<o.value) return -1 ;
        if (value>o.value) return +1 ;
        return 0 ;
    }
    
   public static SortPack [] array (double [] arr) {
        SortPack [] result = new SortPack[arr.length] ;
        for (int i=0 ; i<arr.length ; i++)
            result[i] = new SortPack(arr[i],i) ;
        Arrays.sort(result) ;
        return result ;
    }
}
