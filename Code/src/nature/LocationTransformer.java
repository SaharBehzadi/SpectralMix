/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;


import java.awt.Dimension;
import java.awt.geom.Point2D;
import java.awt.geom.Point2D.Double;
import org.apache.commons.collections15.Transformer;

/**
 *
 * @author claudia.plant
 */
public class LocationTransformer implements Transformer {
    public double[][] locations;
    public Dimension size;
    public int offset = 10;

    public void setLocations(double[][] locations) {
        this.locations = locations;
        
    }

    public void setSize(Dimension size) {
        this.size = size;
    }
    

    

    public Object transform(Object i) {
      return new Point2D.Double(locations[0][((Integer)i).intValue()] * size.width + offset, (1- locations[1][((Integer)i).intValue()]) * size.height + offset);
    }

   

    
    
}
