/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package parseOsm;

import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author claudia.plant
 */
public class NewMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String dir = "export.json";
        JSONParser jp = new JSONParser(dir);
        
        
//        OSM test = new OSM();
//        try {
//           test = OSMParser.parse(dir);
//        } catch (Exception ex) {
//            Logger.getLogger(NewMain.class.getName()).log(Level.SEVERE, null, ex);
//        }
//        Object [] ways = test.getRelations().toArray();
//        for(int i = 0; i < ways.length; i++){
//            Relation nn = (Relation) ways[i];
//            if (nn.tags.containsKey("route"))
//               System.out.println("m");
//        }
        
      //  System.out.println("m");
    }
}
