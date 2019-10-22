/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package attributedEmbedding;

import datautils.Sociopatterns;
import edu.uci.ics.jung.graph.Graph;

/**
 *
 * @author plantc59cs
 */
public class MixedSpectralMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
         DataGenerator dg = new DataGenerator();
       dg.agreementBetweenGraphAndAttributes();
        
        
        Sociopatterns sp = new Sociopatterns();
        sp.readData();
        Graph[] measuredDiary = new Graph[1];
        //measuredDiary[0] = sp.getMeasured();
        //measuredDiary[0] = dg.g;
        measuredDiary[0] = sp.getDiary();
        
        int[][] attributes = new int[sp.getN()][2];
        //int[][] emptyAtt = new int[sp.getN()][1];
        
         //int[][] attributes = new int[dg.g.getVertexCount()][1];
        int[] sex = sp.getGender();
        int[] cl = sp.getSchoolclass();
        
        for (int i = 0; i < attributes.length; i++) {
            attributes[i][0] = sex[i];
            //attributes[i][0] = dg.attribute[i];
            attributes[i][1] = cl[i];
        }

        //public MixedSpectral(Graph[] g, int[][] attributes, int numAtt, int d, int[] classId) {
        boolean[] weighted = new boolean[1];
        weighted[0] = true;
       // weighted[1] = true;
        //public MixedSpectral(Graph[] g, boolean[] weighted, int[][] attributes, int numAtt, int d, int[] classId) {
        MixedSpectral ms = new MixedSpectral(measuredDiary, weighted, attributes, 2, 3, sp.getClassId());
        ms.init(0);
        ms.run();
    }
    
}
