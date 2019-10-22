/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;

/**
 *
 * @author plantc59cs
 */
public class EvaluationLinkPredictionMain {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
       EvaluationLinkPrediction evl = new EvaluationLinkPrediction();
        evl.stackoverflowTrainGraph();
        evl.edgesToPredict();
    }
    
}
