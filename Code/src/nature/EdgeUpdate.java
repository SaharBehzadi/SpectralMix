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
public class EdgeUpdate {
    int endpointI;
    int endpointJ;
    int pNotEdgeI; //not edge considered in previous update
    int pNotEdgeJ;

    public EdgeUpdate(int endpointI, int endpointJ, int pNotEdgeI, int pNotEdgeJ) {
        this.endpointI = endpointI;
        this.endpointJ = endpointJ;
        this.pNotEdgeI = pNotEdgeI;
        this.pNotEdgeJ = pNotEdgeJ;
    }
    
    
    
}
