/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

/**
 *
 * @author claudia.plant
 */
public class ClusterPath {
int pointId;
int reachedFrom;

    public ClusterPath(int pointId, int reachedFrom) {
        this.pointId = pointId;
        this.reachedFrom = reachedFrom;
    }

    public int getPointId() {
        return pointId;
    }

    public int getReachedFrom() {
        return reachedFrom;
    }
    

}
