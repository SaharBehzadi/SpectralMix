/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package parseOsm;

import java.util.Map;

/**
 *
 * @author Willy Tiengo
 */
public class OSMNode extends AbstractNode {

    public String lat;
    public String lon;
    public Map<String, String> tags;

    public OSMNode(String id, String visible, String timestamp, String version, String changeset, String user, String uid, String lat, String lon, Map<String, String> tags) {
        super(id, visible, timestamp, version, changeset, user, uid);
        this.lat = lat;
        this.lon = lon;
        this.tags = tags;
    }

}

