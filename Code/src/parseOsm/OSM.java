/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package parseOsm;

import java.util.Map;
import java.util.Set;

/**
 *
 * @author Willy Tiengo
 */
public class OSM {

    private Set<OSMNode> nodes;
    private Set<Way> ways;
    private Set<Relation> relations;

    public OSM(Set<OSMNode> nodes, Set<Way> ways,
            Set<Relation> relations) {
        this.nodes = nodes;
        this.ways = ways;
        this.relations = relations;
    }

    public OSM() {
    }

  
    public Set<OSMNode> getNodes() {
        return nodes;
    }

    public Set<Relation> getRelations() {
        return relations;
    }

    public Set<Way> getWays() {
        return ways;
    }
}

