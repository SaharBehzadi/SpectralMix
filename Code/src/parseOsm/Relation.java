/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package parseOsm;

import java.util.List;
import java.util.Map;

/**
 *
 * @author Willy Tiengo
 */
public class Relation extends AbstractNode {

    public List<Member> members;
    public Map<String, String> tags;

    public Relation(String id, String visible, String timestamp,
            String version, String changeset, String user,
            String uid, List<Member> members, Map<String, String> tags) {
       
        super(id, visible, timestamp, version, changeset, user, uid);
        this.members = members;
        this.tags = tags;
    }
}

