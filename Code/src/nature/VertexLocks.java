/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import java.util.concurrent.locks.ReentrantLock;

/**
 *
 * @author plantc59cs
 */
public class VertexLocks {
    int num;
    ReentrantLock lock;

    public VertexLocks(int num) {
        this.num = num;
        lock = new ReentrantLock(true);
    }
    
}
