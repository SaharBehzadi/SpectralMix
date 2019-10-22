/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package datautils;

import edu.uci.ics.jung.graph.DirectedSparseMultigraph;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseMultigraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import java.util.HashMap;
import java.util.Vector;
import nature.IO;

/**
 *
 * @author plantc59cs
 */
public class Sociopatterns {

    int[] gender;
    int[] schoolclass;
    int[] classId;
    int[][] facebookLinks; //row: one student, columns: 1: link, 0: not link, -1 dont know
    int[] nodeNameFacebookRow;
    Graph measured;
    Graph diary;
    Graph friendship;
    int[] nodeName;
    HashMap<Integer, Integer> idToNodeName;
    HashMap<Integer, Integer> nodeNameToId;
    MyNode[] nodes;
    int contactsSumOfWeights;
    int genderSumOfWeights; //without -1
    static int n = 329;

    public Sociopatterns() {

    }

    public void readData() {
        initNameMap();
        readFriendship();
        readDiary();
        readContacts();
        readFacebook();
        System.out.println("data read");
    }

    public Graph getFriendship() {
        return friendship;
    }

    public static int getN() {
        return n;
    }
    
    

    public int[] getGender() {
        return gender;
    }

    public int[] getSchoolclass() {
        return schoolclass;
    }

    public int[] getClassId() {
        return classId;
    }

    public int[][] getFacebookLinks() {
        return facebookLinks;
    }

    public Graph getMeasured() {
        return measured;
    }

    public Graph getDiary() {
        return diary;
    }

    public int getContactsSumOfWeights() {
        return contactsSumOfWeights;
    }

    public int getGenderSumOfWeights() {
        return genderSumOfWeights;
    }
    
    

    public void initNameMap() {
        IO ea = new IO();
        double[][] dd = ea.readMatlabMatrix("idToNodename.mat", "idToNodename");
        double[][] cl = ea.readMatlabMatrix("classlabel.mat", "classlabel");
        double[][] gg = ea.readMatlabMatrix("sex.mat", "sex");
        double[][] sc = ea.readMatlabMatrix("schoolclass.mat", "schoolclass");
        idToNodeName = new HashMap<Integer, Integer>();
        nodeNameToId = new HashMap<Integer, Integer>();
        nodes = new MyNode[n];
        classId = new int[n];
        gender = new int[n];
        classId = new int[n];
        schoolclass = new int[n];
        for (int i = 0; i < dd.length; i++) {
            idToNodeName.put((int) dd[i][0], i);
            nodeNameToId.put(i, (int) dd[i][0]);
            nodes[i] = new MyNode(i);
            classId[i] = (int) cl[i][0];
            gender[i] = (int) gg[i][0];
            if(gender[i] != -1)
                genderSumOfWeights++;
            schoolclass[i] = (int) sc[i][0];
        }

    }

    public void readDiary() {
        diary = new DirectedSparseMultigraph<Integer, MyEdge>();
        for (int i = 0; i < nodes.length; i++) {
            diary.addVertex(i);
        }
        IO ea = new IO();
        double[][] dd = ea.readMatlabMatrix("diaryNetwork.mat", "diaryNetwork");
        int edgeCounter = 0;
        for (int i = 0; i < dd.length; i++) {
            diary.addEdge(new MyEdge(edgeCounter, dd[i][2]), idToNodeName.get((int) dd[i][0]), idToNodeName.get((int) dd[i][1]), EdgeType.DIRECTED);
            edgeCounter++;
        }
    }

    public void readContacts() {
        IO ea = new IO();
        measured = new UndirectedSparseMultigraph<Integer, MyEdge>();
        double[][] dd = ea.readMatlabMatrix("measuredContacts.mat", "measuredContacts");
        int[][] numContacts = new int[n][n];
        for (int i = 0; i < dd.length; i++) {
            numContacts[idToNodeName.get((int) dd[i][0])][idToNodeName.get((int) dd[i][1])]++;
            numContacts[idToNodeName.get((int) dd[i][1])][idToNodeName.get((int) dd[i][0])]++;
        }
        for (int i = 0; i < nodes.length; i++) {
            measured.addVertex(i);
        }
        int edgeCounter = 0;
        for (int i = 0; i < n; i++) {
            for (int j = (i + 1); j < n; j++) {
                if (numContacts[i][j] > 0) {
                    measured.addEdge(new MyEdge(edgeCounter, (double) numContacts[i][j]), i, j);
                    edgeCounter++;
                    contactsSumOfWeights += numContacts[i][j];
                    if(measured.getVertexCount() > 329)
                        System.out.println("m");
                }
            }
        }
        System.out.println("m");
    }

    public void readFacebook() {
        IO ea = new IO();
        double[][] dd = ea.readMatlabMatrix("facebookNetwork.mat", "facebookNetwork");
        Vector<Integer> firstNodes = new Vector<Integer>();
        for (int i = 0; i < dd.length; i++) {
            if (!firstNodes.contains((int) dd[i][0])) {
                firstNodes.add((int) dd[i][0]);
            }
        }
        System.out.println(firstNodes.size());
        facebookLinks = new int[firstNodes.size()][n];
        for (int i = 0; i < facebookLinks.length; i++) {
            for (int j = 0; j < n; j++) {
                facebookLinks[i][j] = -1;
            }
        }
        nodeNameFacebookRow = new int[firstNodes.size()];
        HashMap<Integer, Integer> idToIndex = new HashMap<Integer, Integer>();
        for (int i = 0; i < firstNodes.size(); i++) {
            nodeNameFacebookRow[i] = idToNodeName.get(firstNodes.elementAt(i).intValue());
            idToIndex.put(firstNodes.elementAt(i).intValue(), i);
        }
        for (int i = 0; i < dd.length; i++) {
            facebookLinks[idToIndex.get((int) dd[i][0])][idToNodeName.get((int) dd[i][1])] = (int) dd[i][2];
        }

    }

    public void readFriendship() {
        friendship = new DirectedSparseMultigraph<MyNode, MyEdge>();
        for (int i = 0; i < nodes.length; i++) {
            friendship.addVertex(nodes[i]);
        }
        IO ea = new IO();
        double[][] dd = ea.readMatlabMatrix("friendshipNetwork.mat", "friendshipNetwork");
        int edgeCounter = 0;
        for (int i = 0; i < dd.length; i++) {
            friendship.addEdge(new MyEdge(edgeCounter, 1.0), nodes[idToNodeName.get((int) dd[i][0])], nodes[idToNodeName.get((int) dd[i][1])]);
            edgeCounter++;
        }
    }

}
