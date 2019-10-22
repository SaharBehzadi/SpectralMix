/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.EdgeType;

/**
 *
 * @author claudia
 */
public class GraphMDL {

    Graph g;
    int numVertices;
    int numEdges;
    int numPossEdges; //possible number of edges

    public GraphMDL(Graph g) {
        this.g = g;
        numVertices = g.getVertexCount();
        numEdges = g.getEdgeCount() - numVertices;
        numPossEdges = (numVertices * (numVertices - 1)) / 2;
    }

    public void printStatistics() {
        System.out.println("numNodes: " + numVertices + " numEdges: " + numEdges + " percentage set: " + (double) numEdges / (double) numPossEdges + " coding cost per entry: " + Entropy() + " total coding cost: " + Entropy() * numPossEdges);
    }

    public double Entropy() {
        double pEdge = (double) numEdges / (double) numPossEdges;
        double pNoEdge = 1.0 - pEdge;
        double entropy = -(pEdge * lg2(pEdge) + pNoEdge * lg2(pNoEdge));
        return entropy;
        //return entropy * numPossEdges;
    }

    public int[] countTripletFrequencies() {
        int[] counts = new int[8];
        int countCombis = 0;
        for (int i = 0; i < numVertices; i++) {
            for (int j = i + 1; j < numVertices; j++) {
                for (int k = j + 1; k < numVertices; k++) {
                    if (i != j && i != k && j != k) {
                        if (!g.isNeighbor(i, j) && !g.isNeighbor(i, k) & !g.isNeighbor(j, k)) //case a
                        {
                            counts[0]++;
                        }
                        if (g.isNeighbor(i, j) && !g.isNeighbor(i, k) & !g.isNeighbor(j, k))//case b
                        {
                            counts[1]++;
                        }
                        if (!g.isNeighbor(i, j) && g.isNeighbor(i, k) & !g.isNeighbor(j, k))//case c
                        {
                            counts[2]++;
                        }
                        if (!g.isNeighbor(i, j) && !g.isNeighbor(i, k) & g.isNeighbor(j, k))//case d
                        {
                            counts[3]++;
                        }
                        if (g.isNeighbor(i, j) && g.isNeighbor(i, k) & !g.isNeighbor(j, k))//case e
                        {
                            counts[4]++;
                        }
                        if (g.isNeighbor(i, j) && !g.isNeighbor(i, k) & g.isNeighbor(j, k))//case f
                        {
                            counts[5]++;
                        }
                        if (!g.isNeighbor(i, j) && g.isNeighbor(i, k) & g.isNeighbor(j, k))//case g
                        {
                            counts[6]++;
                        }
                        if (g.isNeighbor(i, j) && g.isNeighbor(i, k) & g.isNeighbor(j, k))//case h -- triangle
                        {
                            counts[7]++;
                        }
                        countCombis++;
                    }

                }
            }
        }
        System.out.println("countCombis: " + countCombis);
        System.out.println("case A: " + counts[0] + " case B: " + counts[1] + " case C: " + counts[2] + " case D: " + counts[3] + " case E: " + counts[4] + " case F: " + counts[5] + " case G: " + counts[6] + " case H: " + counts[7]);
        double[] percentages = new double[8];
        for (int i = 0; i < percentages.length; i++) {
            percentages[i] = (double) counts[i] / countCombis;
        }
        System.out.println(percentages[0] + "  " + percentages[1] + " " + percentages[2] + "  " + percentages[3] + "  " + percentages[4] + "  " + percentages[5] + " " + percentages[6] + " " + percentages[7]);
        return counts;
    }

    public double codingWithTrianglesNew(double probEdge, boolean verbose) {
        double[][] codingCosts = new double[numVertices][numVertices];
        double probNoEdge = 1.0 - probEdge;
        double ent = -(probEdge * lg2(probEdge) + probNoEdge * lg2(probNoEdge));
        double firstDiag = (numVertices - 1) * ent; //costs for the first diagonal
        int[] num = countTripletFrequencies();
        int t3 = num[7]; //number of 3-groups with 3 edges
        int t2 = num[6] + num[5] + num[4];        //number of 3-groups with 2 edges
        int t1 = num[3] + num[2] + num[1]; //number of 3-groups with 1 edge
        int t0 = num[0]; //number of 3-groups without edges
        double uf = t1 / 3.0;
        double lf = t0 + uf;
        double pFirstEdge = uf / lf;
        double us = 2.0 / 3.0 * t2;
        double ls = 2.0 / 3.0 * t1 + us;
        double pSecondEdge = us / ls;
        double pClosingTriangle = (double) t3 / (1.0 / 3.0 * t2 + t3);
        double costFirstEdge = -lg2(pFirstEdge);
        double costNotFirstEdge = -lg2(1 - pFirstEdge);
        double costSecondEdge = -lg2(pSecondEdge);
        double costNotSecondEdge = -lg2(1 - pSecondEdge);
        double costClosingTriangle = -lg2(pClosingTriangle);
        double costNotClosingTriangle = -lg2(1 - pClosingTriangle);

//        //DEBUG - ok!
//        double costFirstEdge = -lg2(probEdge);
//        double costNotFirstEdge = -lg2(probEdge);
//        double costSecondEdge = -lg2(probEdge);
//        double costNotSecondEdge = -lg2(probEdge);
//        double costClosingTriangle = -lg2(probEdge);
//        double costNotClosingTriangle = -lg2(probEdge);


        if (verbose) {
            System.out.println("Cost for first edge: " + costFirstEdge + " not existing: " + costNotFirstEdge);
            System.out.println("Cost for second edge: " + costSecondEdge + " not existing: " + costNotSecondEdge);
            System.out.println("Cost for closing triangle: " + costClosingTriangle + " not closing: " + costNotClosingTriangle);

        }


        double cc = 0.0;
        int col = 2;
        int row = 0;
        //DEBUG
        // int numEntries = numVertices-1; - ok!
       // System.out.println("Diagonal: " + col);
        while (col < numVertices) {
            int currCol = col; //start with current column
            row = 0;
            while (currCol < numVertices) {
                boolean triangle = false;  //check if this entry would be the second edge or close a triangle among the stuff we have seen before
                boolean secondEdge = false;
                boolean firstEdge = false;
                int colT = currCol - 1;
                int rowT = currCol - 1;
               System.out.println("entry: " + row + " " + currCol + " " + g.isNeighbor(currCol, row) + " -----------");
                int countTriangle = 0;
                int nodeCounter = 0;
                while (colT > row && rowT > row) {
                    if (g.isNeighbor(colT, row) || g.isNeighbor(rowT, currCol)) {
                        if (g.isNeighbor(colT, row) && g.isNeighbor(rowT, currCol)) {
                            triangle = true;
                            System.out.println(" would form triangle with: " + rowT);
                            countTriangle++;
                        } else {
                            secondEdge = true;
                        }
                    } else {
                        firstEdge = true;
                    }
                    colT--;
                    rowT--;
                    nodeCounter++;
                }
                // System.out.println("contTriangle: " + countTriangle + " examined nodes: " + nodeCounter);
                double aktCost = 0.0;
                //find the cheapest way to code this entry
                if (g.isNeighbor(currCol, row)) {        //edge
                    if(firstEdge && secondEdge && triangle){
                        aktCost = min(costClosingTriangle, costFirstEdge, costSecondEdge);
                    }
                    if(!firstEdge && secondEdge && triangle){
                        aktCost = Math.min(costSecondEdge, costClosingTriangle);
                    }
                    if(firstEdge && !secondEdge && triangle){
                        aktCost = Math.min(costFirstEdge, costClosingTriangle);
                    }
                    if(firstEdge && secondEdge && !triangle){
                        aktCost = Math.min(costFirstEdge, costSecondEdge);
                    }
                    if(firstEdge && !secondEdge && !triangle){
                        aktCost = costFirstEdge;
                    }
                    if(!firstEdge && secondEdge && !triangle){
                        aktCost = costSecondEdge;
                    }
                    if(!firstEdge && !secondEdge && triangle){
                        aktCost = costClosingTriangle;
                    }


                    }

                 else { //no edge
                    if(firstEdge && secondEdge && triangle){
                        aktCost = min(costNotClosingTriangle, costNotFirstEdge, costNotSecondEdge);
                    }
                    if(!firstEdge && secondEdge && triangle){
                        aktCost = Math.min(costNotSecondEdge, costNotClosingTriangle);
                    }
                    if(firstEdge && !secondEdge && triangle){
                        aktCost = Math.min(costNotFirstEdge, costNotClosingTriangle);
                    }
                    if(firstEdge && secondEdge && !triangle){
                        aktCost = Math.min(costNotFirstEdge, costNotSecondEdge);
                    }
                    if(firstEdge && !secondEdge && !triangle){
                        aktCost = costNotFirstEdge;
                    }
                    if(!firstEdge && secondEdge && !triangle){
                        aktCost = costNotSecondEdge;
                    }
                    if(!firstEdge && !secondEdge && triangle){
                        aktCost = costNotClosingTriangle;
                    }

                    
                }

                if(!firstEdge && !secondEdge && !triangle){
                        System.err.println("error when coding entry: " + row + " " + currCol);
                    }


               

                cc += aktCost;
                codingCosts[currCol][row] = aktCost;
                codingCosts[row][currCol] = aktCost;
//numEntries++; //ok
                row++;
                currCol++; //next entry
            }
            col++; //next diagoal
          //  System.out.println("Diagonal: " + col);
        }
        //System.out.println("numEntries: " + numEntries); //ok
        IO ea = new IO();
        ea.writeDoubleToMatlab(codingCosts, "cc.mat");
        int numEntries = (numVertices * numVertices-1)/2;
        double parameterCosts = 0.5 * lg2(numEntries);
        double cCosts = cc + firstDiag;
        double overall = cCosts + parameterCosts;
        System.out.println("coding costs: " + cCosts + " parameter Costs: " + parameterCosts);
        return overall;


    }

    public double min(double a, double b, double c) {
        return Math.min(Math.min(a, b), c);
    }

    public double codingWithTriangles(double probEdge, double probTriangle) {
        double[][] codingCosts = new double[numVertices][numVertices];
        double entropy = -probEdge * lg2(probEdge) + (1 - probEdge) * lg2(1 - probEdge);
        double firstDiag = (numVertices - 1) * entropy; //costs for the first diagonal
        double cc = 0.0;
        int col = 2;
        int row = 0;

        double probOther = 1 - probTriangle;
        // System.out.println("Diagonal: " + col);
        while (col < numVertices) {
            int currCol = col; //start with current column
            row = 0;
            while (currCol < numVertices) {
                boolean triangle = false;
                int colT = currCol - 1;
                int rowT = currCol - 1;
                System.out.println("entry: " + row + " " + currCol + " " + g.isNeighbor(currCol, row) + " -----------");
                int countTriangle = 0;
                int countOther = 0;
                int nodeCounter = 0;
                while (colT > row && rowT > row) {
                    if (g.isNeighbor(colT, row) && g.isNeighbor(rowT, currCol)) {
                        triangle = true;
                        countTriangle++;
                        System.out.println(" would form triangle with: " + rowT);
                    } else {
                        countOther++;

                    }
                    colT--;
                    rowT--;
                    nodeCounter++;
                }
                System.out.println("contTriangle: " + countTriangle + " countOther " + countOther + " examined nodes: " + nodeCounter);
                double aktCost = 0.0;

                //SECOND TRY: 8142 bit with 100 nodes and 0.9 probability for triangles
//                if (triangle) {
//                    if (g.isNeighbor(currCol, row)) {
//                        aktCost = -lg2(probTriangle);
//                    } else {
//                        aktCost = -lg2(probOther);
//                    }
//
//                } else {
//                    aktCost = -lg2(probEdge);
//                }
//                cc += aktCost;
//                codingCosts[row][currCol] = aktCost;
//                codingCosts[currCol][row] = aktCost;


                //FIRST TRY: This does not save anything! 8348 bits with 100 nodes and 0.9 probability for triangles - coding only with entropy: 4950 bit
                if ((triangle && g.isNeighbor(currCol, row)) || (!triangle && !g.isNeighbor(currCol, row))) {
                    double dd = -lg2(probTriangle);
                    cc += -lg2(probTriangle);
                    codingCosts[currCol][row] = dd;
                    codingCosts[row][currCol] = dd;

                } else {
                    double dd = -lg2(probOther);
                    cc += -lg2(probOther);

                    codingCosts[currCol][row] = dd;
                    codingCosts[row][currCol] = dd;
                }
                row++;
                currCol++; //next entry
            }
            col++; //next diagoal
            //System.out.println("Diagonal: " + col);
        }
        IO ea = new IO();
        ea.writeDoubleToMatlab(codingCosts, "cc.mat");
        return cc + firstDiag;


    }

    private double lg2(double d) {
        return Math.log(d) / Math.log(2.0);
    }
}
