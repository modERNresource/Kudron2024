/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rhwlab.chipseq.cluster;

import java.util.Random;
import java.util.concurrent.Callable;

/**
 *
 * @author gevirl
 */
public class CoOccurWorker implements Callable<int[]> {

    Random rnd = new Random();

    String[] clusters;
    String[] tfs;
    short[] peakTF;
    int[] peakCluster;
    int[] pairCounts;

    public CoOccurWorker(CoOccurMultiThread co) {
        this.clusters = co.clusters;
        this.tfs = co.tfs;
        this.peakTF = co.peakTF;  // the TF for each peak
        this.peakCluster = co.peakCluster;  // the cluster designation for each peak
        int nPairs = (tfs.length * (tfs.length - 1)) / 2;
        pairCounts = new int[nPairs];
    }

    @Override
    public int[] call() throws Exception {
//        return allPairsCount(sampleWithReplace(this.peakCluster));
        return allPairsCount(shuffle(this.peakCluster));
    }

    public int[] notRandomCounts() {
        return allPairsCount(this.peakCluster);
    }

    public int[] allPairsCount(int[] clusterIndex) {
        boolean[][] b = tfByCluster(clusterIndex);
        int p = 0;
        for (int tf1 = 0; tf1 < tfs.length - 1; ++tf1) {
            for (int tf2 = tf1 + 1; tf2 < tfs.length; ++tf2) {
                pairCounts[p] = pairCount(b, tf1, tf2);
                ++p;
            }
        }
        return pairCounts;
    }

    public boolean[][] tfByCluster(int[] clusterIndex) {
        boolean[][] ret = new boolean[tfs.length][];
        for (int t = 0; t < tfs.length; ++t) {
            ret[t] = new boolean[clusters.length];
        }

        for (int p = 0; p < clusterIndex.length; ++p) {
            short tfIndex = peakTF[p];
            ret[tfIndex][clusterIndex[p]] = true;
        }
        return ret;
    }

    public int[] sampleWithReplace(int[] v) {

        int[] ret = new int[v.length];
        for (int i = 0; i < ret.length; ++i) {
            int index = rnd.nextInt(v.length);
            ret[i] = v[index];
        }
        return ret;
    }

    public int[] shuffle(int[] v) {
        return shuffle(rnd, v);
    }

    static public int[] shuffle(Random rnd, int[] v) {
        int[] ret = new int[v.length];
        for (int i = 0; i < v.length; ++i) {
            int j = rnd.nextInt(i + 1);
            if (j != i) {
                ret[i] = ret[j];
            }
            ret[j] = v[i];
        }
        return ret;
    }

    static public int pairCount(boolean[][] b, int tf1, int tf2) {
        int count = 0;
        int nClusters = b[0].length;
        for (int c = 0; c < nClusters; ++c) {
            if (b[tf1][c] && b[tf2][c]) {
                ++count;
            }
        }
        return count;
    }

    static public void main(String[] args) {
        Random rnd = new Random();
        int[] v = new int[10];
        for (int i = 0; i < v.length; ++i) {
            v[i] = i;
        }
        int[] s = shuffle(rnd, v);
        for (int i = 0; i < v.length; ++i) {
            System.out.println(s[i]);
        }
        int isadfpoisd = 0;
    }

}
