/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rhwlab.chipseq.cluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.Random;
import java.util.TreeMap;
import org.rhwlab.bedformat.BedBase;

/**
 *
 * @author gevirl
 */
public class CoOccurance {

    String[] clusters;
    String[] tfs;
    short[] peakTF;
    int[] peakCluster;
    TreeMap<String,Integer> peakCounts= new TreeMap<>();
    Random rnd = new Random();

    public CoOccurance(File peakBedFile) throws Exception {
        
        TreeMap<String, Integer> clusterMap = new TreeMap<>();
        TreeMap<String, Short> tfMap = new TreeMap<>();
        ArrayList<Integer> peakClusterList = new ArrayList<>();
        ArrayList<Short> peakTFList = new ArrayList<>();
        

        BufferedReader reader = new BufferedReader(new FileReader(peakBedFile));
        String line = reader.readLine();
        while (line != null) {
            BedBase bed = new BedBase(line);
            String clusterName = bed.getValue(9);
            if (!clusterName.startsWith("Sing")) {
                Integer cIndex = clusterMap.get(clusterName);
                if (cIndex == null) {
                    cIndex = clusterMap.size();
                    clusterMap.put(clusterName, cIndex);
                }
                String tfName = bed.getName().split("_")[0];
                Short tfIndex = tfMap.get(tfName);
                if (tfIndex == null) {
                    tfIndex = (short) tfMap.size();
                    tfMap.put(tfName, tfIndex);
                }
                peakClusterList.add(cIndex);
                peakTFList.add(tfIndex);
                
                Integer count = peakCounts.get(tfName);
                if (count == null){
                    peakCounts.put(tfName, 1);
                } else {
                    peakCounts.put(tfName, count+1);
                }
            }
            line = reader.readLine();
        }

        peakTF = new short[peakTFList.size()];
        int i = 0;
        for (Short s : peakTFList) {
            peakTF[i] = s;
            ++i;
        }

        peakCluster = new int[peakClusterList.size()];
        i = 0;
        for (Integer c : peakClusterList) {
            peakCluster[i] = c;
            ++i;
        }
        reader.close();

        tfs = new String[tfMap.size()];
        for (Entry<String, Short> e : tfMap.entrySet()) {
            tfs[e.getValue()] = e.getKey();
        }

        clusters = new String[clusterMap.size()];
        for (Entry<String, Integer> e : clusterMap.entrySet()) {
            clusters[e.getValue()] = e.getKey();
        }
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

    public int pairCount(boolean[][] b, int tf1, int tf2) {
        int count = 0;
        for (int c = 0; c < clusters.length; ++c) {
            if (b[tf1][c] && b[tf2][c]) {
                ++count;
            }
        }
        return count;
    }

    public void reportPairsCount(PrintStream stream, boolean[][] b) {
        boolean first = true;
        for (int tf1 = 0; tf1 < tfs.length - 1; ++tf1) {
            for (int tf2 = tf1 + 1; tf2 < tfs.length; ++tf2) {
                int count = pairCount(b, tf1, tf2);
                if (!first) {
                    stream.print('\t');
                }
                stream.print(count);
                first = false;
            }
        }
        stream.println();
    }

    public void reportClusterTFs(PrintStream stream, boolean[][] b) {
        for (int c = 0; c < clusters.length; ++c) {
            stream.printf("%s", clusters[c]);
            for (int t = 0; t < tfs.length; ++t) {
                if (b[t][c]) {
                    stream.printf("\t%s", tfs[t]);
                }
            }
            stream.println();
        }

    }

    public void randomize(PrintStream stream, int nIter) {
        int nPairs = (tfs.length * (tfs.length - 1)) / 2;
        Welford[] w = new Welford[nPairs];
        for (int i = 0; i < nPairs; ++i) {
            w[i] = new Welford();
        }

        for (int i = 0; i < nIter; ++i) {
            System.out.println(i);
            int p = 0;
            boolean[][] b = tfByCluster(sampleWithReplace(peakCluster));
            for (int tf1 = 0; tf1 < tfs.length - 1; ++tf1) {
                for (int tf2 = tf1 + 1; tf2 < tfs.length; ++tf2) {
                    w[p].update((double) pairCount(b, tf1, tf2));
                    ++p;
                }
            }
        }

        boolean[][] notRandom = tfByCluster(peakCluster);
        int p = 0;
        for (int tf1 = 0; tf1 < tfs.length - 1; ++tf1) {
            for (int tf2 = tf1 + 1; tf2 < tfs.length; ++tf2) {
                double mu = w[p].getMean();
                double sd = Math.sqrt(w[p].getVariance());
                double c = (double) pairCount(notRandom, tf1, tf2);
                double z = (c - mu) / sd;
                stream.printf("%s\t%d\t%s\t%d\t%f\t%f\t%f\t%f\n", tfs[tf1],peakCounts.get(tfs[tf1]),
                        tfs[tf2],peakCounts.get(tfs[tf2]), mu, sd, c, z);
                ++p;
            }
        }
    }

    static public void main(String[] args) throws Exception {
        PrintStream stream = new PrintStream("/net/waterston/vol9/ChipSeqPipeline/test/AllWormPeaks.TF.noDups.clustered.pairs");
        File peakFile = new File("/net/waterston/vol9/ChipSeqPipeline/test/AllWormPeaks.TF.noDups.clustered.bed");
        CoOccurance co = new CoOccurance(peakFile);
        co.randomize(stream, 50000);

        //       boolean[][] b = co.tfByCluster(co.sampleWithReplace(co.peakCluster));
        //       co.reportPairsCount(stream, b);
        stream.close();
        int uiasd = 0;
    }
}
