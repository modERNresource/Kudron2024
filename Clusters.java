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
import java.util.List;
import java.util.TreeMap;
import org.rhwlab.bedformat.BedBase;
import org.rhwlab.bedformat.Bed_9;
import org.rhwlab.gene.model.Annotation;
import org.rhwlab.gene.model.TSS;
import org.rhwlab.singlecell.expression.PrefixMaps;
import org.rhwlab.singlecell.expression.SingleCellExprMat;

/**
 *
 * @author gevirl
 */
public class Clusters {

    TreeMap<String, TreeMap<Integer, BedBase>> clusterMap = new TreeMap<>();

    public Clusters(File clusterFile) throws Exception {
        BufferedReader reader = new BufferedReader(new FileReader(clusterFile));
        String line = reader.readLine();
        while (line != null) {
            String[] cluster = line.split("\t");
            Bed_9 bed = new Bed_9(cluster);
            int loc = (bed.getThickEnd() + bed.getThickStart()) / 2;
            TreeMap<Integer, BedBase> chromoMap = clusterMap.get(bed.getChromosome());
            if (chromoMap == null) {
                chromoMap = new TreeMap<>();
                clusterMap.put(bed.getChromosome(), chromoMap);
            }
            chromoMap.put(loc, bed);
            line = reader.readLine();
        }
        reader.close();
    }

    // report the distance to the nearest cluster for each TSS
    public void reportNearestCluster(PrintStream stream, TSS tss, TreeMap<String, SingleCellExprMat> mat) {
        TreeMap<String, TreeMap<Integer, List<Annotation>>> allTSS = tss.allTSSs();
        for (String chromo : allTSS.keySet()) {
            TreeMap<Integer, BedBase> clusterBeds = clusterMap.get(chromo);
            TreeMap<Integer, List<Annotation>> chromoTSSMap = allTSS.get(chromo);
            for (Integer tssLoc : chromoTSSMap.keySet()) {
                Integer close = closestCluster(clusterBeds, tssLoc);
                if (close != null) {
                    Annotation annot = chromoTSSMap.get(tssLoc).get(0);
                    String wbGene = ((String) annot.getAttributeValue("Parent")).split(":")[1];

                    Double embExpr = mat.get("Emb").maxExpression(wbGene);
                    Double l2Expr = mat.get("L2").maxExpression(wbGene);
                    if (embExpr != null && l2Expr != null) {

                        String name = (String) annot.getAttributeValue("Name");
                        Object obj = annot.getAttributeValue("locus");
                        if (obj == null) {
                            obj = name;
                        }
                        stream.printf("%s\t%s\t%d\t%f\n", (String) obj, name, close, Math.max(embExpr,l2Expr));
                        int kjasfd = 0;
                    }
                }
                int xidxhf = 0;
            }
        }
    }

    // report the distance to the nearest TSS for each cluster
    public void reportNearestTSS(PrintStream stream, TSS tss) {
        for (String chromo : clusterMap.keySet()) {
            TreeMap<Integer, BedBase> chromoMap = clusterMap.get(chromo);
            for (Integer loc : chromoMap.keySet()) {
                BedBase bed = chromoMap.get(loc);
                Integer close = tss.closestTSS(chromo, loc);
                stream.printf("%s\t%d\n", bed.toString(), close);
            }
        }
    }

    static public Integer closestCluster(TreeMap chromoMap, int loc) {
        if (chromoMap != null) {
            Integer high = (Integer) chromoMap.higherKey(loc);
            Integer floor = (Integer) chromoMap.floorKey(loc);
            if (high != null) {
                if (floor != null) {
                    int right = high - loc;
                    int left = loc - floor;
                    if (right < left) {
                        return right;
                    } else {
                        return -left;
                    }
                } else {
                    return high - loc;
                }
            } else {
                return floor - loc;
            }
        }
        return null;
    }

    public static void main(String[] args) throws Exception {
        File gffFile = new File("/net/waterston/vol9/WS285/c_elegans.PRJNA13758.WS285.annotations.allwormbase.gff3");

        TSS tss = new TSS(gffFile);

        File file = new File("/net/waterston/vol9/ChipSeqPipeline/test/AllWormPeaks.TF.noDups.clusters.bed");
        Clusters cluster = new Clusters(file);

        TreeMap<String, SingleCellExprMat> mat = PrefixMaps.getWormExpressionMatrices();

        PrintStream stream = new PrintStream("/net/waterston/vol9/ChipSeqPipeline/test/nearestClusterDistanceToTSS.tsv");
        cluster.reportNearestCluster(stream, tss, mat);
        stream.close();

        stream = new PrintStream("/net/waterston/vol9/ChipSeqPipeline/test/AllWormPeaks.TF.noDups.clusters.tss.bed");
        cluster.reportNearestTSS(stream, tss);
        stream.close();
        int iausdfh = 0;
    }
}
