/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rhwlab.chipseq.cluster;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;
import org.rhwlab.bedformat.BedBase;
import org.rhwlab.chipseq.PeakCluster;
import org.rhwlab.chipseq.PeakClusterFile;
import org.rhwlab.gene.model.Annotation;
import org.rhwlab.gene.model.GeneIntervals;
import org.rhwlab.gene.model.ModelFromGFF;

/**
 *
 * @author gevirl
 */
public class FormBOW {

    TreeMap<String, TreeMap<String,Integer>> peakMap = new TreeMap<>();
    TreeSet<String> tfs = new TreeSet<>();
    TreeMap<String,Integer> tfIndex = new TreeMap<>();
    int nPeaks;

    public FormBOW(File peakFile,int lowerLimit,int upperLimit) throws Exception {
        BufferedReader reader = new BufferedReader(new FileReader(peakFile));
        String line = reader.readLine();
        while (line != null) {
            BedBase bed = new BedBase(line);
            String tf = bed.getName().split("_")[0];
            String cluster = bed.getValue(10);
            if (!cluster.startsWith("Single")) {
                TreeMap<String,Integer> tfMap = peakMap.get(cluster);
                if (tfMap == null) {
                    tfMap = new TreeMap<>();
                    peakMap.put(cluster, tfMap);
                }
                
                Integer count = tfMap.get(tf);
                if (count == null){
                    tfMap.put(tf, 1);
                } else {
                    tfMap.put(tf, count+1);
                }
            }
            line = reader.readLine();
        }
        reader.close();
        
        // remove any additional single tf clusters
        TreeSet<String> toRemove = new TreeSet<>();
        for (String cluster : peakMap.keySet()){
            TreeMap<String,Integer> tfMap = peakMap.get(cluster);
            if (tfMap.size()== 1){
                toRemove.add(cluster);
            } else if (lowerLimit < tfMap.size() && tfMap.size() > upperLimit){
                toRemove.add(cluster);
            }
        }
        for (String cluster : toRemove){
            peakMap.remove(cluster);
        }
        
        // count the tfs amd peaks
        nPeaks = 0;
        for (String cluster : peakMap.keySet()){
            TreeMap<String,Integer> tfMap = peakMap.get(cluster);
            for (String tf : tfMap.keySet()){
                Integer n = tfMap.get(tf);
                nPeaks = nPeaks + n;
                tfs.add(tf);
            }
        }
        
        int index = 1;
        for (String tf : tfs){
            tfIndex.put(tf,index);
            ++index;
        }
        int huasd=0;
    }
    
    public void saveBOW(File outFile)throws Exception {
        PrintStream bowStream = new PrintStream(outFile);
        PrintStream docStream = new PrintStream(outFile.getPath() + ".clusters");
        bowStream.println(peakMap.size());
        bowStream.println(tfs.size());
        bowStream.println(nPeaks);
        int doc = 1;
        for (String cluster : peakMap.keySet()){
            docStream.println(cluster);
            TreeMap<String,Integer> tfMap = peakMap.get(cluster);
            for (String tf : tfMap.keySet()){
                bowStream.printf("%d\t%d\t%d\n",doc,tfIndex.get(tf),tfMap.get(tf));
            }
            ++doc;
        }
        bowStream.close();
        docStream.close();
        
        PrintStream wordStream = new PrintStream(outFile.getPath() + ".tfs");
        for (String tf : tfs){
            wordStream.println(tf);
        }
        wordStream.close();
    }

    static public void main(String[] args) throws Exception {
        File peakFile = new File("/net/waterston/vol2/home/gevirl/Downloads/OptimalTFClusteredWormPeaks");
        FormBOW bow = new FormBOW(peakFile,1,100);

        File bowFile = new File("/net/waterston/vol9/ChipSeqPipeline/lda/OptimalTFClusteredWormPeaks.bow");
        bow.saveBOW(bowFile);
        
        int asidhf = 0;
    }

    static public void oldClusterFormat(String[] args) throws Exception {
        int pad = 1200;
        File dir = new File("/net/waterston/vol9/ChipSeqPipeline/LDA");

        File gffFile = new File("/net/waterston/vol9/References/WS280/c_elegans.PRJNA13758.WS280.annotations.wormbase.gff3");
        Annotation.remapChromo = true;   // puts chr into chromosome
        ModelFromGFF gff = new ModelFromGFF(gffFile);

        GeneIntervals geneIntervals = new GeneIntervals(gff, pad);

        PeakClusterFile clusterFile = new PeakClusterFile(new File(dir, "PeakClusters.200"));
        IntervalTreeMap map = clusterFile.getTreeMap();

        PrintStream stream = new PrintStream(new File(dir, "PeakClusters.200.bow"));
        PrintStream streamTF = new PrintStream(new File(dir, "PeakClusters.200.tflist"));
        PrintStream streamTarget = new PrintStream(new File(dir, "PeakClusters.200.targetlist"));

        // get all the TFs
        TreeSet<String> tfs = new TreeSet<>();
        for (Object obj : map.values()) {
            PeakCluster cluster = (PeakCluster) obj;
            tfs.addAll(cluster.getTFs());
        }
        ArrayList<String> tfList = new ArrayList<>();

        for (String tf : tfs) {
            streamTF.println(tf);
            tfList.add(tf);
        }
        streamTF.close();

        // count the number of documents and peaks
        int w = 0;
        int doc = 0;
        int noTarget = 0;
        int multiTarget = 0;
        int singleTarget = 0;
        int twoTargets = 0;
        for (Object obj : map.values()) {
            PeakCluster cluster = (PeakCluster) obj;
            Interval key = new Interval(cluster.getChromosome(), cluster.getMeanPosition(), cluster.getMeanPosition());
            Collection genes = geneIntervals.getMap().getOverlapping(key);
            doc = doc + genes.size();
            if (genes.size() == 0) {
                ++noTarget;
            }
            if (genes.size() == 1) {
                ++singleTarget;
            }
            if (genes.size() == 2) {
                ++twoTargets;
            }
            if (genes.size() > 2) {
                ++multiTarget;
            }
            for (String tf : cluster.getTFs()) {
                Short count = cluster.getCount(tf);
                if (count != null) {
                    w = w + count;
                }

            }
        }

        // print the header lines
        stream.println(doc);
        stream.println(tfList.size());
        stream.println(w);

        doc = 1;
        for (Object obj : map.values()) { // loop through the clusters 
            PeakCluster cluster = (PeakCluster) obj;
            Interval key = new Interval(cluster.getChromosome(), cluster.getMeanPosition(), cluster.getMeanPosition());
            Collection genes = geneIntervals.getMap().getOverlapping(key);
            for (Object geneObj : genes) {
                Annotation gene = (Annotation) geneObj;
                streamTarget.println(gene.getGeneName());
                for (String tf : cluster.getTFs()) {
                    Short count = cluster.getCount(tf);
                    if (count != null) {
                        int v = tfList.indexOf(tf);
                        stream.printf("%d %d %d\n", doc, v + 1, count);  // document and vocabulary are one-based
                    }
                }
                ++doc;
            }
        }
        stream.close();
        streamTF.close();
        streamTarget.close();
        System.out.printf("No Target: %d , SingleTarget: %d , Two Targets: %d , MultiTarget: %d\n", noTarget, singleTarget, twoTargets, multiTarget);
    }

}
