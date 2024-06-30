/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rhwlab.chipseq.cluster;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import htsjdk.samtools.util.IntervalTreeMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;
import org.rhwlab.bedformat.BedBase;
import org.rhwlab.bedformat.BedInCluster;
import org.rhwlab.bedformat.Bed_9;
import org.rhwlab.bedformat.NarrowPeakBedRecord;

/**
 *
 * @author gevirl
 */
public class MetaPeakTrees {

    TreeMap<String, PeakTreeBed[]>[] records;  // (chromo -> array of bed records ) [indexed by the cutoff]

    public MetaPeakTrees(List<File> files) throws Exception {
        records = new TreeMap[files.size()];
        int i = 0;
        for (File file : files) {
            records[i] = readBed(file);
            ++i;
        }
        for (i = 1; i < records.length; ++i) {
            linkRecords(records[i - 1], records[i]);
        }
    }

    final public TreeMap<String, PeakTreeBed[]> readBed(File bedFile) throws Exception {
        TreeMap<String, PeakTreeBed[]> ret = new TreeMap<>();

        TreeMap<String, List<PeakTreeBed>> map = new TreeMap<>();
        BufferedReader reader = new BufferedReader(new FileReader(bedFile));
        String line = reader.readLine();
        while (line != null) {
            PeakTreeBed bed = new PeakTreeBed(line);
            bed.shortenNames();

            List<PeakTreeBed> list = map.get(bed.getChromosome());
            if (list == null) {
                list = new ArrayList<>();
                map.put(bed.getChromosome(), list);
            }
            list.add(bed);
            line = reader.readLine();
        }
        reader.close();

        // make lists into sorted arrays for each chromosome
        for (String chromo : map.keySet()) {
            List<PeakTreeBed> list = map.get(chromo);
            // sort the list
            list.sort(new Comparator() {
                @Override
                public int compare(Object o1, Object o2) {
                    PeakTreeBed rec1 = (PeakTreeBed) o1;
                    PeakTreeBed rec2 = (PeakTreeBed) o2;
                    return Integer.compare(rec1.getStart(), rec2.getStart());
                }
            });
            ret.put(chromo, list.toArray(new PeakTreeBed[0]));
        }
        return ret;
    }

    final public void linkRecords(TreeMap<String, PeakTreeBed[]> parent, TreeMap<String, PeakTreeBed[]> child) {
        for (String chromo : parent.keySet()) {
            PeakTreeBed[] childRecs = child.get(chromo);
            if (childRecs != null) {
                linkRecords(parent.get(chromo), childRecs);
            }
        }
    }

    final public void linkRecords(PeakTreeBed[] parent, PeakTreeBed[] child) {
        int c = 0;
        int p = 0;
        while (p < parent.length) {
            if (parent[p].isWithin(child[c])) {
                parent[p].addChild(child[c]);
                ++c;
                if (c == child.length) {
                    return;
                }
            } else {
                ++p;
            }
        }
    }

    // returns all the top records in the tree
    public List<PeakTreeBed> topRecs() {
        List<PeakTreeBed> ret = new ArrayList<>();

        TreeMap<String, PeakTreeBed[]> topMap = records[0];
        for (String chromo : topMap.keySet()) {
            PeakTreeBed[] topRecs = topMap.get(chromo);
            for (PeakTreeBed topRec : topRecs) {
                ret.add(topRec);
            }
        }
        return ret;
    }

    // trims terminal nodes that exist for only one level
    // it is trimmed if a terminal node has a sibling
    public List<PeakTreeBed> trimmedTerminal() {
        List<PeakTreeBed> ret = new ArrayList<>();
        for (PeakTreeBed topRec : topRecs()) {

            List<PeakTreeBed> termRecs = topRec.terminalRecords();
            for (PeakTreeBed rec : termRecs) {
                if (!rec.hasSiblings()) {
                    ret.add(rec);
                } else {
                    rec.unlink();
                }
            }
        }
        return ret;
    }

    public IntervalTreeMap formClusters() throws Exception {
        PrintStream dumpStream = new PrintStream("AllWormPeaks.TF.noDups.macs2.dump");
        IntervalTreeMap ret = new IntervalTreeMap<>();

        // get the trimmed terminal list of nodes
        List<PeakTreeBed> termList = trimmedTerminal();

        // move up the tree to birth parents
        List<PeakTreeBed> birthParents = new ArrayList<>();
        for (PeakTreeBed rec : termList) {
            birthParents.add(rec.birthParent());
        }

        // organize by root parent
        TreeMap<PeakTreeBed, List<PeakTreeBed>> rootMap = new TreeMap<>();
        for (PeakTreeBed birthBed : birthParents) {
            PeakTreeBed root = birthBed.rootParent();
            List<PeakTreeBed> list = rootMap.get(root);
            if (list == null) {
                list = new ArrayList<>();
                rootMap.put(root, list);
            }
            list.add(birthBed);
        }

        // expand the birth parents to fill the space of their root
        for (PeakTreeBed root : rootMap.keySet()) {

            List<PeakTreeBed> nodeList = rootMap.get(root);
            if (nodeList.size() > 1) {
                PeakTreeBed first = nodeList.get(0);
                first.setStart(root.getStart());

                PeakTreeBed lower = first;
                for (int i = 1; i < nodeList.size(); ++i) {
                    PeakTreeBed higher = nodeList.get(i);
                    int mid = (lower.getEnd() + higher.getStart()) / 2;
                    lower.setEnd(mid);
                    higher.setStart(mid + 1);
                    lower = higher;
                }
                lower.setEnd(root.getEnd());
            }

            for (PeakTreeBed bed : nodeList) {
                if (bed.getStart() >= bed.getEnd()) {
                    int iasudhf = 0;
                }
            }
        }

        // color the clusters
        //       String[] colors = {"0,0,255", "0,255,0", "255,0,0"};
        String[] colors = {"0,0,255", "38,155,38", "255,0,0"};
        int c = 0;
        for (PeakTreeBed rec : birthParents) {
            Bed_9 cluster = new Bed_9(rec.getTokens());
            cluster.setItemRgb(colors[c]);
            ++c;
            if (c == colors.length) {
                c = 0;
            }
            dumpStream.printf("%s,%d,%d\t%s,%d,%d\n", rec.getName(), rec.getStart(), rec.getEnd(), cluster.getName(), cluster.getStart(), cluster.getEnd());
            Interval interval = new Interval(cluster.getChromosome(), cluster.getStart(), cluster.getEnd(), false, cluster.getName());
            if (cluster.getName().equals("2_narrowPeak13021")) {
                int hjsdjds = 0;
            }
            ret.put(interval, cluster);
        }
        dumpStream.close();
        return ret;
    }

    public static BedInCluster assignToCluster(NarrowPeakBedRecord peak, IntervalTreeMap clusters, String noClusterColor) {

        Bed_9 bed9Peak = new Bed_9(peak.getTokens());
        bed9Peak.setThickStart(peak.getPeakLocation() - 2);
        bed9Peak.setThickEnd(peak.getPeakLocation() + 2);

        BedInCluster ret = new BedInCluster(bed9Peak.getTokens());

        IntervalTree intervalTree = clusters.debugGetTree(peak.getChromosome());
        if (intervalTree != null) {
            Bed_9 cluster = bestChoice(peak, intervalTree);
            if (cluster != null) {
                ret.setItemRgb(cluster.getItemRgb());
                ret.setCluster(cluster.getName());
            } else {
                ret.setItemRgb(noClusterColor);
                ret.setCluster(null);
            }
        }
        return ret;
    }

    static public Bed_9 bestChoice(NarrowPeakBedRecord peak, IntervalTree intervalTree) {
        Bed_9 cluster = null;
        int apex = peak.getPeakLocation();

        IntervalTree.Node node = intervalTree.minOverlapper(apex, apex);
        if (node != null) {
            cluster = (Bed_9) node.getValue();    // apex in a cluster
        }

        if (cluster == null) {
            // apex not in a cluster - use full width of peak to find all overlapping clusters
            Iterator iter = intervalTree.overlappers(peak.getStart(), peak.getEnd());
            ArrayList<Bed_9> list = new ArrayList<>();
            while (iter.hasNext()) {
                IntervalTree.Node obj = (IntervalTree.Node) iter.next();
                list.add((Bed_9) obj.getValue());
            }
            if (!list.isEmpty()) {
                if (list.size() == 1) {
                    cluster = list.get(0);
                } else {
                    // overlaps more than one cluster - pick closest
                    cluster = list.get(0);
                    int minD = distance(peak, cluster);
                    for (int i = 1; i < list.size(); ++i) {
                        int d = distance(peak, list.get(i));
                        if (d < minD) {
                            cluster = list.get(i);
                            minD = d;
                        }
                    }
                }
            }
        }
        return cluster;
    }

    static int distance(NarrowPeakBedRecord peak, Bed_9 cluster) {
        int l = peak.getPeakLocation();
        return Math.min(Math.abs(l - cluster.getStart()), Math.abs(l - cluster.getEnd()));
    }

    // cluster the TF peaks using the MACS2 tree clustering, fly and worm
    // input is the TF noDups bed files
    // output is three files:
    // 1) clustered peaks bed files
    // 2) ranked and clustered peaks bed files
    // 3) the defined clusters
    static public void main(String[] args) throws Exception {
        String[] bases = {"AllFlyPeaks.TF.noDups", "AllWormPeaks.TF.noDups"};
//             String[] bases = {"AllWormPeaks.TF.noDups"};
        File dir = new File("/data/www/site/waterston/html/ChipSeqPipeline");
//        File dir = new File("/net/waterston/vol9/ChipSeqPipeline");
        for (String fileBase : bases) {
            File peakFile = new File(dir, fileBase + ".bed");
            File clusteredPeaks = new File(peakFile.getPath().replace(".bed", ".clustered.bed"));
            File clusterFile = new File(peakFile.getPath().replace(".bed", ".clusters.bed"));
            File rankedFile = new File(peakFile.getPath().replace(".bed", ".ranked.bed"));

            List<File> files = new ArrayList<>();
            for (int i = 2; i <= 50; ++i) {
                files.add(new File(dir, String.format("%s.macs2.%d.bed", fileBase, i)));
            }
            MetaPeakTrees trees = new MetaPeakTrees(files);

            IntervalTreeMap clusters = trees.formClusters();

            String black = "0,0,0";

            // assign each peak to a cluster or leave it unclustered and report to file
            TreeMap<String, List<BedBase>> clusterMap = new TreeMap<>();
            String singleton = "Singleton";
            ArrayList<Bed_9> singletons = new ArrayList<>();
            int c = 1;
            TreeSet<String> chromosomes = new TreeSet<>();
            PrintStream stream = new PrintStream(clusteredPeaks);
            PrintStream rankedStream = new PrintStream(rankedFile);
            PrintStream tsvrankedStream = new PrintStream(rankedFile.getPath().replace(".bed", ".tsv"));
            String rankedHead = "chrom\tchromStart\tchromEnd\texperiment\tscore\tstrand\tsignalValue\tsignalRank\tpValue\tpeak\tmetaPeak\tTF";
            tsvrankedStream.println(rankedHead);
            
            BufferedReader reader = new BufferedReader(new FileReader(peakFile));
            String line = reader.readLine();
            while (line != null) {
                NarrowPeakBedRecord peak = new NarrowPeakBedRecord(line);
                if (peak.getChromosome().equals("chrIII") && peak.getStart() == 5693510 && peak.getEnd() == 5694000) {
                    int hhh = 0;
                }
                chromosomes.add(peak.getChromosome());
                BedInCluster assignedPeak = assignToCluster(peak, clusters, black);
                if (assignedPeak.getCluster() == null) {
                    String name = String.format("%s_%d", singleton, c);
                    assignedPeak.setCluster(name);
                    ++c;

                    // make a new cluster for the singleton and add it to the tree
                    Bed_9 singBed = new Bed_9(assignedPeak.getTokens());
                    singBed.setName(name);
                    singBed.setScore(1);
                    singletons.add(singBed);

                }
                List<BedBase> bedList = clusterMap.get(assignedPeak.getCluster());
                if (bedList == null) {
                    bedList = new ArrayList<>();
                    clusterMap.put(assignedPeak.getCluster(), bedList);
                }
                bedList.add(assignedPeak);
                stream.println(assignedPeak.toString());

                
                rankedStream.printf("%s\t%s\t%s\n", peak.toString(), assignedPeak.getCluster(),peak.getTF());
                tsvrankedStream.printf("%s\t%s\t%s\n", peak.toString(), assignedPeak.getCluster(),peak.getTF());
                line = reader.readLine();
            }
            reader.close();
            rankedStream.close();
            tsvrankedStream.close();
            stream.close();

            // add in the singletons
            for (Bed_9 singBed : singletons) {
                Interval interval = new Interval(singBed.getChromosome(), singBed.getStart(), singBed.getEnd());
                clusters.put(interval, singBed);
            }

            // report the clusters  
            String clusterHead = "chrom\tchromStart\tchromEnd\tmetaPeak\tnumPeaks\tstrand\tthickStart\tthickEnd\titemRgb";
            PrintStream cStream = new PrintStream(clusterFile);
            PrintStream tsvcStream = new PrintStream(clusterFile.getPath().replace(".bed", ".tsv"));
            tsvcStream.println(clusterHead);
            for (String chromo : chromosomes) {
                IntervalTree contigTree = clusters.debugGetTree(chromo);
                if (contigTree != null) {
                    Iterator iter = contigTree.iterator();
                    while (iter.hasNext()) {
                        Node node = (Node) iter.next();
                        BedBase bed = (BedBase) node.getValue();
                        int offset = Integer.valueOf(bed.getValue(9));

                        List<BedBase> list = clusterMap.get(bed.getName());
                        if (list != null) {
                            int n = list.size();
                            if (n > 1000) {
                                n = 1000;
                            }
                            bed.setScore(n);
                            bed.setValue(6, Integer.toString(bed.getStart() + offset - 2));
                            bed.setValue(7, Integer.toString(bed.getStart() + offset + 2));
                            cStream.println(bed.toString());
                            tsvcStream.println(bed.toString());
                        }

                    }
                }
            }
            cStream.close();
            tsvcStream.close();
        }
    }
}
