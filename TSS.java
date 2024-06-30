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
import java.util.List;
import java.util.TreeMap;
import org.rhwlab.bedformat.BedBase;
import org.rhwlab.gene.model.Annotation;
import org.rhwlab.gene.model.ChromosomeRemap;
import org.rhwlab.gene.model.FlyRemap;
import org.rhwlab.gene.model.ModelGFF;
import org.rhwlab.gene.model.ModelGFF_Fly;
import org.rhwlab.gene.model.ModelGFF_Worm;
import org.rhwlab.gene.model.WormRemap;

/**
 *
 * @author gevirl
 */
public class TSS {

    ModelGFF gff;
    TreeMap<String, TreeMap<Integer, Annotation>> tssMap = new TreeMap<>();  // chromo -> TSS loc -> longest mRNA annotation

    public TSS(ModelGFF gff) {
        this.gff = gff;
        for (Annotation mRNA : gff.getAnnotationsByType("mRNA")) {
            String chromo = mRNA.getChromosome();

            // location of the tss
            int loc = mRNA.getStart();
            if (mRNA.getStrand().equals("-")) {
                loc = mRNA.getEnd();
            }

            TreeMap<Integer, Annotation> chrMap = tssMap.get(chromo);
            if (chrMap == null) {
                chrMap = new TreeMap<>();
                tssMap.put(chromo, chrMap);
            }

            Annotation annot = chrMap.get(loc);
            if ((annot == null) || (mRNA.getLength() > annot.getLength())) {
                chrMap.put(loc, mRNA);
            }
        }

    }

    public Integer closestTSS(String chromo, int loc) {
        Integer ret = null;
        TreeMap<Integer, Annotation> chromoMap = tssMap.get(chromo);
        if (chromoMap != null) {
            Integer high = chromoMap.ceilingKey(loc);
            Integer low = chromoMap.lowerKey(loc);
            if (high == null) {
                ret = low;
            } else if (low == null) {
                ret = high;
            } else {
                if (Math.abs(high - loc) < Math.abs(low - loc)) {
                    ret = high;
                } else {
                    ret = low;
                }
            }
        }
        return ret;
    }

    // return the higher and lower TSS if exists
    public List<Integer> closeTSSs(String chromo, int loc) {
        List<Integer> ret = new ArrayList<>();
        if (chromo == null) {
            int jj = 0;
        }
        TreeMap<Integer, Annotation> chromoMap = tssMap.get(chromo);
        if (chromoMap != null) {
            Integer close = chromoMap.ceilingKey(loc);
            if (close != null) {
                ret.add(close);
            }
            close = chromoMap.lowerKey(loc);
            if (close != null) {
                ret.add(close);
            }
        }
        return ret;
    }

    public Annotation getTranscript(String chromo, int tss) {
        TreeMap<Integer, Annotation> chromoMap = tssMap.get(chromo);
        if (chromoMap != null) {
            return chromoMap.get(tss);
        }
        return null;
    }

    public static int distanceToTSS(int loc, Annotation annot) {
        if (annot == null) {
            return Integer.MAX_VALUE;
        }
        if (annot.getStrand().equals("+")) {
            return loc - annot.getStart();
        }
        return annot.getEnd() - loc;
    }

    // reports the nerest TSS to each bed record (peak or cluster) in a bed file
    // *** this has been replaced with reportTargets below
    public void nearestTSS(File bedFile, ChromosomeRemap remap, File distFile) throws Exception {
        PrintStream distStream = new PrintStream(distFile);
        BufferedReader reader = new BufferedReader(new FileReader(bedFile));
        String line = reader.readLine();
        while (line != null) {
            BedBase bed = new BedBase(line);
            if (bed.getStart() == 519963) {
                int jhh = 0;
            }
            String chromo = bed.getChromosome();
            if (remap != null) {
                chromo = remap.remap(bed.getChromosome());
            }
            int loc = Integer.valueOf(bed.getValue(6)) + 2;
            Integer close = closestTSS(chromo, loc);
            if (close != null) {

                Annotation mRNA = getTranscript(chromo, close);
                String mRNA_name = (String) mRNA.getAttributeValue("Name");
                String geneID = (String) mRNA.getAttributeValue("Parent");
                Annotation geneAnnot = gff.getAnnotation("gene", geneID);
                List<Annotation> all = gff.getChildren(geneID, "mRNA");

                TreeMap<Integer, Annotation> geneTSSMap = gff.uniqueTSS(geneID);

                int i = 0;
                for (Integer t : geneTSSMap.keySet()) {
                    ++i;
                    if (t.intValue() == close.intValue()) {
                        break;
                    }
                }
                if (mRNA.getStrand().equals("-")) {
                    i = geneTSSMap.size() - i + 1;
                }

                int d = loc - mRNA.getStart();
                if (mRNA.getStrand().equals("-")) {
                    d = mRNA.getEnd() - loc;
                }

                String peakStatus = "E";
                if (mRNA.getStart() <= loc && loc <= mRNA.getEnd()) {
                    peakStatus = "I";
                }
                distStream.printf("%s\t%d\t%d\t%s\t%d\t%s\t%s\n", bed.toString(), geneTSSMap.size(), i, peakStatus, d, mRNA_name, mRNA.getStrand());
            } else {
                System.out.println(bed.toString());
            }

            line = reader.readLine();
        }
        reader.close();
        distStream.close();
    }

    public TreeMap<String, PeakTarget> reportTargets(File bedFile, ChromosomeRemap remap, File distFile, File altFile, String head) throws Exception {
        PrintStream distStream = new PrintStream(distFile);
        PrintStream altStream = new PrintStream(altFile);
        PrintStream tsvDistStream = null;
        PrintStream tsvAltStream = null;

        boolean rename = bedFile.getName().contains("clusters");
        
        if (bedFile.getName().contains("ranked") || rename) {
            tsvDistStream = new PrintStream(distFile.getPath().replace(".bed", ".tsv"));
            tsvAltStream = new PrintStream(altFile.getPath().replace(".bed", ".tsv"));
            tsvDistStream.println(head);
            tsvAltStream.println(head);
        }
        TreeMap<String, PeakTarget> ret = new TreeMap<>();
        BufferedReader reader = new BufferedReader(new FileReader(bedFile));
        String line = reader.readLine();
        while (line != null) {
            BedBase bed = new BedBase(line);
            String chromo = bed.getChromosome();

            if (remap != null) {
                chromo = remap.remap(bed.getChromosome());
            }
//            int loc = Integer.valueOf(bed.getValue(6)) + 2;
            int loc = bed.getNarrowPeakApexLoc();
            List<Integer> closes = closeTSSs(chromo, loc);  // get the higher and lower genomic coordinates of the nearest TSSs

            if (closes.size() > 0) {
                if (closes.size() == 1) {
                    saveLocationDistance(bed, chromo, closes.get(0), distStream, tsvDistStream, ret, rename); // only one possible target - loc is at end of chromosome
                } else {
                    PeakTarget lowTarget = getTarget(chromo, closes.get(0), loc);
                    PeakTarget highTarget = getTarget(chromo, closes.get(1), loc);
                    int lowDist = distanceToTSS(loc, lowTarget.getMRNA());
                    int highDist = distanceToTSS(loc, highTarget.getMRNA());
                    if (Math.abs(lowDist) <= Math.abs(highDist)) {
                        saveLocationDistance(bed, chromo, closes.get(0), distStream, tsvDistStream, ret, rename);
                        if (!lowTarget.getGeneID().equals(highTarget.getGeneID())) {
                            saveLocationDistance(bed, chromo, closes.get(1), altStream, tsvAltStream, ret, rename);
                        }
                    } else {
                        saveLocationDistance(bed, chromo, closes.get(1), distStream, tsvDistStream, ret, rename);
                        if (!lowTarget.getGeneID().equals(highTarget.getGeneID())) {
                            saveLocationDistance(bed, chromo, closes.get(0), altStream, tsvAltStream, ret, rename);
                        }
                    }

                }
            }
            line = reader.readLine();
        }
        reader.close();
        distStream.close();
        altStream.close();
        if (tsvDistStream != null) {
            tsvDistStream.close();
            tsvAltStream.close();
        }
        return ret;
    }

    public PeakTarget getTarget(String chromo, int closeTSS, int loc) {
        Annotation mRNA = getTranscript(chromo, closeTSS);
        PeakTarget target = new PeakTarget(mRNA, gff, loc);
        return target;
    }

    public void saveLocationDistance(BedBase bed, String chromo, int closeTSS,
            PrintStream distStream, PrintStream tsvStream, TreeMap<String, PeakTarget> targetMap, boolean rename) {
        
        int loc = bed.getNarrowPeakApexLoc();
        Annotation mRNA = getTranscript(chromo, closeTSS);
        PeakTarget target = new PeakTarget(mRNA, gff, loc);
        targetMap.put(target.getTargetID(), target);

        if (rename) {
            String curName = bed.getName();
            bed.setName(target.getGeneName());
            distStream.printf("%s\t%s\t%d\t%s\n", bed.toString(), target.getTargetID(), distanceToTSS(loc, mRNA), curName);
            if (tsvStream != null){
                tsvStream.printf("%s\t%s\t%d\t%s\n", bed.toString(), target.getTargetID(), distanceToTSS(loc, mRNA), curName);
            }
            bed.setName(curName);
        } else {
            distStream.printf("%s\t%s\t%d\n", bed.toString(), target.getTargetID(), distanceToTSS(loc, mRNA));
            if (tsvStream != null) {
                tsvStream.printf("%s\t%s\t%d\n", bed.toString(), target.getTargetID(), distanceToTSS(loc, mRNA));
            }
        }

    }
// reports the nearest TSS to the records in a bed file

    public static void main(String[] args) throws Exception {

        
        File bedFile = new File("/net/waterston/vol9/ChipSeqPipeline/AllWormPeaks.TF.noDups.clusters.bed");
//        bedFile = new File("/net/waterston/vol9/ChipSeqPipeline/AllFlyPeaks.TF.noDups.clustered.bed");
        if (args.length == 1) {
            bedFile = new File(args[0]);
        }
        
        String head = null;
        if (bedFile.getName().contains("ranked")){
            head = "chrom\tchromStart\tchromEnd\texperiment\tscore\tstrand\tsignalValue\tsignalRank\tpValue\tpeak\tmetaPeak\tTF\ttarget\tdistance";
        }
        if (bedFile.getName().contains("clusters")){
            head = "chrom\tchromStart\tchromEnd\tTF\tnumPeaks\tstrand\tthickStart\tthickEnd\titemRgb\ttarget\tdistance\tmetaPeak";
        }
    
        
        File distFile = new File(bedFile.getPath().replace(".bed", ".distance.bed"));
        File altFile = new File(bedFile.getPath().replace(".bed", ".altdist.bed"));
        File targetFile = new File(bedFile.getPath().replace(".bed", ".peakTSSs"));
        File tsvTargetFile = new File(bedFile.getPath().replace(".bed", ".target.tsv"));

        Annotation.remapChromo = false;
        ModelGFF gff = null;
        if (bedFile.getName().contains("Worm")) {
            gff = new ModelGFF_Worm(new File("/net/waterston/vol9/WS285/c_elegans.PRJNA13758.WS285.annotations.allwormbase.gff3"));
        } else {
            gff = new ModelGFF_Fly(new File("/net/waterston/vol9/dm6/dmel-all-r6.49.flybase.gff3"));
        }

        TSS tss = new TSS(gff);

        ChromosomeRemap remap;
        if (bedFile.getName().contains("Fly")) {
            remap = new FlyRemap(gff);
        } else {
            remap = new WormRemap();
        }
        //       tss.nearestTSS(bedFile,remap,distFile);

        TreeMap<String, PeakTarget> targets = tss.reportTargets(bedFile, remap, distFile, altFile, head);

        // report the target descriptions
        PrintStream stream = new PrintStream(targetFile);
        PrintStream tsvStream = new PrintStream(tsvTargetFile);
        PeakTarget.reportHeader(stream);
        PeakTarget.reportHeader(tsvStream);
        for (PeakTarget target : targets.values()) {
            target.report(stream);
            target.report(tsvStream);
        }
        stream.close();
        tsvStream.close();
        int sdkjfm = 0;
    }
}
