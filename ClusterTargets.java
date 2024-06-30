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
import java.util.Map;
import java.util.TreeMap;
import org.rhwlab.bedformat.BedBase;
import org.rhwlab.gene.model.Annotation;
import org.rhwlab.gene.model.ChromosomeRemap;
import org.rhwlab.gene.model.FlyRemap;
import org.rhwlab.gene.model.ModelGFF;
import org.rhwlab.gene.model.ModelGFF_Fly;
import org.rhwlab.gene.model.ModelGFF_Worm;

/**
 *
 * @author gevirl
 */
public class ClusterTargets {

    ModelGFF gff;
    TreeMap<String, TreeMap<Integer, List<Annotation>>> tssMap = new TreeMap<>();  // chromo -> map<TSS location -> List of mRNA annotation>

    public ClusterTargets(ModelGFF gff) {
        this.gff = gff;
        for (Annotation mRNA : gff.getAnnotationsByType("mRNA")) {
            String chromo = mRNA.getChromosome();

            int loc = mRNA.getStart();
            if (mRNA.getStrand().equals("-")) {
                loc = mRNA.getEnd();
            }

            TreeMap<Integer, List<Annotation>> chrMap = tssMap.get(chromo);
            if (chrMap == null) {
                chrMap = new TreeMap<>();
                tssMap.put(chromo, chrMap);
            }

            List<Annotation> annots = chrMap.get(loc);
            if (annots == null) {
                annots = new ArrayList<>();
                chrMap.put(loc, annots);
            }
            annots.add(mRNA);
            /*            
            if (annots.size() > 4) {
                for (Annotation a : annots){
                    System.out.printf("%s\t%s\t%d\t%d\n", a.getChromosome(),a.getStrand(),a.getStart(),a.getEnd());
                }
                System.out.println();
                int isudhf = 0;
            }
             */
        }

    }

    public static int distanceToTSS(int loc, Annotation annot) {
        if (annot == null) {
            return Integer.MAX_VALUE;
        }
        if (annot.getStrand().equals("+")) {
            return Math.abs(annot.getStart() - loc);
        }
        return Math.abs(annot.getEnd() - loc);
    }

    public void reportTwoTargets(PrintStream stream,String chromo, BedBase bed, int loc) {
        TreeMap<Integer, List<Annotation>> chromoMap = this.tssMap.get(chromo);
        if (chromoMap != null) {
            if (bed.getName().equals("Singleton_14020")){
                int jasbdfh=0;
            }
            String clustName = bed.getName();
            Map.Entry<Integer, List<Annotation>> low = chromoMap.lowerEntry(loc);
            reportTargetAnnotation(stream, bed, low, loc);

            bed.setName(clustName);
            Map.Entry<Integer, List<Annotation>> high = chromoMap.ceilingEntry(loc);
            reportTargetAnnotation(stream, bed, high, loc);
        } else {
            System.out.println(bed.toString());
        }
    }

    public String reportTargetAnnotation(PrintStream stream, BedBase bed, Map.Entry<Integer, List<Annotation>> e, int loc) {
        if (e != null) {
            String clustName = bed.getName();

            // choose the longest transcript
            List<Annotation> annots = e.getValue();
            Annotation mRNA = annots.get(0);
            int dMax = mRNA.getEnd() - mRNA.getStart();
            for (int i = 1; i < annots.size(); ++i) {
                Annotation current = annots.get(i);
                int d = current.getEnd() - current.getStart();
                if (d > dMax) {
                    dMax = d;
                    mRNA = current;
                }
            }

            String geneName = "";
//            String[] geneNam((String) mRNA.getAttributeValue("Alias"));es = gff.nameTripletFromWBGene(mRNA.getGeneName());

            String transcript = ((String) mRNA.getAttributeValue("ID"));
            if (transcript.startsWith("Transcript:")) {
                transcript = transcript.split(":")[1];
            }
            String transcriptParent = ((String) mRNA.getAttributeValue("Parent"));
            Annotation gene = gff.getAnnotation("gene", transcriptParent);

            if (transcriptParent.startsWith("Gene:")) {
                String locus = ((String) gene.getAttributeValue("locus"));
                if (locus != null) {
                    geneName = locus;
                } else {
                    geneName = ((String) gene.getAttributeValue("sequence_name"));
                }
            } else {
                geneName = ((String) gene.getAttributeValue("ID"));
            }

            bed.setName(String.format("%s/%s", geneName, transcript));

            String status = "Out";
            if (mRNA.getStart() <= loc && loc <= mRNA.getEnd()) {
                status = "In";
            }

            String strand = mRNA.getStrand();
            int d = e.getKey() - loc;
            if (strand.equals("-")) {
                d = -d;
            }
            stream.printf("%s\t%s\t%s\t%s\t%d\n", bed.toString(), clustName, strand, status, d);
            return geneName;
        } else {
//            stream.printf("%s\tnull\tnull\tnull\tnull\t0\n",bed.toString());
            return null;
        }
    }

    public static void main(String[] args) throws Exception {
        File bedFile = new File("/net/waterston/vol9/ChipSeqPipeline/test/AllWormPeaks.TF.noDups.clusters.bed");
        bedFile = new File("/net/waterston/vol9/ChipSeqPipeline/test/AllFlyPeaks.TF.noDups.clusters.bed");
        if (args.length == 1) {
            bedFile = new File(args[0]);
        }
        File twoFile = new File(bedFile.getPath().replace(".bed", ".twoTarget.bed"));
        PrintStream twoStream = new PrintStream(twoFile);

        
        ModelGFF gff = null;
        Annotation.remapChromo = false;
        if (bedFile.getName().contains("Worm")) {
            File gffFile = new File("/net/waterston/vol9/WS285/c_elegans.PRJNA13758.WS285.annotations.allwormbase.gff3");
            gff = new ModelGFF_Worm(gffFile);
            Annotation.remapChromo = true;
        } else {
            File gffFile = new File("/net/waterston/vol9/dm6/dmel-all-r6.49.flybase.gff3");
            gff = new ModelGFF_Fly(gffFile);
        }

        

        ChromosomeRemap remap = null;
        if (bedFile.getName().contains("Fly")) {
            remap = new FlyRemap(gff);
        }

        ClusterTargets targets = new ClusterTargets(gff);

        BufferedReader reader = new BufferedReader(new FileReader(bedFile));
        String line = reader.readLine();
        while (line != null) {

            BedBase bed = new BedBase(line);
            String remapChromo = bed.getChromosome();
            if (remap != null) {
                remapChromo = remap.remap(bed.getChromosome());
            }
            int loc = Integer.valueOf(bed.getValue(6)) + 2;
            targets.reportTwoTargets(twoStream,remapChromo, bed, loc);
            line = reader.readLine();
        }
        int sdf = 0;
    }
/*
    // find gene targets for the clusters
    public static void mainSave(String[] args) throws Exception {
        File bedFile = new File("/net/waterston/vol9/ChipSeqPipeline/test/AllWormPeaks.TF.noDups.clusters.bed");
        bedFile = new File("/net/waterston/vol9/ChipSeqPipeline/test/AllFlyPeaks.TF.noDups.clusters.bed");
        if (args.length == 1) {
            bedFile = new File(args[0]);
        }

        File posFile = new File(bedFile.getPath().replace(".bed", ".posTarget.bed"));
        File negFile = new File(bedFile.getPath().replace(".bed", ".negTarget.bed"));
        File bothFile = new File(bedFile.getPath().replace(".bed", ".bothTarget.bed"));
        File noneFile = new File(bedFile.getPath().replace(".bed", ".noneTarget.bed"));
        File consvFile = new File(bedFile.getPath().replace(".bed", ".consvTarget.bed"));
        File p3File = new File(bedFile.getPath().replace(".bed", ".threeP.Target.bed"));
        File twoFile = new File(bedFile.getPath().replace(".bed", ".twoTarget.bed"));
        File distFile = new File(bedFile.getPath().replace(".bed", ".Target.dist"));

        ChromosomeRemap remap = null;
        File gffFile = new File("/net/waterston/vol9/dm6/dmel-all-r6.49.flybase.gff3");
        Annotation.remapChromo = false;
        if (bedFile.getName().contains("Worm")) {
            gffFile = new File("/net/waterston/vol9/WS285/c_elegans.PRJNA13758.WS285.annotations.allwormbase.gff3");
            Annotation.remapChromo = true;
        }
        ModelFromGFF gff = new ModelFromGFF(gffFile);
        if (bedFile.getName().contains("Fly")) {
            //           remap = new ChromosomeRemap(gff, new File("/net/waterston/vol9/dm6/dm6.chromAlias.txt"));
        }
        GeneTargets targets = new GeneTargets(gff);

        int[] failed = new int[3];
        PrintStream posStream = new PrintStream(posFile);
        PrintStream negStream = new PrintStream(negFile);
        PrintStream bothStream = new PrintStream(bothFile);
        PrintStream noneStream = new PrintStream(noneFile);
        PrintStream consvStream = new PrintStream(consvFile);
        PrintStream distStream = new PrintStream(distFile);
        PrintStream p3Stream = new PrintStream(p3File);
        PrintStream twoStream = new PrintStream(twoFile);
        BufferedReader reader = new BufferedReader(new FileReader(bedFile));
        String line = reader.readLine();
        int threePrime = 0;
        int total = 0;
        while (line != null) {
            ++total;
            BedBase bed = new BedBase(line);
            String chromo = bed.getChromosome();
            if (remap != null) {
                String remapChromo = remap.remap(chromo);
                if (remapChromo == null) {
                    int kljasef = 0;
                }
                bed.setChromosome(remapChromo);
            }
            if (gff.getChromosomeSize(bed.getChromosome()) != null) {
                if (bed.getName().equals("Singleton_14827")) {
                    int sjdf = 0;
                }
                int loc = Integer.valueOf(bed.getValue(6)) + 2;
                Annotation posAnnot = targets.posStrandTarget(bed.getChromosome(), loc);
                Annotation negAnnot = targets.negStrandTarget(bed.getChromosome(), loc);

                if (posAnnot != null && negAnnot == null) {
                    targets.reportTargetBed(bed, posAnnot, posStream);
                } else if (negAnnot != null && posAnnot == null) {
                    targets.reportTargetBed(bed, negAnnot, negStream);
                } else if (posAnnot == null && negAnnot == null) {
                    noneStream.println(bed.toString());
                } else {
                    targets.reportTargetBed(bed, posAnnot, bothStream);
                    targets.reportTargetBed(bed, negAnnot, bothStream);
                }

                Annotation consvAnnot = targets.consvTarget(bed.getChromosome(), loc);
                ++failed[targets.getFailed()];
                targets.reportTargetBed(bed, consvAnnot, consvStream);

                if (posAnnot != null || negAnnot != null) {
                    distStream.printf("%s\t%d\t%d\n", bed.getName(), bed.getScore(), Math.min(distanceToTSS(loc, posAnnot), distanceToTSS(loc, negAnnot)));
                }

                if (targets.threePrimeCluster(bed, loc, p3Stream)) {
                    ++threePrime;
                }

                targets.reportTwoTargets(twoStream, bed, loc);
            }
            line = reader.readLine();
        }
        reader.close();
        posStream.close();
        negStream.close();
        bothStream.close();
        noneStream.close();
        consvStream.close();
        distStream.close();
        p3Stream.close();
        twoStream.close();

        System.out.printf("%d\t%d\t%d\n", failed[0], failed[1], failed[2]);
        System.out.printf("3 prime clusters = %d from %d total clusters\n", threePrime, total);
        int iausdhf = 0;
    }
*/
}
