/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rhwlab.chipseq.cluster;

import java.io.PrintStream;
import java.util.List;
import java.util.TreeMap;
import org.rhwlab.gene.model.Annotation;
import org.rhwlab.gene.model.ModelGFF;

/**
 *
 * @author gevirl
 */
public class PeakTarget {

    Annotation mRNA;
    String mRNA_name;
    String mRNA_ID;
    String genePos;
    String geneID;
    String geneName;
    String transcriptPos;
    int geneTSSs;
    int mRNA_No;
    String transcriptLocation;
    int geneTranscripts;
    List<Annotation> geneExons;
    List<Annotation> transcriptExons;

    public PeakTarget(Annotation mRNA, ModelGFF gff, int peakLoc) {
        this.mRNA = mRNA;
        transcriptLocation = String.format("%s:%d:%d:%s", mRNA.getChromosome(), mRNA.getStart(), mRNA.getEnd(), mRNA.getStrand());
        mRNA_name = (String) mRNA.getAttributeValue("Name");
        mRNA_ID = (String) mRNA.getAttributeValue("ID");

        geneID = (String) mRNA.getAttributeValue("Parent");

        geneExons = gff.geneExons(geneID);
        transcriptExons = gff.getChildren(mRNA_ID, "exon");

        Annotation geneAnnot = gff.getAnnotation("gene", geneID);
        geneTranscripts = gff.getChildren(geneID, "mRNA").size();

        if (mRNA.getStrand().equals("+")) {
            genePos = positiveStrandPos(peakLoc,geneAnnot,geneExons);
            transcriptPos = positiveStrandPos(peakLoc,mRNA,transcriptExons);

        } else {
            genePos = negativeStrandPos(peakLoc,geneAnnot,geneExons);
            transcriptPos = negativeStrandPos(peakLoc,mRNA,transcriptExons);
        }

/*
        if (exonNo != null) {
            peakPos = String.format("E%d", exonNo);
        } else {
            Annotation firstExon = exons.get(0);
            Annotation lastExon = exons.get(exons.size() - 1);
            if (peakLoc < firstExon.getStart()) {
                if (mRNA.getStrand().equals("+")) {
                    peakPos = "5P" + genePos;
                } else {
                    peakPos = "3P" + genePos;
                }
            } else if (lastExon.getEnd() < peakLoc) {
                if (mRNA.getStrand().equals("+")) {
                    peakPos = "3P" + genePos;
                } else {
                    peakPos = "5P" + genePos;
                }
            } else {
                List<Annotation> introns = gff.getChildren(mRNA_ID, "intron");
                Integer intronNo = within(peakLoc, introns);
                peakPos = String.format("I%d", intronNo);

            }
        }
*/
        String transcriptParent = ((String) mRNA.getAttributeValue("Parent"));
        Annotation gene = gff.getAnnotation("gene", transcriptParent);
        String id = ((String) gene.getAttributeValue("ID"));
        if (id.contains("WBGene")) {
            String locus = ((String) gene.getAttributeValue("locus"));
            if (locus != null) {
                geneName = locus;
            } else {
                geneName = ((String) gene.getAttributeValue("sequence_name"));
            }
        } else {
            geneName = ((String) gene.getAttributeValue("Name"));
        }

        TreeMap<Integer, Annotation> geneTSSMap = gff.uniqueTSS(geneID);
        geneTSSs = geneTSSMap.size();

        int i = 1;
        for (Annotation annot : geneTSSMap.values()) {
            String annotID = (String) annot.getAttributeValue("ID");
            if (annotID.equals(mRNA_ID)) {
                mRNA_No = i;
                break;
            }
            ++i;
        }
    }
    
    public String getGeneName(){
        return this.geneName;
    }
    public Annotation getMRNA(){
        return this.mRNA;
    }

    public String getGeneID(){
        return this.geneID;
    }
    public static String positiveStrandPos(int peakLoc, Annotation annot, List<Annotation> exons) {
        String peakPos;
        if (peakLoc < annot.getStart()) {
            peakPos = "5P";
        } else if (peakLoc > annot.getEnd()) {
            peakPos = "3P";
        } else {
            Integer e = inAnnotation(peakLoc, exons);
            if (e != null) {
                peakPos = String.format("E%d", e);
            } else {
                peakPos = String.format("I%d", betweenAnnotation(peakLoc, exons));
            }
        }
        return peakPos;
    }

    public static String negativeStrandPos(int peakLoc, Annotation annot, List<Annotation> exons) {
        String peakPos;
        if (peakLoc < annot.getStart()) {
            peakPos = "3P";
        } else if (peakLoc > annot.getEnd()) {
            peakPos = "5P";
        } else {
            Integer e = inAnnotation(peakLoc, exons);
            if (e != null) {
                peakPos = String.format("E%d", 1 + exons.size() - e);
            } else {
                peakPos = String.format("I%d", exons.size() - betweenAnnotation(peakLoc, exons));
            }
        }
        return peakPos;
    }

    public static Integer inAnnotation(int loc, List<Annotation> annots) {
        Integer ret = null;
        int i = 1;
        for (Annotation annot : annots) {
            if (annot.getStart() <= loc && loc <= annot.getEnd()) {
                ret = i;
                break;
            }
            ++i;
        }
        return ret;
    }

    public static Integer betweenAnnotation(int loc, List<Annotation> annots) {
        Integer ret = null;

        for (int i = 1; i < annots.size(); ++i) {
            if (annots.get(i - 1).getEnd() < loc && loc < annots.get(i).getStart()) {
                ret = i;
                break;
            }
        }
        return ret;
    }

    public String getTargetID() {
        return String.format("%s::%s", this.mRNA_name, this.genePos);
    }

    public static void reportHeader(PrintStream stream) {
        stream.println("TargetID\tTranscriptCoord\tTranscript\tTranscriptLoc\tGene\tGeneLoc\tGeneExons\tTranscriptNo\tUniqueGeneTSSs\tGeneTranscripts");
    }

    public void report(PrintStream stream) {
        stream.printf("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n", 
                getTargetID(), transcriptLocation, mRNA_name, transcriptPos, geneName, genePos, geneExons.size(), mRNA_No, geneTSSs, geneTranscripts);
    }

}
