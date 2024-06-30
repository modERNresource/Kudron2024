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
import java.util.TreeMap;
import org.rhwlab.bedformat.BedBase;
import org.rhwlab.gene.model.ChromosomeRemap;
import org.rhwlab.gene.model.ModelGFF;
import org.rhwlab.gene.model.ModelGFF_Worm;
import org.rhwlab.gene.model.WormRemap;

/**
 *
 * @author gevirl
 */
public class AhringerEnhancers {
    static public void toBed(File tsv,File bedFile)throws Exception{
        
        PrintStream stream = new PrintStream(bedFile);
        BufferedReader reader = new BufferedReader(new FileReader(tsv));
        String line = reader.readLine();
        while (line != null){
            
            String[] tokens = line.split("\t");
            String id = tokens[0];
            String chromo = id.split(":")[0];
            String[] coords =  id.split(":")[1].split("-");
            int start = Integer.valueOf(coords[0]);
            int end = Integer.valueOf(coords[1]);
            int loc = (start + end)/2;
            String[] bedTokens = new String[8];
            bedTokens[0] = chromo;
            bedTokens[1] = coords[0];
            bedTokens[2] = coords[1];
            bedTokens[3] = id;
            bedTokens[4] = "0";
            bedTokens[5] = ".";
            bedTokens[6] = String.format("%d",loc-2 );
            bedTokens[7] = String.format("%d",loc+2 );
            BedBase bed = new BedBase(bedTokens);
            stream.println(bed.toString());
            line = reader.readLine();
        }
        reader.close();  
        stream.close();
    }
    static public void main(String[] args)throws Exception{
        File tsv  = new File("/net/waterston/vol9/ChipSeqPipeline/Worm_enhancers_not_overlapped_by_modERN_cluster_nohead.tsv");
        File bedFile = new File(tsv.getPath().replace(".tsv", ".bed"));
        toBed(tsv,bedFile);
        
        File distFile = new File(bedFile.getPath().replace(".bed", ".distance.bed"));
        File altFile = new File(bedFile.getPath().replace(".bed", ".altdist.bed"));
        File targetFile = new File(bedFile.getPath().replace(".bed", ".peakTSSs"));
        
        File gffFile = new File("/net/waterston/vol9/WS285/c_elegans.PRJNA13758.WS285.annotations.allwormbase.gff3");
        ModelGFF gff = new ModelGFF_Worm(gffFile);
        TSS tss = new TSS(gff);
        
        ChromosomeRemap remap = new WormRemap();

        String heads = "Chromo\tStart\tEnd\tEnhancerID\tScore\tStrand\tApexStart\tApexEnd\tTargetID\tDistanceToTSS";

        TreeMap<String, PeakTarget> targets = tss.reportTargets(bedFile, null, distFile, altFile,heads);

        // report the targets
        PrintStream stream = new PrintStream(targetFile);
        PeakTarget.reportHeader(stream);
        for (PeakTarget target : targets.values()) {
            target.report(stream);
        }
        stream.close();

          
    }
}
