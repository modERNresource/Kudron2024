/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rhwlab.chipseq.cluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.TreeMap;
import org.rhwlab.singlecell.expression.PrefixMaps;
import org.rhwlab.singlecell.expression.SingleCellExprMat;

/**
 *
 * @author gevirl
 */
public class FormPriors {
    public static void main(String[] args) throws Exception{
        File dir = new File("/net/waterston/vol9/ChipSeqPipeline/LDA");
        
        ArrayList<String> tfList = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(new File(dir,"PeakClusters.200.tflist")));
        String tf = reader.readLine();
        while (tf != null ){
            tfList.add(tf);
            tf = reader.readLine();
        }
        reader.close();
        
        TreeMap<String,SingleCellExprMat> mats =  PrefixMaps.getWormExpressionMatrices();
        
        // cell type (topic) - 
    }
}
