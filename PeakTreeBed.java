/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rhwlab.chipseq.cluster;

import java.util.ArrayList;
import java.util.List;
import org.rhwlab.bedformat.NarrowPeakBedRecord;

/**
 *
 * @author gevirl
 */
public class PeakTreeBed extends NarrowPeakBedRecord {

    ArrayList<PeakTreeBed> children = new ArrayList<>();
    PeakTreeBed parent = null;

    public PeakTreeBed() {
        super();
    }

    public PeakTreeBed(String rec) {
        super(rec);
    }

    public PeakTreeBed(String[] tokens) {
        super(tokens);
    }

    public void shortenNames() {
        String currentName = this.getName().replace(".out", "");
        String[] tokens = currentName.split("\\.");
        this.setName(tokens[tokens.length - 1]);
    }

    public void addChild(PeakTreeBed child) {
        children.add(child);
        child.setParent(this);
    }

    public void setParent(PeakTreeBed p) {
        this.parent = p;
    }

    // is another record within this record
    public boolean isWithin(PeakTreeBed other) {
        if (other.getStart() >= this.getStart()) {
            if (other.getEnd() <= this.getEnd()) {
                return true;
            }
        }
        return false;
    }

    public boolean hasSiblings() {
        if (parent != null) {
            return this.parent.children.size() > 1;
        }
        return false; // roots do not have sibs
    }

    // return the root of this node
    public PeakTreeBed rootParent(){
        if (this.parent == null){
            return this;
        }
        return this.parent.rootParent();
    }
    // finds the ancestory that first gave rise to this node
    public PeakTreeBed birthParent() {
        if (this.parent == null) {
            return this;
        }
        if (parent.hasSiblings()) {
            return parent;
        }
        return parent.birthParent();
    }

    // unlink a record from its parent
    public void unlink() {
        if (this.parent != null) {
            parent.children.remove(this);
        }
    }

    public List<PeakTreeBed> terminalRecords() {
        ArrayList<PeakTreeBed> ret = new ArrayList<>();
        if (children.isEmpty()) {
            ret.add(this);
            return ret;
        }
        for (PeakTreeBed child : children) {
            ret.addAll(child.terminalRecords());
        }
        return ret;
    }

    public int terminalNodeCount() {
        if (children.isEmpty()) {
            return 1;
        }
        int count = 0;
        for (PeakTreeBed child : children) {
            count = count + child.terminalNodeCount();
        }
        return count;
    }
}
