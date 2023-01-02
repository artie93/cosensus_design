from pyrosetta import *
from Bio import AlignIO
from statistics import mode
import argparse
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation

pyrosetta.init()


scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

def pack_relax(pose, scorefxn):

    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.RestrictToRepacking())
    # Set up a MoveMapFactory
    mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
    mmf.all_bb(setting=True)
    mmf.all_bondangles(setting=True)
    mmf.all_bondlengths(setting=True)
    mmf.all_chi(setting=True)
    mmf.all_jumps(setting=True)
    mmf.set_cartesian(setting=True)

    ## Print informations about structure before apply fast relax
    # display_pose = pyrosetta.rosetta.protocols.fold_from_loops.movers.DisplayPoseLabelsMover()
    # display_pose.tasks(tf)
    # display_pose.movemap_factory(mmf)
    # display_pose.apply(pose) 

    fr = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1)
    fr.cartesian(True)
    fr.set_task_factory(tf)
    fr.set_movemap_factory(mmf)
    fr.min_type("lbfgs_armijo_nonmonotone")
    fr.apply(pose)
    return pose


def consensus_design(pose, consensus, scorefxn, design):
    
    pose = pack_relax(pose, scorefxn)

    
    for position in range(len(consensus)):
        
        if consensus[position] != '-':
            
            print("mutating position"+str(position))
            
            #posi =  pose.pdb_info().number(position+1)
            posi =  position+1
            amino = consensus[position]
            
            #Select Mutate Position
            mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
            mut_posi.set_index(posi)
            #Select Neighbor Position
            nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
            nbr_selector.set_focus_selector(mut_posi)
            nbr_selector.set_include_focus_in_subset(True)
            # Select No Design Area
            not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mut_posi)
            # The task factory accepts all the task operations
            tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
            # These are pretty standard
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
        
            # Disable Packing
            prevent_repacking_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
            prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True )
            tf.push_back(prevent_subset_repacking)
        
            # Disable design
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),not_design))
        
            # Enable design
            aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
            aa_to_design.aas_to_keep(amino)
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(aa_to_design, mut_posi))
        
            # Create Packer
            packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()
            packer.task_factory(tf) 
            packer.apply(pose)
    
    ### Relax
    
    if design == True:
        ################### Design residues with no significant consensus
        
        mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        
        
        ### Calculate N gap
        n_gap = 0
        for position in range(len(consensus)):
            if consensus[position] == '-':
                n_gap = n_gap+1
                #posi = pose.pdb_info().number(position+1)
                posi = position+1
                mut_posi.append_index(posi)
                #### Print selected residues for mut posi
                #print(pyrosetta.rosetta.core.select.get_residues_from_subset(mut_posi.apply(pose)))
        
        if n_gap != 0:
            # Select Neighbor Position
            nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
            nbr_selector.set_focus_selector(mut_posi)
            nbr_selector.set_include_focus_in_subset(True)
            #print(pyrosetta.rosetta.core.select.get_residues_from_subset(nbr_selector.apply(pose)))
            
            # Select No Design Area
            not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mut_posi)
            #### Print residues to NOT design
            #print(pyrosetta.rosetta.core.select.get_residues_from_subset(not_design.apply(pose)))
            
            # The task factory accepts all the task operations
            tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
            
            # These are pretty standard
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
            
            # Disable Packing
            prevent_repacking_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
            prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True )
            tf.push_back(prevent_subset_repacking)
            
            # Disable design
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
                pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),not_design))
            
            # Enable design
            aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
            aa_to_design.aas_to_keep("ACDEFGHIKLMNPQRSTVWY")
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(aa_to_design, mut_posi))
            
            # Create Packer
            packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)
            packer.task_factory(tf)
            
            packer.apply(pose)
            
        ### Second relax
            
    pose = pack_relax(pose, scorefxn)

    
    return pose

def get_consensus(alignment, cv_thresh, pose):
    ### Generate consensus sequence
    ## simple_mode -> doesn't use thresholds
    al = alignment.get_alignment_length()
    consensus = []
    for pos in range(al):
        column = []
        if alignment[0][pos] != '-':
            column = [align[pos] for align in alignment if align[pos] != '-']
            ### verify percentage and thresholds
            mode_count = 0
            for i in range(len(column)):
                if column[i] == mode(column) :
                    mode_count = mode_count + 1
                conservation = mode_count/len(column)
            
            if conservation  > cv_thresh:
                consensus.append(mode(column))
            else:
                consensus.append('-')
    
    ### correct numbering
    pdb_residues = []
    for i in range(1,pose.total_residue()):
        pdb_residues.append(pose.pdb_info().number(i))

    consensus_fixed = []
    for i in pdb_residues:
        consensus_fixed.append(consensus[i-1])
    
    return consensus_fixed



if __name__ == "__main__":
   
    parser = argparse.ArgumentParser(description='')
    
    parser.add_argument('--pdb', type=str, required=True)
    parser.add_argument('--alignment', type=str, required=True)
    parser.add_argument('--out', type=str, required=True)
    parser.add_argument('--cons_thresh', type=str, required=True)
    parser.add_argument('--design', type=str, required=True)

    args = parser.parse_args()

    alignment = AlignIO.read(open(args.alignment), "clustal")
    design = args.design

    al = alignment.get_alignment_length()

    pose = pose_from_pdb(args.pdb)
    

    #### Conservation threshold
    cv_thresh = float(args.cons_thresh)
    
    ### Generate consensus sequence
    consensus = get_consensus(alignment, cv_thresh, pose)
    
    pose = consensus_design(pose, consensus, scorefxn, design)    
    
    # remove
    print(str(consensus))

    pose.dump_pdb(args.out)

