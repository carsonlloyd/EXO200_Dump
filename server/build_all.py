#!/usr/bin/env python
README = """
This script should be run at SLAC where it has access to the EXO200 data.
This script finds alpha like events, and classifies them, and saves them for
further analsysis.

Phase2a started Oct 5th, ended April 16th (runs 2464-3564)
    This is marked by U-wire shaping time change?

AlphaIon run list 
"""

import ROOT
import sys
ROOT.gSystem.Load("libEXOUtilities")
#ROOT.gSystem.Load("libEXOCalibUtilities")
import datetime, pytz

##############################################################################################################
#                                              Constants
##############################################################################################################
fid_cut_cyl = (172.,5.,185.)  # radius, z_min, z_max
fid_cut_hex = (172.,5.,185.)  #apethum, z_min, z_max
cathode_anode_dist = 192. #mm
max_drift_time = cathode_anode_dist/1.711 + 6./2.041
# 1.711 mm/us in drift_vol and 2.041 mm/us in anode region

# Energy calibration database flavors (same database flavors used in 2014 0nu analysis)
fCalibrationDatabaseFlavors = ("2013-0nu-denoised", "2013-0nu-denoised", "vanilla")

##############################################################################################################
#                                                 MAIN
##############################################################################################################
def main():
    """
    This builds a concentrated dataset from the selected initial LB dataset
    Output is a single root file with entries corresponding to the cuts set here.
    """
    # Verbosity: 1=Selection Results, >1 is various debugging information
    verbose = 0
    print "build_all.py running with verbose=%s"%(str(verbose))
    if verbose:
        print "Fiducial Cut: ",fid_cut_hex,"(apethum, z_min, z_max)"
        print "Max Drift Distance = %.4f us"%(max_drift_time)

    tree = get_data_tree(list='All') # Golden All
     
    # We use the EXOFitting processed tree to get high-level physical quantities
    # like the anticorrelated energy, etc. 
    #ptree_file = ROOT.TFile(preprocessed_tree)
    #ptree = ROOT.Get("dataTree")
    #if verbose: print "Indexing EXOFitting PreProcessed tree"
    #ptree.BuildIndex("runNum", "eventNum")
    #if verbose: print " ...done"

    cuts = ""

    #There must be at least 1 scintillation cluster:
    #cuts = "@fScintClusters.size()>=1"
    #cuts = "(fScintClusters.GetCountsOnAPDPlane(0)+fScintClusters.GetCountsOnAPDPlane(1))>20000"

    # The minimum scintinlation counts must be > 20000 and <70000
    # I observe that three peaks presumable alphas are at 38500, 42200, and 55000
    # So Rn222=5.4MeV, Po218=6MeV, Po214=7.7MeV
    # calibrate:: y=mx+b,  m=6167, b=5198
    #cuts = "fScintClusters.fRawEnergy>20000 && fScintClusters.fRawEnergy<70000"
    #cuts += "&& fScintClusters.fRawEnergy>22000 && fScintClusters.fRawEnergy<80000"
    #cuts += " && Sum$(fAPDs.fRawCounts) > 8000"

    # Ignore Noise and Muon tagged events
    cuts +="fEventHeader.fTaggedAsNoise==0 && fEventHeader.fTaggedAsMuon==0"  

    # That's the last of the cuts, lets show the user what the cut looks like
    print "Applying Cuts to data: \n%s"%cuts

    #Draw is the fastest method to apply cuts, in the end what we want is a reduced data list
    # to perform a more targeted analysis...
    tree.Draw(">>+elist_alpha_canidates",cuts,"goff")
    elist_alpha_canidates = ROOT.gDirectory.Get("elist_alpha_canidates")
    print "There are %d events passing the initial cuts"%elist_alpha_canidates.GetN()

    #Now we have to look at events passing the cuts individually
    tf = ROOT.TFile("all.root","RECREATE")
    Rntree = tree.CloneTree(0)
    
    for i in range(elist_alpha_canidates.GetN()):
        # Print Progress
        if i%int(elist_alpha_canidates.GetN()/20) == 0:
            print "%d of %d"%(i,elist_alpha_canidates.GetN())

        #Grab the event data
        tree.GetEntry(elist_alpha_canidates.GetEntry(i))
        #ed = tree.EventBranch
        #if verbose>1: print_event_data(ed,verbose)

        #is_alphaish = check_alpha_like(ed,verbose)
        
        #is the event a fully reconstructed BiPo?
        #is_bipo = check_full_BiPo(ed,verbose)

        # Case1 (position matched Bi-Po)
        #is_case1 = check_case1(ed,verbose)
        #print "BiPo=%s, Case1=%s"%(is_bipo, is_case1) 
        #raw_input('<hit any key to continue>')
        #if is_bipo or is_alphaish:
            # Write the EventData of events which pass any of our selection criteria
            # to ROOT file
        Rntree.Fill()

    Rntree.AutoSave()

##############################################################################################################
#                                              Helper Functions
##############################################################################################################

def get_data_tree(list = 'All', test=False,rebuild=False):
    """
    Returns TChain of requested run list:
      AlphaIon (default) - Custom list made for AlphaIon work, golden + low purity runs
      Golden - The golden run list used for the 2014 0v paper
      GoldenDenoised - The golden run list used for the 2014 0v paper
      All - All data runs since phase 2 started
    """
    run_list = get_run_list(list, rebuild=rebuild)
    #Create TChain of all data runs:
    print "Making TChain..."
    tree = ROOT.TChain("tree")
    if test:
        i=0
        print "Test Mode, only using 1/50th of the data"
    for run in run_list.keys():
        if test:
            i+=1
            if i%50 != 0:
                continue
        for file in run_list[run]['file_paths']:
            tree.Add(file)
            print "tree.Add", run
    print "... Done. There are %d events in this dataset."%tree.GetEntries()
    #Build Index so we can associeate this with the pre-processed tree information
    #if verbose: print "Indexing ROOT tree"
    #tree.BuildIndex("fRunNumber", "fEventNumber")
    #if verbose: print " ...done"
    return tree

def get_run_list(list = 'All', rebuild=False):
    # Get all runs that we can do analysis on. This includes right after feeds where
    # there are more alphas, but detector responce isn't as great.
    if list == 'All':
        run_list = get_all_runs_list()

    # NOTE: THIS ISN'T FUNCTIONING YET
    #The conditions data stores when data_quality group marks parts of data runs BAD
    #Grabbing bad times in-case we want to cut these out...
    #bad_times = get_bad_time_data()
    return run_list


# ---ooo000OOO000ooo------ooo000OOO000ooo------ooo000OOO000ooo------ooo000OOO000ooo------ooo000OOO000ooo---
def get_all_runs_list():
    list = "masked"
    """
    This gets all runs in the data-catalogue that are masked, type=Data-Physics, and run>=7104 (Phase 2)
    """
    all_ds = ROOT.EXORunInfoManager.GetDataSet("Data/Processed/%s"%list,"runType==\"Data-Physics\" && run>=7104")
    run_info = {}
    for i in all_ds:
        run = i.GetRunNumber()
        print "get_all_runs_list", run
        run_info[run] = {}
        run_info[run]['file_paths'] = []
        for afile in i.GetRunFiles():
            run_info[run]['file_paths'].append(afile.GetFileLocation())
        run_data = ROOT.EXORunInfoManager.GetDataRunInfo(run)
        run_info[run]['exposure'] = run_data.FindMetaData('exposure').AsDouble()
        run_info[run]['start_time'] = run_data.FindMetaData('startTime').AsString()  
        run_info[run]['end_time'] = run_data.FindMetaData('endTime').AsString()
        run_info[run]['run_type'] = run_data.FindMetaData('runType').AsString()
    return run_info



if __name__ == '__main__':
    main()
