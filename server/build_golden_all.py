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
    verbose = 1
    print "CLbuild_golden.py running with verbose=%s"%(str(verbose))
    if verbose:
        print "Fiducial Cut: ",fid_cut_hex,"(apethum, z_min, z_max)"
        print "Max Drift Distance = %.4f us"%(max_drift_time)

    tree = get_data_tree(list='Golden')
     
    # We use the EXOFitting processed tree to get high-level physical quantities
    # like the anticorrelated energy, etc. 
    #ptree_file = ROOT.TFile(preprocessed_tree)
    #ptree = ROOT.Get("dataTree")
    #if verbose: print "Indexing EXOFitting PreProcessed tree"
    #ptree.BuildIndex("runNum", "eventNum")
    #if verbose: print " ...done"

    #There must be at least 1 scintillation cluster:
    #cuts = "@fScintClusters.size()>=1"
    cuts = "(fScintClusters.GetCountsOnAPDPlane(0)+fScintClusters.GetCountsOnAPDPlane(1))>20000"

    # The minimum scintinlation counts must be > 20000 and <70000
    # I observe that three peaks presumable alphas are at 38500, 42200, and 55000
    # So Rn222=5.4MeV, Po218=6MeV, Po214=7.7MeV
    # calibrate:: y=mx+b,  m=6167, b=5198
    #cuts = "fScintClusters.fRawEnergy>20000 && fScintClusters.fRawEnergy<70000"
    #cuts += "&& fScintClusters.fRawEnergy>22000 && fScintClusters.fRawEnergy<80000"
    #cuts += " && Sum$(fAPDs.fRawCounts) > 8000"

    # Ignore Noise and Muon tagged events
    cuts +=" && fEventHeader.fTaggedAsNoise==0 && fEventHeader.fTaggedAsMuon==0"  

    # That's the last of the cuts, lets show the user what the cut looks like
    print "Applying Cuts to data: \n%s"%cuts

    #Draw is the fastest method to apply cuts, in the end what we want is a reduced data list
    # to perform a more targeted analysis...
    tree.Draw(">>+elist_alpha_canidates",cuts,"goff")
    elist_alpha_canidates = ROOT.gDirectory.Get("elist_alpha_canidates")
    print "There are %d events passing the initial cuts"%elist_alpha_canidates.GetN()

    #Now we have to look at events passing the cuts individually
    tf = ROOT.TFile("GoldenMaskedData.root","RECREATE")
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

def get_data_tree(list = 'Golden', test=False,rebuild=False):
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
    print "... Done. There are %d events in this dataset."%tree.GetEntries()
    #Build Index so we can associeate this with the pre-processed tree information
    #if verbose: print "Indexing ROOT tree"
    #tree.BuildIndex("fRunNumber", "fEventNumber")
    #if verbose: print " ...done"
    return tree

def get_run_list(list = 'Golden', rebuild=False):
    # This retrieves a list of golden data runs from the RESTfull interface (not denoised)
    if list == 'Golden':
        run_list = get_golden_data_list()
    #   ---- or ----
    #
    # This gets the LB data from the EXO_Fitting DeNoised Run2abc file list
    elif list == 'EXO_Fitting':
        run_list = get_golden_data_from_file('BuildRunList/fileListRun2abcDenoised.dat')
    #   ---- or ----
    #
    # This gets the runs from the exofitting list, but then grabs the files from the database
    elif list == 'EXO_Fitting_runs':
        run_list = get_runs_from_fitting('BuildRunList/fileListRun2abcDenoised.dat',list='masked')
    #
    # Get all runs that we can do analysis on. This includes right after feeds where
    # there are more alphas, but detector responce isn't as great.
    elif list == 'All':
        run_list = get_all_runs_list()
    #
    #  ---- or ----
    # A run list of all Physics runs that are in 0v or were >1hour and passed hand-scan
    # This includes runs with low purity, but not with poor operational conditions
    elif list == 'AlphaIon':
        run_list = get_alphaion_run_list(rebuild)
    #
    # ---- or ----
    # Run 1 (runs 2464-3564) from EXOFitting run list
    # http://java.freehep.org/svn/repos/exo/show/EXO_Fitting/tags/EXO_Fitting_2ndAnalysis_v1.0/analysis/fileListUnmasked.dat?revision=HEAD
    elif list == 'Run1':
        run_list = get_runs_from_2nd_analysis_fitting('BuildRunList/EXOFitting_2ndAnalysis_fileListUnmasked.dat',list='masked')
    #
    # A run_list is a dictionary object that cointains at a minimum:
    #     run_info[run]['file_paths'] = []

    # NOTE: THIS ISN'T FUNCTIONING YET
    #The conditions data stores when data_quality group marks parts of data runs BAD
    #Grabbing bad times in-case we want to cut these out...
    #bad_times = get_bad_time_data()
    return run_list

# ---ooo000OOO000ooo------ooo000OOO000ooo------ooo000OOO000ooo------ooo000OOO000ooo------ooo000OOO000ooo---
def get_golden_data_list(list='masked',**kwargs):
    """
    This function gets the golden run list from the RESTful interface, but...
    we dont use the !!&$(@ RESTful interface to describe our run lists for analysis
    (YET)
    """
    if 'verbose' in kwargs.keys():
        verbose = int(kwargs['verbose'])
    run_start = 2464 # 6385 technically
    run_end   = 9097
    run_cut_str = "quality==\"GOLDEN\" && runType==\"Data-Physics\""
    if 'run_min' in kwargs.keys():
        run_start = int(kwargs['run_min'])
    if 'run_max' in kwargs.keys():
        run_end = int(kwargs['run_max'])
    run_cut_str += " && run>= %d && run<= %d"%(run_start,run_end)
    if 'verbose' in kwargs:
        print "Getting Run Data List:"
        print "  ",run_cut_str
    golden_ds = ROOT.EXORunInfoManager.GetDataSet("Data/Processed/%s"%list,run_cut_str)
    run_info = {}
    for i in golden_ds:
        run = i.GetRunNumber()
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
