#This file must be in sync with runAnalyzer3.py
#differences must only be for the Radius  parameter
import FWCore.ParameterSet.Config as cms
OutputFileName = "Analyzer3.root"
numEventsToRun = 10
# Cone radii to calibrate
calib_cone_radii = (0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0)
#####change radius cone
def adjust_cone_radius(R,pfJetsJetrecoPrototype):
#
	rtag = "R%d" % (R*100.0)
#
	newPFReco = pfJetsJetrecoPrototype.clone(
			recoScaleCalcPeak = cms.PSet(
				Class = cms.string("ConstDouble"),
				value = cms.double(R)
				),
			recoScaleCalcJet = cms.PSet(
				Class = cms.string("ConstDouble"),
				value = cms.double(R)
				)
			)
	scale = R*0.17/0.5
	if scale<0.09:
			scale=0.09
	newPFReco.fixedScale = cms.double(scale)
	newPFReco_module = "fftPFJet" + rtag
#
	newPFCorr = configure_fftjet_correction_producer(
								 (jetCorrectionSequenceTag,jetCorrectionSequenceTag2), 
									newPFReco_module)
	newPFCorr.calculatePileup = cms.bool(False)
	newPFCorr.subtractPileup  = cms.bool(False)
#
	newPFCorr_module = "fftPFJetCorr" + rtag
#
	return (
			(newPFReco, newPFReco_module),
			(newPFCorr, newPFCorr_module),
			)
##### import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *
##---------  Load standard Reco modules ------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
##----- Global tag: conditions database ------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START53_V7E::All'
############################################
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(numEventsToRun)
)
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
'/store/data/Run2012A/Photon/AOD/22Jan2013-v1/20000/007A5A14-1069-E211-8BDD-0025905964B4.root',
) )
#############fftjet config####################################
# apply corrections
database = 'sqlite_file:fftjet_corr.db'
tableName = "L2L3"
tableCategory = "PF"
jetCorrectionSequenceTag = "PF0"
tableName2 = "FFT"
tableCategory2 = "L2Res"
jetCorrectionSequenceTag2 = "PF1"
# Configure database access
from JetMETCorrections.FFTJetModules.fftjetcorrectionesproducer_cfi import *
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.CondDBCommon.connect = database
# Configure event setup to work with jet corrections
configure_fftjet_pooldbessource(process, jetCorrectionSequenceTag)
configure_fftjet_pooldbessource(process, jetCorrectionSequenceTag2)
config, esProducer = configure_L2L3_fftjet_esproducer(jetCorrectionSequenceTag,
                                                      tableName, tableCategory)
setattr(process, esProducer, config)
config, esProducer = configure_L2Res_fftjet_esproducer(jetCorrectionSequenceTag2,
                                                      tableName2, tableCategory2)
setattr(process, esProducer, config)
# Configure the jet corrections module
from RecoJets.FFTJetProducers.fftjetcorrectionproducer_cfi import *
# Do charged hadron subtraction
from RecoJets.FFTJetProducers.fftjetpfpileupcleaner_cfi import *
#Using the table by FFTJetProducer:
#----------------------------------
from JetMETCorrections.FFTJetModules.fftjetlookuptableesproducer_cfi import *
##pucorrections 
etaSequenceTag = "EtaFlatteningFactors"
rhoSequenceTag = "PileupRhoCalibration"
configure_fftjetlut_pooldbessource(process, etaSequenceTag)
config, esProducer = configure_fftjetlut_esproducer(etaSequenceTag)
setattr(process, esProducer, config)
configure_fftjetlut_pooldbessource(process, rhoSequenceTag)
config, esProducer = configure_fftjetlut_esproducer(rhoSequenceTag)
setattr(process, esProducer, config)

fftjet_pileup_grid_pf_calibrated = cms.PSet(
    nEtaBins = cms.uint32(256),
    etaMin = cms.double(-math.pi*2),
    etaMax = cms.double(math.pi*2),
    nPhiBins = cms.uint32(128),
    phiBin0Edge = cms.double(0.0),
    title = cms.untracked.string("FFTJet Pileup Grid")
)
from RecoJets.FFTJetProducers.fftjetpileupprocessor_pfprod_cfi import *
from RecoJets.FFTJetProducers.fftjetpileupestimator_pf_cfi import *
#Using the table by FFTJetPileupProcessor (on top of normal config):
#-------------------------------------------------------------------
fftjetPileupProcessorPf.src = cms.InputTag('fftjetcleanup')
fftjetPileupProcessorPf.etaFlatteningFactors = cms.vdouble()
fftjetPileupProcessorPf.convolverMinBin = cms.uint32(22)
fftjetPileupProcessorPf.convolverMaxBin = cms.uint32(234)
fftjetPileupProcessorPf.pileupEtaPhiArea = cms.double(2*math.pi*4*math.pi*((234-22)/256.0))
fftjetPileupProcessorPf.nScales = cms.uint32(1)
fftjetPileupProcessorPf.minScale = cms.double(0.1)
fftjetPileupProcessorPf.maxScale = cms.double(0.1)
fftjetPileupProcessorPf.GridConfiguration = fftjet_pileup_grid_pf_calibrated
fftjetPileupProcessorPf.flatteningTableRecord = cms.string(fftjet_lut_types[etaSequenceTag].LUTRecord)
fftjetPileupProcessorPf.flatteningTableName = cms.string("FFTEtaFlatteningFactorsTable")
fftjetPileupProcessorPf.flatteningTableCategory = cms.string("Pileup")
fftjetPileupProcessorPf.loadFlatteningFactorsFromDB = cms.bool(True)

#Using the table by FFTJetPileupEstimator (on top of normal config):
#-------------------------------------------------------------------
fftjetPileupEstimatorPf.inputLabel = cms.InputTag("pileupprocessor", "FFTJetPileupPF")
fftjetPileupEstimatorPf.cdfvalue = cms.double(0.4)
fftjetPileupEstimatorPf.uncertaintyZones = cms.vdouble()
fftjetPileupEstimatorPf.calibrationCurve = cms.PSet(
    Class = cms.string("Polynomial"),
    c0 = cms.double(0.0)
)
fftjetPileupEstimatorPf.uncertaintyCurve = cms.PSet(
    Class = cms.string("Polynomial"),
    c0 = cms.double(0.0)
)
fftjetPileupEstimatorPf.calibTableRecord = cms.string(fftjet_lut_types[rhoSequenceTag].LUTRecord)
fftjetPileupEstimatorPf.calibTableCategory = cms.string("Pileup")
fftjetPileupEstimatorPf.uncertaintyZonesName = cms.string("FFTPileupRhoUncertaintyZones")
fftjetPileupEstimatorPf.calibrationCurveName = cms.string("FFTPileupRhoCalibrationTable")
fftjetPileupEstimatorPf.uncertaintyCurveName = cms.string("FFTPileupRhoUncertaintyTable")
fftjetPileupEstimatorPf.loadCalibFromDB = cms.bool(True)

from RecoJets.FFTJetProducers.fftjetpatrecoproducer_cfi import *
# Pattern recognition for PFJets
fftjetPFPatReco = fftjetPatrecoProducer.clone()
fftjetPFPatReco.src = cms.InputTag('fftjetcleanup')
fftjetPFPatReco.jetType = cms.string('PFJet')
#fftjetPFPatReco.insertCompleteEvent = cms.bool(False)
fftjetPFPatReco.InitialScales = fftjet_patreco_scales_50
fftjetPFPatReco.calculateClusterRadii = cms.bool(True)
fftjetPFPatReco.calculateClusterSeparations = cms.bool(True)
# The Jet producer module
from RecoJets.FFTJetProducers.fftjetproducer_cfi import *
fftjetPFJetMaker = fftjetJetMaker.clone()
fftjetPFJetMaker.InitialScales = fftjet_patreco_scales_50
fftjetPFJetMaker.src = cms.InputTag('fftjetcleanup')
fftjetPFJetMaker.jetType = cms.string('PFJet')
fftjetPFJetMaker.resolution = cms.string('fixed')
fftjetPFJetMaker.maxIterations = cms.uint32(1000)
fftjetPFJetMaker.fixedScale = cms.double(0.1)
fftjetPFJetMaker.PeakSelectorConfiguration = fftjet_peak_selector_allpass
fftjetPFJetMaker.treeLabel = cms.InputTag("fftjetPFPatReco", "FFTJetPatternRecognition")
fftjetPFJetMaker.calculatePileup = cms.bool(True)
fftjetPFJetMaker.subtractPileup = cms.bool(True)
##subtract pile up 
fftjetPFJetMaker.PileupGridConfiguration = fftjet_pileup_grid_pf_calibrated
fftjetPFJetMaker.pileupDensityCalc = cms.PSet(
    Class = cms.string("PileupGrid2d"),
    Grid2d = fftjet_pileup_grid_pf_calibrated,
    rhoFactor = cms.double(0.0)
)
fftjetPFJetMaker.pileupTableRecord = cms.string(fftjet_lut_types[etaSequenceTag].LUTRecord)
fftjetPFJetMaker.pileupTableName = cms.string("FFTPileupRhoEtaDependenceTable")
fftjetPFJetMaker.pileupTableCategory = cms.string("Pileup")
fftjetPFJetMaker.loadPileupFromDB = cms.bool(True)
#############end fftjet config################################
###declare process 
process.fftjetcleanup    = fftjetPfPileupCleaner
process.pileupprocessor  = fftjetPileupProcessorPf
process.pileupestimator  = fftjetPileupEstimatorPf
process.fftjetPFPatReco  = fftjetPFPatReco

process.myseq  = cms.Sequence (
									process.fftjetcleanup*
									process.pileupprocessor*
									process.pileupestimator*
									process.fftjetPFPatReco)

for R in calib_cone_radii:
    for config, modulename in adjust_cone_radius(R,fftjetPFJetMaker):
        setattr(process, modulename, config)
        process.myseq *= getattr(process, modulename)
process.p = cms.Path(process.myseq) 
###write the outputstream 
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(OutputFileName),
#    outputCommands = cms.untracked.vstring('drop *')
    outputCommands = cms.untracked.vstring(
										 'keep *_*_*_*'
     									)                   
	)
#process.options.wantSummary = cms.untracked.bool(False)
## let it run
process.p = cms.Path(process.myseq)

