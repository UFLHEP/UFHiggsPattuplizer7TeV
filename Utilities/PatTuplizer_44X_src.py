######################################################################
#                                                                    #
# This version works with CMSSW_4_4_X                                #
#                                                                    #
######################################################################

import FWCore.ParameterSet.Config as cms

triggerProc = cms.string("REDIGI42X")

from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *

###################### Options #######################
isMC = bool(True)

isHorZZ = bool(True)

isZSkim = bool(False)

isSync = bool(False)

maxEvents = -1

MCGlobalTag = "START44_V13::All"
DataGlobalTag = "GR_R_44_V15C::All"

smearingRatio = 0.607 #0.0=HCP only, 1.0=2012D only, 0.607=whole dataset

doDebug = bool(False)
######################################################
    
### Info ###
if isMC :
    print "Sample Type: MC"
else :
    print "Sample Type: Data"

if isHorZZ :
    print "Sample is Higgs or ZZ --- No Skimming"
else :
    print "Sample is not Higgs or ZZ --- 3L Skim"

if isZSkim :
    print "******************************************"
    print "\n"
    print "WARNING:  YOU HAVE SET isZSkim TO TRUE!!!"
    print "THIS MEANS YOU WILL SKIM FOR A TIGHT Z"
    print "\n"
    print "******************************************"


###########


process = cms.Process("UFHiggsPatTuplizer")

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(maxEvents))


# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'ERROR' # Options: INFO, WARNING, ERROR
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.suppressWarning.append('patTrigger') # 1 warning per event on old runs!
process.MessageLogger.suppressWarning.append('classByHitsGlb') # kill stupid RPC hit associator warning
process.MessageLogger.cerr.FwkJob.limit=1
process.MessageLogger.cerr.ERROR = cms.untracked.PSet( limit = cms.untracked.int32(1))
                                                       
######################### Frontier Conditions #########################
# Conditions data for calibration and alignment                       #
# are defined in the Offline Conditions Database (ORCOF),             #
# which is read in CMSSW applications via Frontier caching servers.   #
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions  #
#######################################################################
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

if isMC:
    process.GlobalTag.globaltag=MCGlobalTag     
else:
    process.GlobalTag.globaltag=DataGlobalTag


#----------------------------------------------------------------------
# PAT
#----------------------------------------------------------------------
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerMatchEmbedder_cfi")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.helpers import *
from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.trigTools import *
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.muonTools import addMuonUserIsolation
from PhysicsTools.PatAlgos.tools.electronTools import addElectronUserIsolation
from PhysicsTools.PatAlgos.patEventContent_cff import *
from PhysicsTools.PatAlgos.recoLayer0.electronIsolation_cff import *
#from CommonTools.ParticleFlow.ParticleSelectors.pfCandsForIsolation_cff import *
from CommonTools.ParticleFlow.Isolation.tools_cfi import *
from CommonTools.ParticleFlow.pfNoPileUp_cff import *
from CommonTools.ParticleFlow.Tools.pfIsolation import *


process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('MYOUTPUTFILE'),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p_trilep') ),
                               outputCommands = cms.untracked.vstring('keep *')

                               )


##From mangano
process.pfMuId = cms.EDProducer('MuonPFIDMapProd',
                                muonLabel = cms.untracked.InputTag("muons"),
                                pfLabel = cms.untracked.InputTag("particleFlow"),
                                )


process.patMuons.userData.userInts.src = cms.VInputTag( cms.InputTag("pfMuId"), )


#----------------------------------------------------------------------
# PAT
#----------------------------------------------------------------------
if not isMC:
    removeMCMatching(process)
#'Photons', 'Electrons', 'Muons', 'Taus', 'Jets', 'METs', 'All', 'PFAll', 'PFElectrons', 'PFTaus', 'PFMuons'


#Isolation
from UFHiggsPattuplizer7TeV.CustomIso.addPatPfIso import * 
#addMuonUserIsolation(process)
#addElectronUserIsolation(process)
startPFIso(process)
#process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
#process.muIsoSequence = setupPFMuonIso(process, 'muons')
#process.phIsoSequence = setupPFPhotonIso(process,'photons')

#process.isoSequence = cms.Sequence(    process.pfParticleSelectionSequence *
#                                       process.eleIsoSequence *
#                                       process.muIsoSequence *
#                                       process.phIsoSequence
 #                                      )
                
#process.patPFCandidateIsoDepositSelection = cms.Sequence(
#    pfNoPileUpSequence *
#    ( pfAllNeutralHadrons +
#      pfAllChargedHadrons +
#      pfAllPhotons 
#      )
#    )




#MET
addPfMET(process, 'PF')


#To embed the Trigger Matching in the patMuons & patElectrons 
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
#process.patTrigger.processName = triggerProc
#process.patTrigger.onlyStandAlone = True 
process.patTrigger.processName  = '*' 
process.patTriggerEvent.processName = '*'

#----------------------------------------------------------------------
# pfPileUp/pfNoPileUp; move them outside of patDefaultSequence
#----------------------------------------------------------------------
process.pfPileUp.PFCandidates = "particleFlow"
process.pfNoPileUp.bottomCollection = "particleFlow"
process.patDefaultSequence.remove( process.pfPileUp )
process.patDefaultSequence.remove( process.pfNoPileUp )

#----------------------------------------------------------------------
# Rho correction and jets
#----------------------------------------------------------------------
# Needed to put AK5GenJets into the event
process.load("RecoJets/Configuration/RecoGenJets_cff")
process.load("RecoJets/Configuration/GenJetParticles_cff")
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJets = process.kt4PFJets.clone()
process.kt6PFJets.rParam = 0.6
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.doAreaFastjet = True

process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets','rho')
process.kt6PFJetsAA = process.kt6PFJets.clone( voronoiRfact = -0.9 )
process.kt6PFJetsForIso = process.kt6PFJets.clone( Rho_EtaMax = cms.double(2.5), Ghost_EtaMax = cms.double(2.5) )
process.kt6PFJetsForIsolation = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True, doAreaFastjet = True )#EGamma POG
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5) #EGamma POG
process.kt6PFJetsForIso3p1 = process.kt6PFJets.clone( Rho_EtaMax = cms.double(2.5), Ghost_EtaMax = cms.double(3.1) )
process.kt6PFJetsNoPU = process.kt6PFJets.clone( src = 'pfNoPileUp' )
process.kt6PFJetsNoPUAA = process.kt6PFJetsNoPU.clone( voronoiRfact = -0.9 )
process.kt6PFJetsForIsoNoPU = process.kt6PFJetsForIso.clone( src = 'pfNoPileUp' )
process.kt6PFJetsCentralNeutral = process.kt6PFJets.clone( src = cms.InputTag("pfAllNeutralHadronsAndPhotons"),
                                                           Ghost_EtaMax = cms.double(3.1), Rho_EtaMax = cms.double(2.5),
                                                           inputEtMin = cms.double(0.5) )


#process.patJetCorrFactorsFastJet = process.patJetCorrFactors.clone()
#process.patJetCorrFactorsFastJetAA = process.patJetCorrFactorsAA.clone()

# Re-cluster jets starting from pfNoPileUp
process.ak5PFJetsNoPU = process.ak5PFJets.clone( src = 'pfNoPileUp' )



def addFastJetCorrection(process,label,seq="patDefaultSequence",thisRho="kt6PFJets"):
    corrFact = getattr(process,"patJetCorrFactors"+label)
    setattr(process,"patJetCorrFactorsFastJet"+label,corrFact.clone())
    getattr(process,"patJetCorrFactorsFastJet"+label).levels[0] = 'L1FastJet'
    getattr(process,"patJetCorrFactorsFastJet"+label).rho = cms.InputTag(thisRho,"rho")
    getattr(process,"patJetCorrFactorsFastJet"+label).useRho = cms.bool(True)
    
    getattr(process,seq).replace(
        getattr(process,"patJetCorrFactors"+label),
        getattr(process,"patJetCorrFactors"+label) +
        getattr(process,"patJetCorrFactorsFastJet"+label) #+
 #       getattr(process,"patJetCorrFactorsFastJet"+label+"AA") 
        )
    
    getattr(process,"patJets"+label).jetCorrFactorsSource = cms.VInputTag(
#        cms.InputTag("patJetCorrFactorsFastJet"+label+"AA") ,
        cms.InputTag("patJetCorrFactorsFastJet"+label) ,
        cms.InputTag("patJetCorrFactors"+label) 
        )


# Jets &  B-Tagging
process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
process.load("RecoBTag.Configuration.RecoBTag_cff")
process.ak5PFJetTracksAssociatorAtVertex = process.ak5JetTracksAssociatorAtVertex.clone()
process.ak5PFJetTracksAssociatorAtVertex.jets = "ak5PFJets"

process.pfImpactParameterTagInfos = process.impactParameterTagInfos.clone()
process.pfImpactParameterTagInfos.jetTracks = "ak5PFJetTracksAssociatorAtVertex"
process.pfTrackCountingHighEffBJetTags = process.trackCountingHighEffBJetTags.clone()
process.pfTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("pfImpactParameterTagInfos") )
process.pfTrackCountingHighPurBJetTags = process.trackCountingHighPurBJetTags.clone()
process.pfTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("pfImpactParameterTagInfos") )
process.pfJetProbabilityBJetTags = process.jetProbabilityBJetTags.clone()
process.pfJetProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("pfImpactParameterTagInfos") )
process.pfJetBProbabilityBJetTags = process.jetBProbabilityBJetTags.clone()
process.pfJetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("pfImpactParameterTagInfos") )

process.pfSecondaryVertexTagInfos = process.secondaryVertexTagInfos.clone()
process.pfSecondaryVertexTagInfos.trackIPTagInfos = "pfImpactParameterTagInfos"
process.pfSimpleSecondaryVertexBJetTags = process.simpleSecondaryVertexBJetTags.clone()
process.pfSimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("pfSecondaryVertexTagInfos") )
process.pfCombinedSecondaryVertexBJetTags = process.combinedSecondaryVertexBJetTags.clone()
process.pfCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("pfImpactParameterTagInfos"), cms.InputTag("pfSecondaryVertexTagInfos") )
process.pfCombinedSecondaryVertexMVABJetTags = process.combinedSecondaryVertexMVABJetTags.clone()
process.pfCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("pfImpactParameterTagInfos"), cms.InputTag("pfSecondaryVertexTagInfos") )
process.jetsAndTagging = cms.Sequence(
                    process.ak5PFJetTracksAssociatorAtVertex     +
                    process.pfImpactParameterTagInfos            +
                    process.pfSecondaryVertexTagInfos            +
                    process.pfTrackCountingHighEffBJetTags       +
                    process.pfTrackCountingHighPurBJetTags       +
                    process.pfJetProbabilityBJetTags             +
                    process.pfJetBProbabilityBJetTags            +
                    process.pfSimpleSecondaryVertexBJetTags      +
                    process.pfCombinedSecondaryVertexBJetTags    +
                    process.pfCombinedSecondaryVertexMVABJetTags
                )

#---------------------------------------------------------------------
# Regression
#---------------------------------------------------------------------
process.load('EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi')
process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('patElectrons')
process.eleRegressionEnergy.energyRegressionType = cms.uint32(2)
process.eleRegressionEnergy.rhoCollection = cms.InputTag("kt6PFJetsForIso","rho")

#---------------------------------------------------------------
# preLeptonSequence: calibratedGSF + rho + jets from pfNoPileUp 
#---------------------------------------------------------------
process.load("EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi")
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedPatElectrons = cms.PSet(initialSeed = cms.untracked.uint32(1),
                                                                                     engineName = cms.untracked.string('TRandom3')
                                                                                     ),
                                                   calibratedOnlyPatElectrons = cms.PSet(initialSeed = cms.untracked.uint32(1),
                                                                                         engineName = cms.untracked.string('TRandom3')
                                                                                         ),
                                                   regressedOnlyPatElectrons = cms.PSet(initialSeed = cms.untracked.uint32(1),
                                                                                        engineName = cms.untracked.string('TRandom3')
                                                                                        ),
                                                   regressedNewPatElectrons = cms.PSet(initialSeed = cms.untracked.uint32(1),
                                                                                        engineName = cms.untracked.string('TRandom3')
                                                                                        ),
                                                   )

#process.calibratedPatElectrons.debug = cms.bool(True)
process.calibratedPatElectrons.isMC = cms.bool(isMC)
process.calibratedPatElectrons.synchronization = cms.bool(isSync)
process.calibratedPatElectrons.isAOD = cms.bool(True)
process.calibratedPatElectrons.inputPatElectronsTag = cms.InputTag("eleRegressionEnergy")
process.calibratedPatElectrons.updateEnergyError = cms.bool(True)
process.calibratedPatElectrons.correctionsType = cms.int32(2)# 0-No Corrections, 1-Scale+old regression, 2-Scale+new regression, 3-Scale+no regression
process.calibratedPatElectrons.combinationType = cms.int32(3)# 0-No comb, 1-Standard comb w/smeared+corr regression, 2-std comb unsmeared+uncorr regression, 3-new regression comb smeared+corr regression 
process.calibratedPatElectrons.lumiRatio = cms.double(smearingRatio)#0.0=HCP only, 1.0=2012D only, 0.607=whole dataset



if isMC:
    process.calibratedPatElectrons.inputDataset = cms.string("Fall11") 
else:
    process.calibratedPatElectrons.inputDataset = cms.string("Jan16ReReco")


process.calibratedOnlyPatElectrons = process.calibratedPatElectrons.clone(inputPatElectronsTag = cms.InputTag("patElectrons"),
                                                                          correctionsType = cms.int32(1),
                                                                          combinationType = cms.int32(0),
                                                                          lumiRatio = cms.double(smearingRatio)
                                                                          )

process.regressedOnlyPatElectrons = process.calibratedPatElectrons.clone(inputPatElectronsTag = cms.InputTag("eleRegressionEnergy"),
                                                                         correctionsType = cms.int32(0),
                                                                         combinationType = cms.int32(1),
                                                                         lumiRatio = cms.double(smearingRatio)
                                                                         )

process.regressedNewPatElectrons = process.calibratedPatElectrons.clone(inputPatElectronsTag = cms.InputTag("eleRegressionEnergy"),
                                                                         correctionsType = cms.int32(0),
                                                                         combinationType = cms.int32(3),
                                                                         lumiRatio = cms.double(smearingRatio)
                                                                         )




process.preLeptonSequence = cms.Sequence(
    process.ak5PFJets * 
    process.kt6PFJets *
    process.kt6PFJetsAA *
    process.kt6PFJetsForIso *
    process.kt6PFJetsForIsolation *
#    process.kt6PFJetsCentralNeutral *
    process.kt6PFJetsForIso3p1 *
    process.pfPileUp *
    process.pfNoPileUp * (process.ak5PFJetsNoPU *
                          process.kt6PFJetsNoPU * 
                          process.kt6PFJetsNoPUAA * 
                          process.kt6PFJetsForIsoNoPU ) 
    
    )




#----------------------------------------------------------------------
# PHOTONS PATH
#----------------------------------------------------------------------
process.load("PhysicsTools.PatAlgos.producersLayer1.photonProducer_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.photonSelector_cfi")

#process.cleanElePhotonOverlap = cms.EDProducer("PATPhotonCleaner",
#    ## Input collection of Photons
#    src = cms.InputTag("patPhotons"),
#
#    # preselection (any string-based cut for pat::Photon)
#    preselection = cms.string(''),
#
#    # overlap checking configurables
#    checkOverlaps = cms.PSet(
#        electrons = cms.PSet(
#           src       = cms.InputTag("cleanPatElectrons"),
#           algorithm = cms.string("bySuperClusterSeed"),
#           requireNoOverlaps = cms.bool(True), # mark photons that overlap with electrons
#                                                # for further studies, but DO NOT discard
#                                                # them
#        ),
#    ),
#
#    # finalCut (any string-based cut for pat::Photon)
#    finalCut = cms.string(''),
#
#)


#----------------------------------------------------------------------
# ELECTRONS
#----------------------------------------------------------------------
#process.patElectrons.embedPFCandidate = True
#process.patElectrons.embedSuperCluster = True
#process.patElectrons.embedTrack = True
#process.patElectrons.addElectronID = True
#process.electronMatch.matched = "prunedGen"

####### PF & DetBased Iso
#isoKinds4Electrons='DetBased+PFlow'
###OLD CUSTOM
addElectronPFIso(process,process.patElectrons)



####### MVA Electron ID
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )


######Code for the CIC eID
# Modules for the Cut-based Electron ID in the VBTF prescription
import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi as vbtfid
process.eidVBTFRel95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95relIso' )
process.eidVBTFRel80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80relIso' )
process.eidVBTFCom95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95cIso'   )
process.eidVBTFCom80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80cIso'   )

process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi")
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi")

import RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi as cic
process.eidCiCVeryLoose = cic.eidVeryLoose.clone()
process.eidCiCLoose = cic.eidLoose.clone()
process.eidCiCMedium = cic.eidMedium.clone()
process.eidCiCTight = cic.eidTight.clone()

import RecoEgamma.ElectronIdentification.cutsInCategoriesHZZElectronIdentificationV06_cfi as cicHZZ
process.eidCiCHZZVeryLoose = cicHZZ.eidHZZVeryLoose.clone()
process.eidCiCHZZLoose = cicHZZ.eidHZZLoose.clone()
process.eidCiCHZZMedium = cicHZZ.eidHZZMedium.clone()
process.eidCiCHZZTight = cicHZZ.eidHZZTight.clone()
process.eidCiCHZZSuperTight = cicHZZ.eidHZZSuperTight.clone()
process.eidCiCHZZHyperTight1 = cicHZZ.eidHZZHyperTight1.clone()


process.eidSequence = cms.Sequence(
            process.eidVBTFRel95  +
            process.eidVBTFRel80  +
            process.eidVBTFCom95  +
            process.eidVBTFCom80  +
            process.eidVeryLoose  +
            process.eidLoose      +
            process.eidMedium     +
            process.eidTight      +
            process.eidSuperTight +
            process.eidHyperTight1 +
            process.eidCiCVeryLoose +
            process.eidCiCLoose +
            process.eidCiCMedium +
            process.eidCiCTight +
            process.eidCiCHZZVeryLoose +
            process.eidCiCHZZLoose +
            process.eidCiCHZZMedium +
            process.eidCiCHZZTight +
            process.eidCiCHZZSuperTight +
            process.eidCiCHZZHyperTight1 +
            process.mvaID
            
        )


# Electron Selection
process.patElectrons.electronIDSources = cms.PSet(
            eidVBTFRel95   = cms.InputTag("eidVBTFRel95"),
            eidVBTFRel80   = cms.InputTag("eidVBTFRel80"),
            eidVBTFCom95   = cms.InputTag("eidVBTFCom95"),
            eidVBTFCom80   = cms.InputTag("eidVBTFCom80"),
            eidVeryLoose   = cms.InputTag("eidVeryLoose"),
            eidLoose       = cms.InputTag("eidLoose"),
            eidMedium      = cms.InputTag("eidMedium"),
            eidTight       = cms.InputTag("eidTight"),
            eidSuperTight  = cms.InputTag("eidSuperTight"),
            eidCiCVeryLoose = cms.InputTag("eidCiCVeryLoose"),
            eidCiCLoose = cms.InputTag("eidCiCLoose"),
            eidCiCMedium = cms.InputTag("eidCiCMedium"),
            eidCiCTight = cms.InputTag("eidCiCTight"),
            eidCiCHZZVeryLoose = cms.InputTag("eidCiCHZZVeryLoose"),
            eidCiCHZZLoose = cms.InputTag("eidCiCHZZLoose"),
            eidCiCHZZMedium = cms.InputTag("eidCiCHZZMedium"),
            eidCiCHZZTight = cms.InputTag("eidCiCHZZTight"),
            eidCiCHZZSuperTight = cms.InputTag("eidCiCHZZSuperTight"),
            eidCiCHZZHyperTight1 = cms.InputTag("eidCiCHZZHyperTight1"),
            mvaTrigV0 = cms.InputTag("mvaTrigV0"),
            mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
            )


######Code for the conversion
#process.convValueMapProd = cms.EDProducer('ConvValueMapProd',
#    gsfLabel = cms.untracked.InputTag("gsfElectrons"),
#    tkLabel = cms.untracked.InputTag("generalTracks")
#)



######Code for trigger matching
#process.eleTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRDPtLessByR",
#    src     = cms.InputTag( "patElectrons" ),
#    matched = cms.InputTag( "patTrigger" ),
#    andOr          = cms.bool( False ),
#    matchedCuts = cms.string('coll("hltL1IsoRecoEcalCandidate")||coll("hltL1NonIsoRecoEcalCandidate")'),
#    maxDPtRel = cms.double( 0.5 ),
#    maxDeltaR = cms.double( 0.5 ),
#    resolveAmbiguities    = cms.bool( True ),
#    resolveByMatchQuality = cms.bool( True )
#)


########Code for Cross Clean --- cone size 0.01
#process.cleanPatElectrons = cms.EDProducer("PATElectronCleaner",
#    ## pat electron input source
#    src = cms.InputTag("patElectrons"), 
#    # preselection (any string-based cut for pat::Electron)
#    preselection = cms.string(''),
#    # overlap checking configurables
#    checkOverlaps = cms.PSet(
#        muons = cms.PSet(
#           src       = cms.InputTag("patMuons"),
#           algorithm = cms.string("byDeltaR"),
#           preselection        = cms.string("isGlobalMuon"),  # don't preselect the muons
#           deltaR              = cms.double(0.01),  
#           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
#           pairCut             = cms.string(""),
 #          requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
#        )
#    ),
#    # finalCut (any string-based cut for pat::Electron)
#    finalCut = cms.string(''),
#)

# The electron selector
process.softElectrons = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("patElectrons"),
    cut = cms.string("pt>5"),
)


process.selectedPatElectrons.cut = (
    "abs(eta) < 2.5 "
    #    "& (isEE || isEB) & !isEBEEGap " 
    #    "& (electronID('eidVBTFRel95') == 7 || electronID('eidVBTFRel95') == 5 " +
    #    "|| electronID('eidVBTFCom95') == 7 || electronID('eidVBTFCom95') == 5 " +
    #    "|| electronID('eidVBTFRel80') == 7 || electronID('eidVBTFRel80') == 5 " +
    #    "|| electronID('eidVBTFCom80') == 7 || electronID('eidVBTFCom80') == 5 " +
    #    "|| electronID('eidLoose') == 7 || electronID('eidLoose') == 5 " +
    #    "|| electronID('eidLoose') == 13 || electronID('eidLoose') == 15 " +
    #    "|| electronID('eidVeryLoose') == 7 || electronID('eidVeryLoose') == 5 " +
    #    "|| electronID('eidVeryLoose') == 13 || electronID('eidVeryLoose') == 15 " +
    #    "|| electronID('eidTight') == 7 || electronID('eidTight') == 5 " +
    #    "|| electronID('eidTight') == 13 || electronID('eidTight') == 15  " +
    #    "|| electronID('eidSuperTight') == 7 || electronID('eidSuperTight') == 5 " +
    #    "|| electronID('eidSuperTight') == 13 || electronID('eidSuperTight') == 15) "
    )


process.preElectronSequence = cms.Sequence(
    process.eidSequence +
    process.patTrigger 
)


process.patDefaultSequence.replace(
    process.patElectrons,
    process.patElectrons * process.softElectrons
    )



# ----------------------------------------------------------------------
# MUON PATH  
# ----------------------------------------------------------------------
process.patMuons.embedPFCandidate = True
process.patMuons.embedTrack = True
#process.muonMatch.matched = "prunedGen"

####### PF & DetBased Iso
#isoKinds4Muons = 'DetBased+PFlow'
#skipVetoLeptonSeqence=True
#skipSelectingPfCands=False
### OLD CUSTOM
addMuonPFIso(process,process.patMuons)


if isMC: 
    if False: ## Turn this on to get extra info on muon MC origin, on GEN-SIM-RECO
        process.load("MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi")
        from MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi import addUserData as addClassByHits
        addClassByHits(process.patMuons, labels=['classByHitsGlb'], extraInfo = True)
        process.muonClassificationByHits = cms.Sequence(process.mix * process.trackingParticlesNoSimHits * process.classByHitsGlb )
        
        process.preMuonSequence = cms.Sequence(process.muonMatch +
                                               process.patTrigger +
                                               process.muonClassificationByHits +
                                               process.MuonIsolationMakeIso)

    else:
        
        process.preMuonSequence = cms.Sequence(process.muonMatch +
                                               process.patTrigger #+
                                               #process.MuonIsolationMakeIso
                                               )
else:
    process.preMuonSequence = cms.Sequence(process.patTrigger #+
                                           #process.MuonIsolationMakeIso
                                           )



# The Muon Selector
process.softMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("(isGlobalMuon || isTrackerMuon) && pt>3")
)

process.selectedPatMuons.cut = (
            "(isGlobalMuon || isTrackerMuon) & abs(eta) < 2.5" 
#            "innerTrack().hitPattern().numberOfValidTrackerHits > 10 && "                      +
#            "innerTrack().hitPattern().numberOfValidPixelHits > 0 && "                         +
#            "globalTrack().hitPattern().numberOfValidMuonHits > 0 && "                         
#            "dB < 0.2 && "                                                                     +
#            "trackIso + caloIso < 0.15 * pt && "                                               +
#            "numberOfMatches > 0 && abs(eta) < 2.4"
        )


#---> USING THE PAT DEFAULT SEQUENCE
process.patDefaultSequence.replace(
    process.patMuons,
    process.patMuons * process.softMuons * process.selectedPatMuons
    )

##clean muons by segments --- from Gio
process.boostedMuons = cms.EDProducer("PATMuonCleanerBySegments",
                                      src = cms.InputTag("cleanPatMuons"),
                                      preselection = cms.string("track.isNonnull"),
                                      passthrough  = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                      fractionOfSharedSegments = cms.double(0.499),
                                      )


process.rochesterMuons = cms.EDProducer('RochesterPATMuonCorrector',
                                        src = cms.InputTag("boostedMuons"),
                                        )


muscleIdentifier = ''
if isMC:
    muscleIdentifier = 'Fall11_START42'
else:
    muscleIdentifier = 'Data2011_42X'
    
process.muscleMuons = cms.EDProducer('MuScleFitPATMuonCorrector',
                                     src = cms.InputTag("boostedMuons"),
                                     debug = cms.bool(doDebug),
                                     identifier = cms.string(muscleIdentifier),
                                     applySmearing = cms.bool(isMC),
                                     fakeSmearing = cms.bool(isSync)
                                     )
    
    



#----------------------------------------------------------------------
# Jet sequences
#----------------------------------------------------------------------

# Jet energy corrections to use:
##MC
myCorrLabels = cms.vstring()

if isMC:
    myCorrLabels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute')
else:
    myCorrLabels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual')
                

switchJetCollection(
    process,
    cms.InputTag('ak5PFJets'),
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF',myCorrLabels),
    doType1MET   = True,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True
    )

addJetCollection(
    process,
    cms.InputTag("ak5PFJetsNoPU"),
    algoLabel    = "NoPU",
    typeLabel    = "",
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF',myCorrLabels),
    doL1Cleaning = False,
    doL1Counters = True,                 
    doType1MET   = True,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True,
    jetIdLabel   = 'ak5',
)


process.patJets.tagInfoSources  = cms.VInputTag(
        cms.InputTag("secondaryVertexTagInfosAOD"),
            )


# Some stuff to save space
#process.patJets.embedCaloTowers = False
#process.patJetsNoPU.embedCaloTowers = False
#process.patJets.addTagInfos = False
#process.patJetsNoPU.addTagInfos = False
#process.patJets.embedPFCandidates = False
#process.patJetsNoPU.embedPFCandidates = False
#process.patJets.addAssociatedTracks = False
#process.patJetsNoPU.addAssociatedTracks = False

#Jet Selection
process.selectedPatJets.cut = (
            " et > 5 && abs(eta) < 5.0"
            )


# Not set up correctly by PAT:
process.cleanPatJetsNoPU = process.cleanPatJets.clone( src = cms.InputTag("selectedPatJetsNoPU") )
process.patDefaultSequence.replace(
    process.cleanPatJets,
    process.cleanPatJets +
    process.cleanPatJetsNoPU 
)

# Add the fast jet correction:
addFastJetCorrection(process,"")
#addFastJetCorrection(process,"NoPU","patDefaultSequence","kt6PFJetsNoPU")


process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
                                            src = cms.InputTag('offlinePrimaryVertices'),
                                            cut = cms.string('!isFake && isValid && ndof >= 4.0 && position.Rho < 2.0 && abs(z) < 24'),
                                            filter = cms.bool(True)
                                        )


#### PU Jet ID for VBF Selection ###
process.load("CMGTools.External.pujetidsequence_cff")
process.puJetId.vertexes = cms.InputTag("goodOfflinePrimaryVertices")
process.puJetMva.vertexes = cms.InputTag("goodOfflinePrimaryVertices")


process.source = cms.Source ("PoolSource",
                             # Disable duplicate event check mode because the run and event -numbers
                             # are incorrect in current Madgraph samples (Dec 16, 2008)
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                             
                             fileNames = cms.untracked.vstring(),      
                             )

if isMC:
    process.source.fileNames = cms.untracked.vstring(
        #"/store/mc/Fall11/GluGluToHToZZTo4L_M-120_7TeV-powheg-pythia6/AODSIM/PU_S6_START42_V14B-v1/0000/40E86BD8-0BF0-E011-BA16-00215E21D5C4.root"

        )
else:
    process.source.fileNames = cms.untracked.vstring(
 
           

        )
    
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False)
                                     #,SkipEvent = cms.untracked.vstring('ProductNotFound')
                                     )



process.out.outputCommands.extend(cms.untracked.vstring('drop *'))
process.out.outputCommands.extend(patEventContent)
#process.out.outputCommands.extend(patExtraAodEventContent)
process.out.outputCommands.extend(patTriggerL1RefsEventContent)
process.out.outputCommands.extend(patEventContentTriggerMatch)
process.out.outputCommands.extend(cms.untracked.vstring(
    #####DROP#####
    'drop patTaus_*_*_*',
    'drop patMETs_patMETs_*_*',
    ##### doubles #####
    'keep *_*_sigmas_*',
    'keep *_*_rhos_*',
    'keep *_*_sigma_*',
    'keep *_*_rho_*',
#    'keep *_kt6PF*_rho_'+process.name_(),
#    'keep double_*PFlow*_*_PAT',
    'keep doubleedmValueMap_*_*_*',
    'keep recoIsoDepositedmValueMap_iso*_*_*',
#    'drop double_*_*_RECO',
    ##### vertices #####
    'keep *_goodOfflinePrimaryVertices*_*_*',
    'keep *_offlinePrimaryVertices*_*_*',
    'keep *_offlineBeamSpot_*_*',
    ##### Trigger ######
    'keep *_patTrigger_*_*',
    'keep *_hltTriggerSummaryAOD_*_*',
    'keep *_TriggerResults_*_*',
    'keep *_vertexMapProd_*_*',
    'keep l1extraL1EmParticles_*_*_*',
    'keep l1extraL1MuonParticles_*_*_*',
    'keep L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap_*_*',
    ##### counter ######
    'keep edmMergeableCounter_*_*_*',
    ##### leptons ######
    'keep *_muons_*_*',
    'keep *_gsfElectrons_*_*',
    'keep *_cleanPatMuons_*_*',
    'keep *_calibratedPatElectrons_*_*',
    'keep *_calibratedOnlyPatElectrons_*_*',
    'keep *_regressedOnlyPatElectrons_*_*',
    'keep *_regressedNewPatElectrons_*_*',
    'keep *_boostedMuons_*_*',
    'keep *_rochesterMuons_*_*',
    'keep *_muscleMuons_*_*',
    'keep *_cleanPatElectrons_*_*',
    'keep *_cleanPatPhotons_*_*',
    ##### PFlow #####
    'keep recoTracks_generalTracks_*_*',
    'keep *_*fsrPhoton*_*_*',
    'keep *_*boostedFsrPhoton*_*_*',
    'keep *_particleFlow_*_*',
    'keep *_pfMET_*_*',
    'keep patMETs_patMETsPF*_*_*',
    ##### Jets #####
    'keep patJets_selectedPatJets_*_*',
    'keep patJets_selectedPatJetsPFlow_*_*',
    "keep *_puJetId_*_*", # input variables
    "keep *_puJetMva_*_*" # final MVAs and working point flags        
    )
                                  )


if isMC:
    process.out.outputCommands.extend(cms.untracked.vstring(
        'keep *_*PileupInfo_*_*',
        'keep LHERunInfoProduct_source_*_*',
        'keep LHEEventProduct_source_*_*',
        'keep recoGenParticles_genParticles*_*_*',
        'keep GenEventInfoProduct_*_*_*',
        'keep GenRunInfoProduct_*_*_*'
        
        
        )
                                     )



### Physics Declared Filter (for data)
process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

### Technical Trigger Filters (for data)
#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
#process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
#process.L1T1=process.hltLevel1GTSeed.clone()
#process.L1T1.L1TechTriggerSeeding = cms.bool(True)
#process.L1T1.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
#process.L1T1.L1SeedsLogicalExpression = cms.string('0')
#process.bscnobeamhalo = cms.Path(process.L1T1)

### No scraping
process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(doDebug),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.2)
                                  )


# HB + HE noise filtering
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')
process.HBHENoiseFilter.minIsolatedNoiseSumE        = 999999.
process.HBHENoiseFilter.minNumIsolatedNoiseChannels = 999999
process.HBHENoiseFilter.minIsolatedNoiseSumEt       = 999999.


############################## S K I M M I N G #############################

######### Tri-lepton #############
muonsInDataString = "rochesterMuons"
muonInputTag = cms.InputTag(muonsInDataString)
elecInDataString = "calibratedPatElectrons"
electronInputTag = cms.InputTag(elecInDataString)

## Declaring cut objects
process.goodMuonsLowPt = cms.EDFilter("CandViewSelector",
                                       src = muonInputTag,
                                       cut = cms.string("pt > 3.0 & abs(eta) < 2.4 & (isGlobalMuon = 1 || isTrackerMuon = 1)")
                                       )

process.goodElectronsLowPt = cms.EDFilter("CandViewSelector",
                                           src = electronInputTag,
                                           cut = cms.string("pt > 4.0 & abs(eta) < 2.5")
                                           )

process.goodMuonsHighPt = cms.EDFilter("CandViewSelector",
                                        src = cms.InputTag("goodMuonsLowPt"),
                                        cut = cms.string("pt > 7.0 & abs(eta) < 2.4 & (isGlobalMuon = 1 || isTrackerMuon = 1)")
                                        )

process.goodElectronsHighPt = cms.EDFilter("CandViewSelector",
                                            src = cms.InputTag("goodElectronsLowPt"),
                                            cut = cms.string("pt > 7.0 & abs(eta) < 2.5")
                                            )

process.goodLeptonsLowPt = cms.EDProducer("CandViewMerger",
                                          src = cms.VInputTag( "goodElectronsLowPt", "goodMuonsLowPt"),
                                          )

process.goodLeptonsHighPt = cms.EDProducer("CandViewMerger",
                                           src = cms.VInputTag( "goodElectronsHighPt", "goodMuonsHighPt"),
                                           )

process.leptonFilterLowPt = cms.EDFilter("CandViewCountFilter",
                                    src = cms.InputTag("goodLeptonsLowPt"),
                                    minNumber = cms.uint32(3)
                                    )

process.leptonFilterHighPt = cms.EDFilter("CandViewCountFilter",
                                         src = cms.InputTag("goodLeptonsHighPt"),
                                         minNumber = cms.uint32(2)
                                         )



############# Z Candle ###############
# cuts
MUON_CUT=("pt > 10 && abs(eta)<2.5 && (isGlobalMuon || isTrackerMuon)")
ELECTRON_CUT=("pt > 10 && abs(eta)<2.5")
DIMUON_CUT=("mass > 40")
DIELECTRON_CUT=("mass > 40")
EMU_CUT=("mass > 40")

# single lepton selectors
process.goodHzzMuons = cms.EDFilter("MuonRefSelector",
                            src = cms.InputTag("muons"),
                            cut = cms.string(MUON_CUT)
                            )

process.goodHzzElectrons = cms.EDFilter("GsfElectronRefSelector",
                                src = cms.InputTag("gsfElectrons"),
                                cut = cms.string(ELECTRON_CUT)
                                )

# dilepton selectors
process.diHzzMuons = cms.EDProducer("CandViewShallowCloneCombiner",
                            decay       = cms.string("goodHzzMuons goodHzzMuons"),
                            checkCharge = cms.bool(False),
                            cut         = cms.string(DIMUON_CUT)
                            )

process.diHzzElectrons = cms.EDProducer("CandViewShallowCloneCombiner",
                                decay       = cms.string("goodHzzElectrons goodHzzElectrons"),
                                checkCharge = cms.bool(False),
                                cut         = cms.string(DIELECTRON_CUT)
                                )

process.crossHzzLeptons  = cms.EDProducer("CandViewShallowCloneCombiner",
                                  decay       = cms.string("goodHzzMuons goodHzzElectrons"),
                                  checkCharge = cms.bool(False),
                                  cut         = cms.string(EMU_CUT)
                                  )

process.dilep = cms.EDProducer("CandViewMerger",
                               src = cms.VInputTag(cms.InputTag("diHzzMuons"), cms.InputTag("diHzzElectrons"),cms.InputTag("crossHzzLeptons"),)
                               )


process.ZSkim = cms.EDFilter("CandViewCountFilter",
                              src = cms.InputTag("dilep"),
                              minNumber = cms.uint32(1)
                              )


process.ZCandleSkim = cms.Sequence( process.goodHzzMuons * process.goodHzzElectrons * process.diHzzMuons * process.diHzzElectrons * process.crossHzzLeptons * process.dilep * process.ZSkim )


process.ZCandleHistosM = cms.EDAnalyzer("CandViewHistoAnalyzer",
                                            src = cms.InputTag("diHzzMuons"),
                                            histograms = cms.VPSet(
                    cms.PSet(
                            name = cms.untracked.string("mass"),
                                                    description = cms.untracked.string("Z Candidate Mass"),
                                                    nbins = cms.untracked.int32(100),
                                                    min = cms.untracked.double(35.0),
                                                    max = cms.untracked.double(120.0),
                                                    plotquantity = cms.untracked.string("mass")
                                            ),
                    
                    )
)
process.ZCandleHistosE = cms.EDAnalyzer("CandViewHistoAnalyzer",
                                        src = cms.InputTag("diHzzElectrons"),
                                        histograms = cms.VPSet(
    cms.PSet(
    name = cms.untracked.string("mass"),
    description = cms.untracked.string("Z Candidate Mass"),
    nbins = cms.untracked.int32(100),
    min = cms.untracked.double(35.0),
    max = cms.untracked.double(120.0),
    plotquantity = cms.untracked.string("mass")
    ),
    
    )
)

process.ZCandleHistosEM = cms.EDAnalyzer("CandViewHistoAnalyzer",
                                        src = cms.InputTag("crossHzzLeptons"),
                                        histograms = cms.VPSet(
    cms.PSet(
    name = cms.untracked.string("mass"),
    description = cms.untracked.string("Z Candidate Mass"),
    nbins = cms.untracked.int32(100),
    min = cms.untracked.double(35.0),
    max = cms.untracked.double(120.0),
    plotquantity = cms.untracked.string("mass")
    ),
    
    )
)


################# 2L SKim ###################
process.muons4skim = cms.EDFilter("CandViewSelector",
                                  src = cms.InputTag("rochesterMuons"),
                                  cut = cms.string("(isTrackerMuon||isGlobalMuon) && abs(eta) < 2.4 && pt > 10"),
                                  )

process.electrons4skim = cms.EDFilter("CandViewSelector",
                                      src = cms.InputTag("calibratedPatElectrons"),
                                      cut = cms.string("abs(eta) < 2.5 && pt > 10"),
                                      )

process.leptons4skim = cms.EDProducer("CandViewMerger",
                                      src = cms.VInputTag( cms.InputTag("muons4skim"),
                                                           cms.InputTag("electrons4skim"), )
                                      )

process.dileptons4skim = cms.EDProducer("CandViewShallowCloneCombiner",
                                        decay = cms.string('leptons4skim leptons4skim'),
                                        cut = cms.string('deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) > 0.01'),
                                        checkCharge = cms.bool(False)
                                        )

process.skim2010 = cms.EDFilter("CandViewSelector",
                               src = cms.InputTag("dileptons4skim"),
                               cut = cms.string('min(daughter(0).pt,daughter(1).pt) > 10 && max(daughter(0).pt,daughter(1).pt) > 20'),
                               filter = cms.bool(True),
                               )

process.skimCharge = cms.EDFilter("CandViewSelector",
                                  src = cms.InputTag("dileptons4skim"),
                                  cut = cms.string('daughter(0).charge * daughter(1).charge == -1'),
                                  filter = cms.bool(True),
                                  )



process.skim40NoOF  = cms.EDFilter("CandViewSelector",
                                   src = cms.InputTag("dileptons4skim"),
                                   cut = cms.string('mass > 40 && abs(daughter(0).pdgId) == abs(daughter(1).pdgId)'), ## and SF only
                                   filter = cms.bool(True),
                                   )

process.skim2L = cms.Sequence(  process.muons4skim
                                + process.electrons4skim
                                + process.leptons4skim
                                + process.dileptons4skim
                                + process.skim2010
                                + process.skimCharge
                                + process.skim40NoOF
                                )


#process.source.eventsToProcess = cms.untracked.VEventRange("1:54066")

#EventCount                                                                                                                             
process.nEventsTotal = cms.EDProducer("EventCountProducer")
process.nEventsTriLep = cms.EDProducer("EventCountProducer")
process.nEvents2l2nu = cms.EDProducer("EventCountProducer")
process.nEvents2LSkim = cms.EDProducer("EventCountProducer")

process.load('UFHiggsPattuplizer7TeV.FSRPhotons.fsrPhotons_cff')

process.p_trilep = cms.Path(
    process.nEventsTotal
    *process.goodOfflinePrimaryVertices
    *process.isoSequence
    *process.pfMuId
    *process.preLeptonSequence
    *(process.preMuonSequence * process.preElectronSequence)
    *process.noscraping
    *process.patDefaultSequence
    *process.selectedPatJets 
    ##   VBF Jets
    *process.puJetIdSqeuence
    ##   FSR
    *process.fsrPhotonSequence
    ###  Regression
    *process.eleRegressionEnergy
    *process.calibratedPatElectrons
    *process.calibratedOnlyPatElectrons
    *process.regressedOnlyPatElectrons
    *process.regressedNewPatElectrons
    *process.boostedMuons
    *process.rochesterMuons
    *process.muscleMuons
    ###  Skim
    *process.goodMuonsLowPt
    *process.goodElectronsLowPt
    *process.goodMuonsHighPt
    *process.goodElectronsHighPt
    *process.goodLeptonsLowPt
    *process.goodLeptonsHighPt
    *process.leptonFilterLowPt
    *process.leptonFilterHighPt
    *process.nEventsTriLep
    )


#print process.p_trilep

if isMC:
    process.load("GeneratorInterface.GenFilters.TotalKinematicsFilter_cfi")
    process.totalKinematicsFilter.tolerance=5.0
    process.p_trilep.replace(process.goodOfflinePrimaryVertices,
                              process.goodOfflinePrimaryVertices *
                              process.totalKinematicsFilter)
    process.p_trilep.replace(process.ak5PFJets,
                             process.genParticlesForJets
                             *process.ak5GenJets
                             *process.ak5PFJets)

                             
if isHorZZ:
    process.p_trilep.remove(process.leptonFilterLowPt)
    process.p_trilep.remove(process.leptonFilterHighPt)


if isZSkim :
    process.p_trilep.remove(process.leptonFilterLowPt)
    process.p_trilep.remove(process.leptonFilterHighPt)
    process.p_trilep.remove(process.goodMuonsLowPt)
    process.p_trilep.remove(process.goodElectronsLowPt)
    process.p_trilep.remove(process.goodMuonsHighPt)
    process.p_trilep.remove(process.goodElectronsHighPt)
    process.p_trilep.remove(process.goodLeptonsLowPt)
    process.p_trilep.remove(process.goodLeptonsHighPt)
    process.p_trilep.replace(process.nEventsTriLep,
                             process.skim2L
                             *process.nEvents2LSkim)


#Turn pat::TriggerEvent on                                                                                                                          
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, hltProcess = '*' )
switchOnTriggerMatching( process, hltProcess = '*' )
process.patTrigger.addL1Algos = cms.bool( True )

process.out_step = cms.EndPath(process.out)

process.schedule = cms.Schedule(
    process.p_trilep,
    process.out_step
    )


