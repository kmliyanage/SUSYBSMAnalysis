
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'ana_datamc_Wantitop'
config.General.workArea = 'crab_2017_MC'
#config.General.transferLogs = 'True'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_crab.py'   
#config.JobType.priority = 1
config.Data.inputDataset =  '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
config.Data.inputDBS = 'global'

#config.Data.splitting = 'EventAwareLumiBased'
config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 20
    
config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.outputDatasetTag = 'ana_datamc_Wantitop'
config.Data.outLFNDirBase = '/store/user/kaliyana'
config.Data.ignoreLocality = True 
config.Site.whitelist = ["T2_IT_Bari"]
config.Site.storageSite = 'T2_IT_Bari'
#config.Debug.scheddName = 'crab3-5@vocms059.cern.ch'
