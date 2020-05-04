# SUSYBSMAnalysis

# For MC 

Files needed: 
	/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/histos.py
	/SUSYBSMAnalysis/Zprime2muAnalysis/python/MCSamples.py
	/SUSYBSMAnalysis/Zprime2muAnalysis/python/ METFilterMiniAOD_cfi.py (only to change between 	data and MC)

Changes to be done:

	histos.py

•	process.source.fileNames =[#'file:./pat.root'
		<das root file> (only if needed to run in local)	]

•	process.GlobalTag.globaltag =<Golab tag of the MC sample>

•	samples = [('dy50to120', '/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'), ……………..  ]

•	def ntuplify(process, fill_gen_info=True)

•	TriggerResults_src = cms.InputTag('TriggerResults', '', 'PAT'),	#mc

•	config.General.workArea = 'crab_<sample name>'

•	config.Data.outLFNDirBase = '/store/user/kaliyana'

•	config.Site.whitelist = ["T2_IT_Bari"]

•	config.Site.storageSite = 'T2_IT_Bari'

	MCSamples.py

•	samples = [
	sample('dy50to120',   'DY50to120', '/ZToMuMu_NNPDF31_13TeV-	powheg_M_50_120/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-	v2/MINIAODSIM', 2982000, 209 , 1., 2112.904,   k_factor=1.),

			

To submit:	 python histos.py submit no_data

# For Data

Files needed: 
	/SUSYBSMAnalysis/Zprime2muAnalysis/test/DataMCSpectraComparison/histos.py
	/SUSYBSMAnalysis/Zprime2muAnalysis/python/METFilterMiniAOD_cfi.py (only to change between 	data and MC)
	/SUSYBSMAnalysis/Zprime2muAnalysis/python/HistosFromPAT_cfi.py (fill_gen_info = cms.bool(False) # for DATA)
	/SUSYBSMAnalysis/Zprime2muAnalysis/python/goodlumis.py(to add JSON files)

        histos.py

•	process.source.fileNames =[#'file:./pat.root'
		<das root file> (only if needed to run in local)	]

•	process.GlobalTag.globaltag =<Golab tag of the dataset>

•	def ntuplify(process, fill_gen_info=False)

•	TriggerResults_src = cms.InputTag('TriggerResults', '', 'RECO'),	#data

•	dataset_details = [
						('SingleMuonRun2018A_v1-PromptReco', '/SingleMuon/Run2018A-PromptReco-v1/MINIAOD'),


To submit:	 python histos.py submit no_mc

