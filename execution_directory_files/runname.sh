mv ALLUVIAL.DAT $1.ALV
mv OUTELEV.DAT	$1.DAT
mv SUMMARY.DAT	 $1.SUM
mv GRADIENT.DAT $1.GRA
mv BISTABLE.DAT $1.BIS
mv RESIST.OUT $1.RES
mv REGOLITH.DAT $1.REG
mv OUTBASE.DAT	$1.BSE
mv DEFORM.DAT $1.DFM
mv BEDROCK.DAT	$1.BED
mv KINDEX.DAT $1.KND
mv EROSION_DEPTH_INDEX.DAT $1.DPT
mv TOPO.RAW $1TOP.RAW
mv TOPO.DAT $1.TOP
mv GLEN_FLUX.DAT $1.GFX
mv GLACIAL_FLUX_GLEN.LST $1.GLS
rename BSHADE B_$1 BSHADE????.PGM
rename SUBMRG S_$1 SUBMRG????.RAW
rename DSCHRG Q_$1 DSCHRG????.RAW 
rename RELELE E_$1 RELELE????.RAW
rename SEDFLX F_$1 SEDFLX????.RAW
mv CHANNEL.DAT	$1.CHN
mv SOURCE.DAT $1.SRC
mv RECORD.DAT $1.REC
mv REPORT.PRN $1.REP
mv DEBUG.PRN $1.DBG
mv QQ.DAT $1.QQQ
mv INITIAL_ELEVATION.DAT $1.ENN
mv EOLIAN.DAT $1.EDT
mv EWATER.DAT $1.WAT
mv LAVA.DAT $1.LAV
mv LACTIVE.DAT	$1.LCT
mv LAGE.DAT $1.LGE
mv LAGE.RAW $1LAGE.RAW
mv BASIN.LST $1.LST
mv AVALANCHE.DAT $1.AVL
mv CRATER.DAT $1.CRT
mv SUBMERGED.DAT $1.SUB
mv SEDBUDGET.PRN $1.SBG
mv CHPROP.PRN $1.CHP
mv SEDTEMP.DAT $1.SDT
mv STATISTICS.PRN $1.PRN
mv DISCHARGES.DAT $1.DSC
mv GLACIAL_FLUX_BINGHAM.LST $1.BLS
mv BINGHAM_FLUX.DAT $1.BFX
cat MARSSIM_INITIAL_BOUNDARY_CONDITIONS.PRM \
ALLUVIAL_CHANNEL_PARAMETERS.PRM \
BEDROCK_CHANNEL_PARAMETERS.PRM \
FLOW_PARAMETERS.PRM \
WEATHERING_PARAMETERS.PRM \
MASS_WASTING_PARAMETERS.PRM \
MASS_FLOW_PARAMETERS.PRM \
CRATERING_PARAMETERS.PRM \
EOLIAN_PARAMETERS.PRM \
GROUNDWATER_PARAMETERS.PRM \
ACCRETION_ABLATION_PARAMETERS.PRM \
GRAVEL_MIXTURE.PRM \
LAVA_FLOW_PARAMETERS.PRM \
 > $1.PARAMS
cp MARSSIM_INITIAL_BOUNDARY_CONDITIONS.PRM $1_MARSSIM_INITIAL_BOUNDARY_CONDITIONS.PRM
cp ALLUVIAL_CHANNEL_PARAMETERS.PRM $1_ALLUVIAL_CHANNEL_PARAMETERS.PRM
cp BEDROCK_CHANNEL_PARAMETERS.PRM $1_BEDROCK_CHANNEL_PARAMETERS.PRM
cp FLOW_PARAMETERS.PRM $1_FLOW_PARAMETERS.PRM 
cp WEATHERING_PARAMETERS.PRM $1_WEATHERING_PARAMETERS.PRM
cp MASS_WASTING_PARAMETERS.PRM $1_MASS_WASTING_PARAMETERS.PRM
cp MASS_FLOW_PARAMETERS.PRM $1_MASS_FLOW_PARAMETERS.PRM
cp CRATERING_PARAMETERS.PRM $1_CRATERING_PARAMETERS.PRM
cp EOLIAN_PARAMETERS.PRM $1_EOLIAN_PARAMETERS.PRM
cp GROUNDWATER_PARAMETERS.PRM $1_GROUNDWATER_PARAMETERS.PRM
cp ACCRETION_ABLATION_PARAMETERS.PRM $1_ACCRETION_ABLATION_PARAMETERS.PRM
cp GRAVEL_MIXTURE.PRM $1_GRAVEL_MIXTURE.PRM
cp LAVA_FLOW_PARAMETERS.PRM $1_LAVA_FLOW_PARAMETERS.PRM
cp CRATER_EVENTS.PRM $1.CVT
cp REAL_CRATERS.TXT $1.RCT
cp LAVA_SOURCES.TXT $1.LST
cp GRAVEL_VALUES.TXT $1.GVT
cp SPATIAL_VARIATION.DAT $1.SPV
mv CUMULATIVE_ELEVATION_CHANGE.DAT $1.CEC.DAT
mv CUMULATIVE_FLUVIAL_EROSION.DAT $1.CFE.DAT
mv CUMULATIVE_MASS_WASTING.DAT $1.CMW.DAT
mv CUMULATIVE_SEDIMENTATION.DAT $1.CSD.DAT
mv CUMULATIVE_WEATHERING.DAT $1.CWR.DAT
mv CUMULATIVE_EOLIAN_CHANGE.DAT $1.CEO.DAT
mv CUMULATIVE_EJECTA_DEPOSITION.DAT $1.CED.DAT
mv CUMULATIVE_CRATER_EXCAVATION.DAT $1.CCE.DAT
mv CUMULATIVE_CRATERING_CHANGE.DAT $1.CCC.DAT
mv CUMULATIVE_LAVA_CHANGE.DAT $1.CLC.DAT
mv CUMULATIVE_AVALANCHE_EROSION.DAT $1.ALV.DAT
mv CUMULATIVE_MASS_FLOW_CHANGE.DAT $1.CMF.DAT
mv TOTAL_ELEVATION_CHANGE.DAT $1.TOT.DAT
cp FINAL_SEDCOVER.DAT $1.FINSED.DAT
cp FINAL_SEDBASE.DAT $1.FINBSE.DAT
cp FINAL_BEDROCK.DAT $1.FINROCK.DAT
cp FINAL_REGOLITH.DAT $1.FINREG.DAT
cp FINAL_ELEVATION.DAT $1.FINELEV.DAT
cp FINAL_MASS_FLUX.DAT $1.FINFLUX.DAT
cp FINAL_TIME.DAT $1.FINTIME.DAT
cp INRATES.DAT $1.RTE
cp INELEV.DAT $1.EIN
cp EVENTS.PRM $1.EVT
cp BINGHAM_EVENTS.PRM $1_BINGEVENTS.BEV
mv WORK_STATISTICS.DAT $1.WRK
zip -m $1 $1.*
mv BINARY_ELEVATION.DAT $1_BINARY_ELEVATION.DAT
mv BINARY_CRATERS.DAT $1_BINARY_CRATERS.DAT
mv BINARY_SUBMERGE.DAT $1_BINARY_SUBMERGE.DAT
mv BINARY_DISCHARGES.DAT $1_BINARY_DISCHARGES.DAT
cp FINAL_SEDCOVER.DAT $1.FINAL_SEDCOVER.DAT
cp FINAL_SEDBASE.DAT $1.FINAL_SEDBASE.DAT
cp FINAL_BEDROCK.DAT $1.FINAL_BEDROCK.DAT
cp FINAL_REGOLITH.DAT $1.FINAL_REGOLITH.DAT
cp FINAL_ELEVATION.DAT $1.FINAL_ELEVATION.DAT
cp FINAL_MASS_FLUX.DAT $1.FINAL_MASS_FLUX.DAT
cp FINAL_TIME.DAT $1.FINAL_TIME.DAT