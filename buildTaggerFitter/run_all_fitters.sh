#!/bin/bash

condor_submit condorConfig_maxHits_auc.sub
condor_submit condorConfig_maxHits_maxEff.sub
condor_submit condorConfig_maxHits_sigmaInit.sub
condor_submit condorConfig_maxHits_sigmaL0.sub
condor_submit condorConfig_maxHits_sigmaL1.sub

condor_submit condorConfig_maxLayers_auc.sub
condor_submit condorConfig_maxLayers_maxEff.sub
condor_submit condorConfig_maxLayers_sigmaInit.sub
condor_submit condorConfig_maxLayers_sigmaL0.sub
condor_submit condorConfig_maxLayers_sigmaL1.sub

condor_submit condorConfig_maxLength_auc.sub
condor_submit condorConfig_maxLength_maxEff.sub
condor_submit condorConfig_maxLength_sigmaInit.sub
condor_submit condorConfig_maxLength_sigmaL0.sub
condor_submit condorConfig_maxLength_sigmaL1.sub

condor_submit condorConfig_nHelices_auc.sub
condor_submit condorConfig_nHelices_maxEff.sub
condor_submit condorConfig_nHelices_sigmaInit.sub
condor_submit condorConfig_nHelices_sigmaL0.sub
condor_submit condorConfig_nHelices_sigmaL1.sub