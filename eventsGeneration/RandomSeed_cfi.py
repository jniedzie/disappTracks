import FWCore.ParameterSet.Config as cms

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper


def customizeRandomSeed(process):
  
  randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
  randSvc.populate()

  return(process)
