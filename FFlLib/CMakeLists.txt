cmake_minimum_required(VERSION 3.22)

project(FFlLib LANGUAGES C CXX)

set(FFA_USE_ROUND "OFF" CACHE STRING "Round to 10 digits in Nastran parser?")
set_property(CACHE FFA_USE_ROUND PROPERTY STRINGS OFF ON NEW)

ifem_add_library(
  NAME
    FFlLib
  SOURCES
    FFlAttributeBase.C
    FFlCrossSection.C
    FFlElementBase.C
    FFlFEAttributeRefs.C
    FFlFEAttributeSpec.C
    FFlFEElementTopSpec.C
    FFlFENodeRefs.C
    FFlFieldBase.C
    FFlGroup.C
    FFlLinkHandler.C
    FFlLoadBase.C
    FFlMemPool.C
    FFlNamedPartBase.C
    FFlNastranProcessor.C
    FFlNastranReader.C
    FFlNastranResolver.C
    FFlPartBase.C
    FFlUtils.C
    FFlFEParts/FFlAllFEParts.C
    FFlFEParts/FFlBeams.C
    FFlFEParts/FFlBUSH.C
    FFlFEParts/FFlCMASS.C
    FFlFEParts/FFlCurvedFace.C
    FFlFEParts/FFlLoads.C
    FFlFEParts/FFlNode.C
    FFlFEParts/FFlPBEAMECCENT.C
    FFlFEParts/FFlPBEAMPIN.C
    FFlFEParts/FFlPBEAMSECTION.C
    FFlFEParts/FFlPBUSHCOEFF.C
    FFlFEParts/FFlPBUSHECCENT.C
    FFlFEParts/FFlPCOMP.C
    FFlFEParts/FFlPCOORDSYS.C
    FFlFEParts/FFlPEFFLENGTH.C
    FFlFEParts/FFlPMASS.C
    FFlFEParts/FFlPMAT.C
    FFlFEParts/FFlPNSM.C
    FFlFEParts/FFlPORIENT.C
    FFlFEParts/FFlPRBAR.C
    FFlFEParts/FFlPRGD.C
    FFlFEParts/FFlPSPRING.C
    FFlFEParts/FFlPTHICK.C
    FFlFEParts/FFlPWAVGM.C
    FFlFEParts/FFlRBAR.C
    FFlFEParts/FFlRGD.C
    FFlFEParts/FFlShells.C
    FFlFEParts/FFlSolids.C
    FFlFEParts/FFlSPRING.C
    FFlFEParts/FFlWAVGM.C
    FFaLib/FFaAlgebra/FFaAlgebra.C
    FFaLib/FFaAlgebra/FFaCheckSum.C
    FFaLib/FFaAlgebra/FFaMat33.C
    FFaLib/FFaAlgebra/FFaMat34.C
    FFaLib/FFaAlgebra/FFaMath.C
    FFaLib/FFaAlgebra/FFaTensor1.C
    FFaLib/FFaAlgebra/FFaTensor2.C
    FFaLib/FFaAlgebra/FFaTensor3.C
    FFaLib/FFaAlgebra/FFaTensorTransforms.C
    FFaLib/FFaAlgebra/FFaTensorTransforms_F.C
    FFaLib/FFaAlgebra/FFaVec3.C
    FFaLib/FFaAlgebra/FFaVolume.C
    FFaLib/FFaDefinitions/FFaMsg.C
    FFaLib/FFaOS/FFaFilePath.C
    FFaLib/FFaPatterns/FFaMemPool.C
  HEADERS
    FFlAttributeBase.H
    FFlCrossSection.H
    FFlElementBase.H
    FFlFEAttributeRefs.H
    FFlFEAttributeSpec.H
    FFlFEElementTopSpec.H
    FFlFENodeRefs.H
    FFlFEResultBase.H
    FFlFieldBase.H
    FFlField.H
    FFlGroup.H
    FFlLinkCSMask.H
    FFlLinkHandler.H
    FFlLoadBase.H
    FFlMemPool.H
    FFlNamedPartBase.H
    FFlNastranReader.H
    FFlPartBase.H
    FFlReaderBase.H
    FFlReaders.H
    FFlReference.H
    FFlTypeInfoSpec.H
    FFlUtils.H
    FFlFEParts/FFlAllFEParts.H
    FFlFEParts/FFlBeams.H
    FFlFEParts/FFlBUSH.H
    FFlFEParts/FFlCMASS.H
    FFlFEParts/FFlCurvedFace.H
    FFlFEParts/FFlLoads.H
    FFlFEParts/FFlNode.H
    FFlFEParts/FFlPBEAMECCENT.H
    FFlFEParts/FFlPBEAMPIN.H
    FFlFEParts/FFlPBEAMSECTION.H
    FFlFEParts/FFlPBUSHCOEFF.H
    FFlFEParts/FFlPBUSHECCENT.H
    FFlFEParts/FFlPCOMP.H
    FFlFEParts/FFlPCOORDSYS.H
    FFlFEParts/FFlPEFFLENGTH.H
    FFlFEParts/FFlPMASS.H
    FFlFEParts/FFlPMAT.H
    FFlFEParts/FFlPNSM.H
    FFlFEParts/FFlPORIENT.H
    FFlFEParts/FFlPRBAR.H
    FFlFEParts/FFlPRGD.H
    FFlFEParts/FFlPSPRING.H
    FFlFEParts/FFlPTHICK.H
    FFlFEParts/FFlPWAVGM.H
    FFlFEParts/FFlRBAR.H
    FFlFEParts/FFlRGD.H
    FFlFEParts/FFlShells.H
    FFlFEParts/FFlSolids.H
    FFlFEParts/FFlSPRING.H
    FFlFEParts/FFlWAVGM.H
    FFaLib/FFaAlgebra/FFaAlgebra.H
    FFaLib/FFaAlgebra/FFaCheckSum.H
    FFaLib/FFaAlgebra/FFaMat33.H
    FFaLib/FFaAlgebra/FFaMat34.H
    FFaLib/FFaAlgebra/FFaMath.H
    FFaLib/FFaAlgebra/FFaTensor1.H
    FFaLib/FFaAlgebra/FFaTensor2.H
    FFaLib/FFaAlgebra/FFaTensor3.H
    FFaLib/FFaAlgebra/FFaTensorTransforms.H
    FFaLib/FFaAlgebra/FFaUnitCalculator.H
    FFaLib/FFaAlgebra/FFaVec3.H
    FFaLib/FFaAlgebra/FFaVolume.H
    FFaLib/FFaDefinitions/FFaMsg.H
    FFaLib/FFaDynCalls/FFaDynCB.H
    FFaLib/FFaOS/FFaFilePath.H
    FFaLib/FFaOS/FFaFortran.H
    FFaLib/FFaPatterns/FFaGenericFactory.H
    FFaLib/FFaPatterns/FFaMemPool.H
    FFaLib/FFaPatterns/FFaSingelton.H
)

target_compile_definitions(FFlLib PUBLIC $<$<COMPILE_LANGUAGE:CXX>:HAS_FFLLIB=1 FFL_REFCOUNT=1>)
target_compile_definitions(FFlLib PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:HAS_FFLLIB=1>)

if(FFL_USE_MEMPOOL)
  target_compile_definitions(FFlLib PUBLIC FT_USE_MEMPOOL=1)
endif()
if(IFEM_USE_ANDES)
  target_compile_definitions(FFlLib PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:HAS_ANDES=1>)
endif()

if(VERBOSE_DEBUG GREATER 12)
  math(EXPR FFL_DEBUG "${VERBOSE_DEBUG} - 10")
elseif(VERBOSE_DEBUG GREATER 1)
  set(FFL_DEBUG 2)
else()
  set(FFL_DEBUG 0)
endif()

if(FFL_DEBUG GREATER 0)
  target_compile_options(FFlLib PRIVATE $<$<CONFIG:Debug>:FFA_DEBUG=1 FFL_DEBUG=${FFL_DEBUG}>)
endif()

if("${FFA_USE_ROUND}" MATCHES "NEW")
  target_compile_definitions(FFlLib PRIVATE FFA_USE_ROUND=1)
elseif("${FFA_USE_ROUND}" MATCHES "OFF")
  target_compile_definitions(FFlLib PRIVATE FFA_NO_ROUND=1)
endif()
