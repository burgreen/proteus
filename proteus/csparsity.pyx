# A type of -*- python -*- file
import numpy as np
cimport numpy as np
from csparsity cimport SparsityInfo

cdef class PySparsityInfo:
    cdef SparsityInfo cpp
    def __cinit__(self):
        pass
    def __dealloc__(self):
        pass
        #self.cpp.columnIndecesMap.clear()
        #self.cpp.columnOffsetsMap.clear()
    def findNonzeros(self,
                     int nElements_global,
                     int nDOF_test_element,
                     int nDOF_trial_element,
                     np.ndarray nFreeDOF_test,
                     np.ndarray freeGlobal_test,
                     np.ndarray nFreeDOF_trial,
                     np.ndarray freeGlobal_trial,
                     int offset_test,
                     int stride_test,
                     int offset_trial,
                     int stride_trial,
                     int hasNumericalFlux,
                     int hasDiffusionInMixedForm,
                     int needNumericalFluxJacobian,
                     int nElementBoundaries_element,
                     np.ndarray elementNeighborsArray,
                     int nInteriorElementBoundaries_global,
                     np.ndarray interiorElementBoundariesArray,
                     np.ndarray elementBoundaryElementsArray,
                     np.ndarray elementBoundaryLocalElementBoundariesArray,
                     int hasFluxBoundaryConditions,
                     int nExteriorElementBoundaries_global,
                     np.ndarray exteriorElementBoundariesArray,
                     int hasOutflowBoundary,
                     int needOutflowJacobian):
        self.cpp.findNonzeros(nElements_global,
                              nDOF_test_element,
                              nDOF_trial_element,
                              <int*> nFreeDOF_test.data,
                              <int*> freeGlobal_test.data,
                              <int*> nFreeDOF_trial.data,
                              <int*> freeGlobal_trial.data,
                              offset_test,
                              stride_test,
                              offset_trial,
                              stride_trial,
                              hasNumericalFlux,
                              hasDiffusionInMixedForm,
                              needNumericalFluxJacobian,
                              nElementBoundaries_element,
                              <int*> elementNeighborsArray.data,
                              nInteriorElementBoundaries_global,
                              <int*> interiorElementBoundariesArray.data,
                              <int*> elementBoundaryElementsArray.data,
                              <int*> elementBoundaryLocalElementBoundariesArray.data,
                              hasFluxBoundaryConditions,
                              nExteriorElementBoundaries_global,
                              <int*> exteriorElementBoundariesArray.data,
                              hasOutflowBoundary,
                              needOutflowJacobian)
    def getOffsets_CSR(self,
                       int nElements_global,
                       int nDOF_test_element,
                       int nDOF_trial_element,
                       np.ndarray nFreeDOF_test,
                       np.ndarray freeGlobal_test,
                       np.ndarray nFreeDOF_trial,
                       np.ndarray freeGlobal_trial,
                       int offset_test,
                       int stride_test,
                       int offset_trial,
                       int stride_trial,
                       int hasNumericalFlux,
                       int hasDiffusionInMixedForm,
                       int needNumericalFluxJacobian,
                       int nElementBoundaries_element,
                       np.ndarray elementNeighborsArray,
                       int nInteriorElementBoundaries_global,
                       np.ndarray interiorElementBoundariesArray,
                       np.ndarray elementBoundaryElementsArray,
                       np.ndarray elementBoundaryLocalElementBoundariesArray,
                       int hasFluxBoundaryConditions,
                       int nExteriorElementBoundaries_global,
                       np.ndarray exteriorElementBoundariesArray,
                       int hasOutflowBoundary,
                       int needOutflowJacobian,
                       np.ndarray rowptr,
                       np.ndarray csrRowIndeces,
                       np.ndarray csrColumnOffsets,
                       np.ndarray csrColumnOffsets_eNebN,
                       np.ndarray csrColumnOffsets_eb,
                       np.ndarray csrColumnOffsets_eb_eNebN):
        self.cpp.getOffsets_CSR(nElements_global,
                                nDOF_test_element,
                                nDOF_trial_element,
                                <int*> nFreeDOF_test.data,
                                <int*> freeGlobal_test.data,
                                <int*> nFreeDOF_trial.data,
                                <int*> freeGlobal_trial.data,
                                offset_test,
                                stride_test,
                                offset_trial,
                                stride_trial,
                                hasNumericalFlux,
                                hasDiffusionInMixedForm,
                                needNumericalFluxJacobian,
                                nElementBoundaries_element,
                                <int*> elementNeighborsArray.data,
                                nInteriorElementBoundaries_global,
                                <int*> interiorElementBoundariesArray.data,
                                <int*> elementBoundaryElementsArray.data,
                                <int*> elementBoundaryLocalElementBoundariesArray.data,
                                hasFluxBoundaryConditions,
                                nExteriorElementBoundaries_global,
                                <int*> exteriorElementBoundariesArray.data,
                                hasOutflowBoundary,
                                needOutflowJacobian,
                                <int*> rowptr.data,
                                <int*> csrRowIndeces.data,
                                <int*> csrColumnOffsets.data,
                                <int*> csrColumnOffsets_eNebN.data,
                                <int*> csrColumnOffsets_eb.data,
                                <int*> csrColumnOffsets_eb_eNebN.data)
    def getCSR(self):
        self.cpp.getCSR()
        return (np.asarray(<int[:self.cpp.nrows]>self.cpp.rowptr),
                np.asarray(<int[:self.cpp.nnz]>self.cpp.colind),
                self.cpp.nnz,
                np.asarray(<double[:self.cpp.nnz]>self.cpp.nzval))