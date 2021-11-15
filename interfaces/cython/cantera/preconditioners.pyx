cdef class PreconditionerBase:
    """
    Common base class for reactors and reservoirs.
    """
    precon_type = "PreconditionerBase"
    def __cinit__(self, *args, **kwargs):
        self.pbase = newPreconditioner(stringify(self.precon_type))

    def __dealloc__(self):
        del self.pbase

cdef class AdaptivePreconditioner(PreconditionerBase):
    precon_type = "AdaptivePreconditioner"

    def __cinit__(self, *args, **kwargs):
        self.preconditioner = <CxxAdaptivePreconditioner*>(self.pbase)

    def getThreshold(self):
        return self.preconditioner.getThreshold()

    def setThreshold(self, val):
        self.preconditioner.setThreshold(val)

    def printContents(self):
        self.preconditioner.printPreconditioner()
