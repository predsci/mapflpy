Using the Tracer Class
======================

If you need more control over the tracing process, you can use the :class:`~mapflpy.tracer.Tracer`
class directly. As mentioned throughout these examples, the :class:`~mapflpy.tracer.Tracer` class
enforces a singleton pattern to manage issues that arise from the underlying :mod:`mapflpy_fortran` object
not being thread-safe. As a result, it is recommended to use the :class:`~mapflpy.tracer.Tracer` class
in single-threaded contexts only *viz.* instantiating one instance of the class at a time.

The benefit of using the :class:`~mapflpy.tracer.Tracer` class directly is that it provides direct
access to the :mod:`mapflpy_fortran` object – allowing for faster and more flexible tracing operations.
Alternatively, the :class:`~mapflpy.tracer.TracerMP` (while thread-safe) branches new instances of the
:mod:`mapflpy_fortran` object into discrete processes and communicates with these instances through
python's :py:func:`~multiprocessing.Pipe` protocol. Although this approach is more robust in multi-threaded
contexts, it incurs a performance penalty due to the overhead of inter-process communication.
