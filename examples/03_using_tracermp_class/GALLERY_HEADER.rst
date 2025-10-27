Using the TracerMP Class
========================

The :class:`~mapflpy.tracer.TracerMP` class is a multi-processing, thread-safe variant of the
:class:`~mapflpy.tracer.Tracer` class. It enables users to create multiple (distinct) instances of the
:mod:`mapflpy_fortran` object in separate processes, allowing for concurrent tracing operations
across multiple threads. the :class:`~mapflpy.tracer.TracerMP` branches new instances of the
:mod:`mapflpy_fortran` object into discrete processes and communicates with these instances through
python's :py:func:`~multiprocessing.Pipe` protocol.

The :py:func:`~multiprocessing.Pipe` protocol requires that
the :py:class:`multiprocessing.connection.Connection` is properly opened and closed to
avoid resource leaks. Therefore, it is recommended to use the :class:`~mapflpy.tracer.TracerMP`
class within a context manager (the ``with`` statement) to ensure that resources are
properly managed (or, alternatively, one can explicitly open and close the connection using the
:meth:`~mapflpy.tracer.TracerMP.connect` and :meth:`~mapflpy.tracer.TracerMP.disconnect` methods).