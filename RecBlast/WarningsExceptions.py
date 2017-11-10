
class RecBlastWarning(Warning):
    pass


class RecBlastException(Exception):
    pass


class RecBlastWriteError(RecBlastException):
    pass


class SearchError(RecBlastException):
    pass


class FetchSeqError(RecBlastException):
    pass


class SearchEngineNotImplementedError (SearchError, NotImplementedError):
    pass


class StopRecBlast(RecBlastException):
    pass


class IDError(RecBlastException):
    pass


class DatabaseNotFoundError(RecBlastException, FileNotFoundError):
    pass


class BLATServerError(RecBlastException):
    pass


class CheckSumError(RecBlastException):
    pass


class NoHitsError(RecBlastException):
    pass
